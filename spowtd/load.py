"""Load data into Spowtd data file

"""

import csv as csv_mod
import datetime as datetime_mod
import logging
import os

import numpy as np
import pytz


ISO_8601_FORMAT = '%Y-%m-%d %H:%M:%S'
SCHEMA_PATH = os.path.join(
    os.path.dirname(__file__),
    'schema.sql')
LOG = logging.getLogger('spowtd.load')


def load_data(connection,
              precipitation_data_file,
              evapotranspiration_data_file,
              water_level_data_file,
              time_zone_name='Africa/Lagos'):
    """Load data into Spowtd data file

    """
    time_zone = pytz.timezone(time_zone_name)

    connection.execute("PRAGMA foreign_keys = 1")
    cursor = connection.cursor()
    with open(SCHEMA_PATH, 'rt') as schema_file:
        cursor.executescript(schema_file.read())
    # Format: Datetime, precipitation intensity (mm / h)
    precip_csv = csv_mod.reader(precipitation_data_file,
                                delimiter=',')
    header = next(precip_csv)
    assert header[0].lower().startswith('datetime'), header
    # Strategy: load into staging tables first, then establish time
    # grid
    cursor.executemany("""
    INSERT INTO rainfall_intensity_staging
      (epoch, rainfall_intensity_mm_h)
    VALUES (?, ?)""", generate_timestamped_rows(precip_csv,
                                                time_zone))
    del precip_csv

    # Format: Datetime, evapotranspiration (mm / h)
    et_csv = csv_mod.reader(evapotranspiration_data_file,
                            delimiter=',')
    header = next(et_csv)
    assert header[0].lower().startswith('datetime'), header
    cursor.executemany("""
    INSERT INTO evapotranspiration_staging
      (epoch, evapotranspiration_mm_h)
    VALUES (?, ?)""", generate_timestamped_rows(et_csv,
                                                time_zone))
    del et_csv

    water_level_csv = csv_mod.reader(water_level_data_file,
                                     delimiter=',')
    # Format: Datetime, water level (mm)
    header = next(water_level_csv)
    assert header[0].lower().startswith('datetime'), header
    cursor.executemany("""
    INSERT INTO water_level_staging
      (epoch, zeta_mm)
    VALUES (?, ?)""", generate_timestamped_rows(water_level_csv,
                                                time_zone))
    time_grid, time_step = populate_grid_time(cursor)
    populate_rainfall_intensity(cursor, time_grid, time_step)
    populate_water_level(cursor, time_grid, time_step)
    cursor.close()
    connection.commit()


def populate_water_level(cursor, time_grid, time_step):
    """Interpolate water level onto precipitation time grid

    """
    cursor.execute("""
    SELECT epoch, zeta_mm
    FROM water_level_staging""")
    zeta_t, zeta_mm = zip(*cursor.fetchall())
    zeta_on_grid = np.interp(time_grid[:-1], zeta_t, zeta_mm)
    cursor.executemany("""
    INSERT INTO water_level (epoch, zeta_mm)
    VALUES (?, ?)""", zip(time_grid[:-1], zeta_on_grid))


def populate_grid_time(cursor):
    """Determine and populate grid_time

    Identifies interval with both precipitation and water level data,
    and returns the time grid array and time step as a tuple.

    """
    time_grid = [
        epoch for epoch, in cursor.execute("""
        WITH a AS (
          SELECT min(epoch) AS min_t_zeta,
                 max(epoch) AS max_t_zeta
          FROM water_level_staging
        )
        SELECT epoch
        FROM rainfall_intensity_staging AS ris
        JOIN a
          ON ris.epoch >= min_t_zeta
          AND ris.epoch <= max_t_zeta
        ORDER BY epoch""")]
    delta_t = sorted(set(np.diff(time_grid)))
    if len(delta_t) != 1:
        raise ValueError(
            'Nonuniform time steps in rainfall data: {} s'
            .format(delta_t))
    time_step = int(delta_t[0])
    del delta_t
    cursor.execute("""
    INSERT INTO time_grid (time_step_s)
    VALUES (?)""", (time_step,))
    # Add a grid time for the end of the last rainfall interval
    time_grid.append(
        time_grid[-1] + time_step)
    cursor.executemany("""
    INSERT INTO grid_time (epoch)
    VALUES (?)""", [(epoch,) for epoch in time_grid])
    return (time_grid, time_step)


def populate_rainfall_intensity(cursor, time_grid, time_step):
    """Populate rainfall intensity on target grid

    """
    cursor.execute("""
    INSERT INTO rainfall_intensity
      (from_epoch, thru_epoch, rainfall_intensity_mm_h)
    SELECT ris.epoch, ris.epoch + ?, rainfall_intensity_mm_h
    FROM rainfall_intensity_staging AS ris
    JOIN grid_time AS gt
      USING (epoch)
    WHERE ris.epoch <= ?""", (time_step,
                              time_grid[-2]))


def generate_timestamped_rows(rows, tz):
    """Generate rows with the first value replaced by a UNIX timestamp

    The first item in each row is assumed to be a text datetime
    in ISO 8601 format.

    The pytz object passed as the second argument is used to convert
    the datetime to UTC (if necessary) and convert to a UNIX timestamp
    (seconds since 1970-01-01 00:00:00).

    """
    for row in rows:
        local_datetime = tz.localize(
            datetime_mod.datetime.strptime(
                row[0], ISO_8601_FORMAT))
        is_aware = (local_datetime.tzinfo.utcoffset(local_datetime)
                    is not None)
        assert is_aware
        epoch = local_datetime.timestamp()
        if not epoch.is_integer():
            raise ValueError('Non-integer seconds in datetime {}'
                             .format(row[0]))
        yield [int(epoch)] + row[1:]
