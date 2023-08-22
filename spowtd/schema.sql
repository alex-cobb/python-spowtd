-- Data schema for Spowtd

PRAGMA foreign_keys = 1;


-- Staging tables, prior to regridding
CREATE TABLE rainfall_intensity_staging (
  epoch integer NOT NULL PRIMARY KEY,
  rainfall_intensity_mm_h double precision NOT NULL
);


CREATE TABLE water_level_staging (
  epoch integer NOT NULL PRIMARY KEY,
  zeta_mm double precision NOT NULL
);


CREATE TABLE evapotranspiration_staging (
  epoch integer NOT NULL PRIMARY KEY,
  evapotranspiration_mm_h double precision NOT NULL
);

-- Gridded tables

CREATE TABLE time_grid (
  time_step_s integer NOT NULL,
  source_time_zone text NOT NULL,
  -- Singleton
  is_valid integer NOT NULL PRIMARY KEY
    CHECK (is_valid = 1)
    DEFAULT 1
);


CREATE TABLE grid_time (
  epoch integer NOT NULL PRIMARY KEY,
  data_interval integer NULL
);


CREATE TABLE thresholds (
  storm_rain_threshold_mm_h double precision NOT NULL,
  rising_jump_threshold_mm_h double precision NOT NULL,
  -- Singleton
  is_valid integer NOT NULL PRIMARY KEY
    CHECK (is_valid = 1)
    DEFAULT 1
);


CREATE TABLE grid_time_flags (
  start_epoch integer NOT NULL PRIMARY KEY
    REFERENCES grid_time(epoch),
  is_jump boolean NOT NULL,
  is_mystery_jump boolean NOT NULL,
  is_interstorm boolean NOT NULL
);


CREATE TABLE rainfall_intensity (
  from_epoch integer NOT NULL
    REFERENCES grid_time (epoch)
    CHECK (from_epoch < thru_epoch),
  thru_epoch integer NOT NULL
    REFERENCES grid_time (epoch)
    CHECK (from_epoch < thru_epoch),
  rainfall_intensity_mm_h double precision NOT NULL,
  PRIMARY KEY (from_epoch)
);


CREATE TABLE evapotranspiration (
  from_epoch integer NOT NULL
    REFERENCES grid_time (epoch)
    CHECK (from_epoch < thru_epoch),
  thru_epoch integer NOT NULL
    REFERENCES grid_time (epoch)
    CHECK (from_epoch < thru_epoch),
  evapotranspiration_mm_h double precision NOT NULL,
  PRIMARY KEY (from_epoch)
);


CREATE TABLE water_level (
  epoch integer NOT NULL
    REFERENCES grid_time (epoch),
  zeta_mm double precision NOT NULL,
  PRIMARY KEY (epoch)
);


CREATE TABLE storm (
  start_epoch integer NOT NULL
    REFERENCES grid_time (epoch)
    CHECK (start_epoch < thru_epoch),
  thru_epoch integer NOT NULL
    REFERENCES grid_time (epoch)
    CHECK (start_epoch < thru_epoch),
  PRIMARY KEY (start_epoch)
);


CREATE TABLE zeta_interval (
  start_epoch integer NOT NULL PRIMARY KEY
    REFERENCES water_level (epoch),
  interval_type text NOT NULL
    CHECK (interval_type in ('storm', 'interstorm')),
  thru_epoch integer NOT NULL
    CHECK (start_epoch < thru_epoch)
    REFERENCES water_level (epoch),
  -- For FK
  UNIQUE (start_epoch, interval_type)
);


CREATE TABLE zeta_interval_storm (
  interval_start_epoch integer NOT NULL PRIMARY KEY,
  interval_type text NOT NULL CHECK (interval_type = 'storm'),
  storm_start_epoch integer NOT NULL UNIQUE,
  FOREIGN KEY (interval_start_epoch, interval_type)
    REFERENCES zeta_interval (start_epoch,
                              interval_type),
  FOREIGN KEY (storm_start_epoch)
    REFERENCES storm (start_epoch)
);


CREATE TABLE zeta_grid (
  -- Enforce singleton
  id boolean PRIMARY KEY DEFAULT TRUE CHECK (id = TRUE),
  grid_interval_mm double precision NOT NULL
);


CREATE TABLE discrete_zeta (
  zeta_number integer PRIMARY KEY,
  zeta_grid boolean NOT NULL DEFAULT TRUE REFERENCES zeta_grid(id)
);


CREATE TABLE rising_interval (
  start_epoch integer NOT NULL PRIMARY KEY
    REFERENCES zeta_interval_storm (interval_start_epoch),
  interval_type text NOT NULL DEFAULT 'storm'
    CHECK (interval_type = 'storm'),
  rain_depth_offset_mm double precision NOT NULL
);


CREATE TABLE recession_interval (
  start_epoch integer NOT NULL PRIMARY KEY,
  interval_type text NOT NULL DEFAULT 'interstorm'
    CHECK (interval_type = 'interstorm'),
  time_offset_s double precision NOT NULL,
  FOREIGN KEY (start_epoch, interval_type)
    REFERENCES zeta_interval (start_epoch, interval_type)
);


CREATE TABLE rising_interval_zeta (
  start_epoch integer NOT NULL
    REFERENCES rising_interval (start_epoch),
  zeta_number integer NOT NULL REFERENCES discrete_zeta(zeta_number),
  mean_crossing_depth_mm double precision NOT NULL,
  PRIMARY KEY (start_epoch, zeta_number)
);


CREATE TABLE recession_interval_zeta (
  start_epoch integer NOT NULL
    REFERENCES recession_interval (start_epoch),
  zeta_number integer NOT NULL REFERENCES discrete_zeta(zeta_number),
  mean_crossing_time interval NOT NULL,
  PRIMARY KEY (start_epoch, zeta_number)
);


CREATE TABLE curvature (
  curvature_m_km2 double precision NOT NULL,
  -- Singleton
  is_valid integer NOT NULL PRIMARY KEY
    CHECK (is_valid = 1)
    DEFAULT 1
);


-- Views

CREATE VIEW storm_total_rain_depth AS
SELECT s.start_epoch AS storm_start_epoch,
       SUM(ri.rainfall_intensity_mm_h *
           (ri.thru_epoch - ri.from_epoch)
	   / 3600.
           ) AS total_depth_mm
FROM storm AS s
JOIN rainfall_intensity AS ri
  ON ri.from_epoch >= s.start_epoch
  AND ri.thru_epoch <= s.thru_epoch
GROUP BY s.start_epoch;


CREATE VIEW average_recession_time AS
SELECT zeta_number * zg.grid_interval_mm AS zeta_mm,
       AVG(time_offset_s + mean_crossing_time)
	 AS elapsed_time_s
FROM recession_interval AS ri
JOIN recession_interval_zeta
  USING (start_epoch)
JOIN discrete_zeta AS dz
  USING (zeta_number)
JOIN zeta_grid AS zg
  ON zg.id = dz.zeta_grid
GROUP BY zeta_number, grid_interval_mm;


CREATE VIEW average_rising_depth AS
SELECT zeta_number * zg.grid_interval_mm AS zeta_mm,
       AVG(rain_depth_offset_mm + mean_crossing_depth_mm)
	 AS mean_crossing_depth_mm
FROM rising_interval AS ri
JOIN rising_interval_zeta
  USING (start_epoch)
JOIN discrete_zeta AS dz
  USING (zeta_number)
JOIN zeta_grid AS zg
  ON zg.id = dz.zeta_grid
GROUP BY zeta_number, grid_interval_mm;


CREATE VIEW storm_total_rise AS
SELECT zis.storm_start_epoch,
       zis.interval_start_epoch,
       iz.zeta_mm AS initial_zeta_mm,
       fz.zeta_mm AS final_zeta_mm
FROM zeta_interval_storm AS zis
JOIN zeta_interval AS zi
  ON zis.interval_start_epoch = zi.start_epoch
JOIN water_level AS iz
  ON iz.epoch = zi.start_epoch
JOIN water_level AS fz
  ON fz.epoch = zi.thru_epoch;


CREATE VIEW rising_curve_line_segment AS
SELECT interval_start_epoch,
       ri.rain_depth_offset_mm,
       total_depth_mm AS rain_total_depth_mm,
       initial_zeta_mm,
       final_zeta_mm,
       (ABS(final_zeta_mm - initial_zeta_mm) / 100.) AS factor_stretch,
       (total_depth_mm/(ABS(final_zeta_mm - initial_zeta_mm) / 100.)) AS stretched_rain_total_depth_mm,
       (((initial_zeta_mm+final_zeta_mm)/ 2.) - 50.0) AS stretched_initial_zeta_mm,
       (((initial_zeta_mm+final_zeta_mm)/ 2.) + 50.0) AS stretched_final_zeta_mm
FROM storm_total_rise AS str
JOIN storm_total_rain_depth AS strd
  USING (storm_start_epoch)
JOIN rising_interval AS ri
  ON ri.start_epoch = str.interval_start_epoch;
-- SA! additional variables are created for stretching and these have to be manually adjusted, if still included
