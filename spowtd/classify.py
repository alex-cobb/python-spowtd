"""Classify data into storm and interstorm intervals

Match intervals of rapidly increasing water level (“rises”) to intervals of
heavy rain (“storms”) in such a way that each storm is matched to no more than
one rise and each rise is matched to no more than one storm. This matching is
performed in two steps. First, all storms and rises that overlap in time are
matched. This first step may result in matching from a single storm to multiple
rises and vice versa. This step is followed by an arbitration step based on a
variant of the Gale-Shapley deferred acceptance algorithm for the stable
matching problem: it finds a set of matches between storms and rises that is
stable in the sense that, by switching a pair of matches between storms and
rises, one cannot improve the agreement in duration and start time for both
matches.

"""

import datetime
import logging

from collections import defaultdict

import numpy as np
import pandas as pd

LOG = logging.getLogger("spowtd.classify")


def classify_intervals(
    connection, storm_rain_threshold_mm_h=4.0, rising_jump_threshold_mm_h=8.0
):
    """Classify data into storm and interstorm intervals"""
    cursor = connection.cursor()
    cursor.execute(
        """
    INSERT INTO thresholds
      (storm_rain_threshold_mm_h, rising_jump_threshold_mm_h)
    VALUES
      (:storm_rain_threshold_mm_h, :rising_jump_threshold_mm_h)""",
        {
            "storm_rain_threshold_mm_h": storm_rain_threshold_mm_h,
            "rising_jump_threshold_mm_h": rising_jump_threshold_mm_h,
        },
    )
    cursor.execute(
        """
    SELECT DISTINCT data_interval
    FROM grid_time
    WHERE data_interval IS NOT NULL
    ORDER BY data_interval"""
    )
    data_intervals = [row[0] for row in cursor.fetchall()]
    if not data_intervals:
        raise ValueError("No valid data intervals found")
    for data_interval in data_intervals:
        populate_zeta_interval(
            cursor,
            data_interval,
            storm_rain_threshold_mm_h,
            rising_jump_threshold_mm_h,
        )
    cursor.close()
    connection.commit()


def populate_zeta_interval(
    cursor,
    data_interval,
    storm_rain_threshold_mm_h,
    rising_jump_threshold_mm_h,
):
    """Identify storm and interstorm intervals"""
    classify_interstorms(cursor, data_interval, rising_jump_threshold_mm_h)
    match_all_storms(
        cursor,
        data_interval,
        storm_rain_threshold_mm_h,
        rising_jump_threshold_mm_h,
    )

# SA! Rainfall intensity threshold set to 0.5 mm/h to cope with IMERG drizzling
# SA! Still manual adjustment needed when used with in situ precipitation/throughfall
def classify_interstorms(cursor, data_interval, rising_jump_threshold_mm_h):
    """Populate interstorm intervals"""
    (epoch, zeta_mm, is_raining) = (
        np.array(v)
        for v in zip(
            *cursor.execute(
                """
         SELECT water_level.epoch,
                zeta_mm,
                rainfall_intensity_mm_h > 0.5 AS is_raining
         FROM grid_time
         JOIN rainfall_intensity
           ON rainfall_intensity.from_epoch = grid_time.epoch
           AND grid_time.data_interval = ?
         JOIN water_level
           ON rainfall_intensity.from_epoch = water_level.epoch
         ORDER BY from_epoch""",
                (data_interval,),
            )
        )
    )
    assert len(epoch), epoch.shape
    check_for_uniform_time_steps(epoch)
    hour = epoch / 3600.0
    is_raining = is_raining.astype(bool)
    assert np.isfinite(zeta_mm).all()
    # Look for jumps in head much bigger than noise, which could
    # indicate the onset of rain, and mark everything after the jump
    # until the next rain as a "mystery jump".
    rates = np.concatenate(([0], (zeta_mm[1:] - zeta_mm[:-1]) / (hour[1:] - hour[:-1])))
    is_jump = (rates > rising_jump_threshold_mm_h).astype(bool)
    is_mystery_jump = get_mystery_jump_mask(is_jump, is_raining, zeta_mm) #SA! added zeta_mm to stop mystery_jumps sooner
    is_interstorm = (~is_mystery_jump) & (~is_raining)
    interval_mask = is_interstorm
    del is_raining

    cursor.executemany(
        """
    INSERT INTO grid_time_flags
      (start_epoch, is_jump, is_mystery_jump, is_interstorm)
    VALUES
      (?, ?, ?, ?)""",
        zip(
            (int(t) for t in epoch),
            (int(b) for b in is_jump),
            (int(b) for b in is_mystery_jump),
            (int(b) for b in is_interstorm),
        ),
    )
    del is_jump, is_mystery_jump, is_interstorm

    masks = get_true_interval_masks(interval_mask)

    mask = None
    series_indices = [np.nonzero(mask)[0] for mask in masks]
    del mask, masks
    indices = None
    series_indices = [indices for indices in series_indices if len(indices) > 1]
    del indices

    LOG.info("%s series found", len(series_indices))

    #SA! Only inlcude recession events that last at least 24 hours - include diurnal effect on recession rate
    #SA! except if the initial zeta_mm value of the recession event belongs to the highest 10% of WL, then include without cropping
    #SA! Round/crop recession events to full days to even night and day time recession contribution
    duration_recession = 1

    if duration_recession == 1:
        for indices in series_indices:
            if (epoch[indices[-1]] - epoch[indices[0]] > 86400) | (zeta_mm[indices[0]] >= np.quantile(zeta_mm, 0.90)): # full day 86400  20 hours 72000
                remainder = (epoch[indices[-1]] - epoch[indices[0]]) % 86400
                if (epoch[indices[-1]] - epoch[indices[0]] > 86400):
                    epoch[indices[-1]] = epoch[indices[-1]] - remainder
                cursor.execute(
                    """
                    INSERT INTO zeta_interval
                        (start_epoch, interval_type, thru_epoch)
                    SELECT
                        :start_epoch, :interval_type, :thru_epoch
                    """,
                    {
                        "interval_type": "interstorm",
                        "start_epoch": int(epoch[indices[0]]),
                        "thru_epoch": int(epoch[indices[-1]]),
                    },
                )
    else:
        for indices in series_indices:
		cursor.execute(
                """
                INSERT INTO zeta_interval
                    (start_epoch, interval_type, thru_epoch)
                SELECT
                    :start_epoch, :interval_type, :thru_epoch
                """,
                {
                    "interval_type": "interstorm",
                    "start_epoch": int(epoch[indices[0]]),
                    "thru_epoch": int(epoch[indices[-1]]),
                },
            )
    del hour, zeta_mm
    del series_indices


def match_all_storms(
    cursor,
    data_interval,
    storm_rain_threshold_mm_h,
    rising_jump_threshold_mm_h,
):
    """Match intervals of increasing head with rainstorms

    Calls match_storms on rain data and head data from the database.

    Thresholds are given in mm / h.

    """
    (epoch, zeta_mm, rainfall_intensity_mm_h) = (
        np.array(v)
        for v in zip(
            *cursor.execute(
                """
         SELECT water_level.epoch,
                zeta_mm,
                rainfall_intensity_mm_h
         FROM grid_time
         JOIN rainfall_intensity
           ON rainfall_intensity.from_epoch = grid_time.epoch
           AND grid_time.data_interval = ?
         JOIN water_level
           ON rainfall_intensity.from_epoch = water_level.epoch
         ORDER BY from_epoch""",
                (data_interval,),
            )
        )
    )
    check_for_uniform_time_steps(epoch)
    (time_step_h,) = cursor.execute(
        """
    SELECT CAST(time_step_s AS double precision) / 3600.
    FROM time_grid"""
    ).fetchone()
    jump_delta_threshold = rising_jump_threshold_mm_h * time_step_h
    (rain_intervals, jump_intervals) = match_storms(
        rainfall_intensity_mm_h,
        zeta_mm,
        storm_rain_threshold_mm_h,
        jump_delta_threshold,
    )

    is_storm = rainfall_intensity_mm_h > storm_rain_threshold_mm_h
    for i, (rain_start, rain_stop) in enumerate(rain_intervals):
        assert is_storm[rain_start:rain_stop].all(), (
            "Storm interval [{}, {}] includes only raining time steps"
            # XXX Convert to datetime for error message?
            .format(epoch[rain_start], epoch[rain_stop - 1])
        )
        jump_start, jump_stop = jump_intervals[i]
        assert (np.diff(zeta_mm[jump_start:jump_stop]) > jump_delta_threshold).all()

        # Times associated with rainfall intensities are *start* times
        # for the interval, so the time slice that *starts* at the
        # storm thru_epoch is not included; but the slice that goes
        # thru the thru_epoch is.
        storm_start_epoch = int(epoch[rain_start])
        storm_thru_epoch = int(epoch[rain_stop])
        # On the other hand, heads are instantaneous values, so the
        # epoch of the end of the jump interval is the one to use.
        jump_start_epoch = int(epoch[jump_start])
        jump_thru_epoch = int(epoch[jump_stop - 1])
        # Duplicates are possible if multiple jumps match the same
        # storm

        #SA! if the rise event is increasing but at a lower rate than the threshold
        #SA! and the storm event is longer than the rise event, extend to the storm event while increasing

        #for size in range(int((jump_start_epoch-storm_start_epoch)/1800)):
        #    if zeta_mm[jump_start] > zeta_mm[jump_start-size]:
        #        jump_start_epoch = int(jump_start_epoch-size*1800)

        #for size in range(int((storm_thru_epoch-jump_thru_epoch)/1800)):
        #    if zeta_mm[jump_stop] > zeta_mm[jump_stop+size]:
        #        jump_thru_epoch = int(jump_thru_epoch+size*1800)

        # SA! Ensure minimum time window for P contribution to jump is the time window of that jump
        storm_start_epoch = int(np.min([storm_start_epoch,jump_start_epoch]))
        storm_thru_epoch = int(np.max([storm_thru_epoch-1800,jump_thru_epoch]))

        # SA! If the 2 hours before and/or after are still raining (P >= 0.5 mm/h) then extend the P time window
	    # SA! Still manual adjustment needed when used with in situ precipitation/throughfall.
        drizzling = 1
        if drizzling == 1:
            for i in range(4):
                if (rainfall_intensity_mm_h[jump_start-i-1]) >= (0.5):
                    i=i
                else:
                    break

            if ((i > 0) and (i<3)):
                storm_start_epoch = int(epoch[jump_start-(i)])
            elif (i == 3):
                storm_start_epoch = int(epoch[jump_start - (i+1)])
            else:
                storm_start_epoch = storm_start_epoch

            for j in range(4):
                if (rainfall_intensity_mm_h[jump_stop+(j)]) >= (0.5):
                    j=j
                else:
                    break

            if storm_thru_epoch < int(epoch[jump_stop+(j-1)]):
                storm_thru_epoch = int(epoch[jump_stop+(j)])
            else:
                storm_thru_epoch = storm_thru_epoch


        already_seen = cursor.execute(
            """
        SELECT EXISTS (
          SELECT 1 FROM storm
          WHERE start_epoch = :start_epoch
          AND thru_epoch = :thru_epoch
        )""",
            {"start_epoch": storm_start_epoch, "thru_epoch": storm_thru_epoch},
        ).fetchone()[0]
        assert not already_seen, (
            f"Storm at {convert_epoch_to_datetime_text(storm_start_epoch)}"
            f"--{convert_epoch_to_datetime_text(storm_thru_epoch)} UTC "
            "matched with more than one rise"
        )
        cursor.execute(
            """
        INSERT INTO storm (start_epoch, thru_epoch)
        SELECT :start_epoch, :thru_epoch""",
            {
                "start_epoch": storm_start_epoch,
                "thru_epoch": storm_thru_epoch,
            },
        )
        cursor.execute(
            """
        INSERT INTO zeta_interval
          (start_epoch, interval_type, thru_epoch)
        VALUES
          (:start_epoch, :interval_type, :thru_epoch)""",
            {
                "interval_type": "storm",
                "start_epoch": jump_start_epoch,
                "thru_epoch": jump_thru_epoch,
            },
        )
        cursor.execute(
            """
        INSERT INTO zeta_interval_storm
          (interval_start_epoch, interval_type, storm_start_epoch)
        VALUES
          (:interval_start_epoch, :interval_type, :storm_start_epoch)""",
            {
                "interval_type": "storm",
                "interval_start_epoch": jump_start_epoch,
                "storm_start_epoch": storm_start_epoch,
            },
        )


def match_storms(rain, head, rain_threshold, jump_threshold):
    """Match intervals of increasing head with rainstorms

    For contiguous intervals of head increasing faster than jump_threshold, look
    for overlapping contiguous periods of rain above rain_threshold.  Then, use
    a variant of the Gale-Shapley algorithm to find a set of matches between
    storms and rises that is stable in the sense that, by switching a pair of
    matches between storms and rises, one cannot improve the agreement in
    duration and start time for both matches.

    Returns (rain_intervals, head_intervals), two sequences of equal
    length containing the matching contiguous intervals.

    """
    is_raining = rain > rain_threshold
    rain_masks = list(get_true_interval_masks(is_raining))
    LOG.info(
        "%s intervals of rain above threshold %s mm / h",
        len(rain_masks),
        rain_threshold,
    )
    # initialize to -1; this will indicate interstorms
    storm_indices = np.zeros(len(is_raining), np.int64) - 1
    rain_mask = None
    for i, rain_mask in enumerate(rain_masks):
        storm_indices[rain_mask] = i
    del rain_mask
    assert ((storm_indices != -1) == is_raining).all()

    # Semantics: is_jump is True on an interval if the head was higher
    # at the end of the interval than it was at the beginning.
    # In other words, in the sequence [0, 1, 2, 2],
    # is_jump == [1, 1, 0] because the head values are recorded
    # at the boundaries between time intervals:
    #              *  *
    #           *
    # head = *
    #        |  |  |  |
    #       t0 t1 t2 t3
    # and head increased on the intervals (t0, t1] and
    # (t1, t2], but not on the interval (t2, t3].
    # In this case the (start, stop) slice corresponding to
    # the jump interval is (0, 2), suggesting head increase
    # from rain on the time interval (t0, t2]
    head_increments = np.diff(head)
    is_jump = head_increments > jump_threshold
    jump_masks = list(get_true_interval_masks(is_jump))
    LOG.info(
        "%s intervals of water level increment above threshold %s mm",
        len(jump_masks),
        jump_threshold,
    )
    rain_intervals = []
    head_intervals = []
    for i, jump_mask in enumerate(jump_masks):
        # No head was recorded at the end of the last interval, so
        # jump mask provides no information about rain in last
        # interval.
        intersection = is_raining[:-1] & jump_mask
        matching_storms = set(storm_indices[np.nonzero(intersection)[0]])
        assert -1 not in matching_storms, matching_storms
        for storm_index in matching_storms:
            rain_interval, head_interval = get_candidate_match_intervals(
                head, jump_threshold, is_raining, rain_masks, jump_mask, storm_index
            )
            rain_intervals.append(rain_interval)
            head_intervals.append(head_interval)
    rain_intervals, head_intervals = disambiguate_matching(
        rain_intervals, head_intervals
    )
    assert_equal(
        len(set(rain_intervals)),
        len(rain_intervals),
        "Rain intervals in matches not unique",
    )
    assert_equal(
        len(set(head_intervals)),
        len(head_intervals),
        "Head intervals in matches not unique",
    )
    return (rain_intervals, head_intervals)


def get_candidate_match_intervals(
    head, jump_threshold, is_raining, rain_masks, jump_mask, storm_index
):
    """Identify rain and head intervals corresponding to a match between jump and storms"""
    rain_indices = np.nonzero(rain_masks[storm_index])[0]
    assert is_raining[rain_indices].all(), "Storm includes only raining time steps"
    rain_start = rain_indices[0]
    rain_stop = rain_indices[-1] + 1
    assert (
        rain_start == 0 or not is_raining[rain_start - 1]
    ), "No heavy rain just before slice"
    assert not is_raining[rain_stop], "No heavy rain at end of slice"
    jump_indices = np.nonzero(jump_mask)[0]
    jump_start = jump_indices[0]
    jump_stop = jump_indices[-1] + 2
    assert (
        np.diff(head[jump_start:jump_stop]) > jump_threshold
    ).all(), "Head pair [{}, {}) meets jump threshold {} mm".format(
        jump_start, jump_stop, jump_threshold
    )
    assert (
        jump_start == 0 or head[jump_start] - head[jump_start - 1] <= jump_threshold
    ), "Jump <= {} starts at jump_start: {}, {}".format(
        jump_threshold, jump_start, head[jump_start] - head[jump_start - 1]
    )
    assert (
        jump_stop == len(head)
        or head[jump_stop] - head[jump_stop - 1] <= jump_threshold
    ), "Jump > {} ends at jump_stop: {}, {}".format(
        jump_threshold, jump_stop, head[jump_stop] - head[jump_stop - 1]
    )
    # At a minimum, all intervals include 1 rainfall intensity
    # value and 2 head values (note that these correspond to the
    # same duration: one time slice.
    assert rain_stop - rain_start >= 1
    assert jump_stop - jump_start >= 2
    return ((rain_start, rain_stop), (jump_start, jump_stop))


def disambiguate_matching(rain_intervals, jump_intervals):
    """Disambiguate matching between rain and jump intervals

    The input sequences of (start, stop) intervals for rain and head
    are of the same length, and in general constitute a many-to-many
    relation between rain intervals and jump intervals.  This function
    finds a stable one-to-one matching between rain intervals and jump
    intervals, returning them in the same data structure.

    Start and stop are integers.

    """
    assert len(rain_intervals) == len(jump_intervals)
    storm_stops = {rain_start: rain_stop for rain_start, rain_stop in rain_intervals}
    jump_stops = {jump_start: jump_stop for jump_start, jump_stop in jump_intervals}
    candidate_matches = [
        (rain_start, jump_start)
        for (rain_start, _), (jump_start, _) in zip(rain_intervals, jump_intervals)
    ]
    duration_differences = {
        (rain_start, jump_start): float(
            (rain_stop - rain_start) - (jump_stop - jump_start)
        )
        for (rain_start, rain_stop), (jump_start, jump_stop) in zip(
            rain_intervals, jump_intervals
        )
    }

    storms_dict = defaultdict(list)
    jumps_dict = defaultdict(list)
    for rain_start, jump_start in candidate_matches:
        storms_dict[rain_start].append(jump_start)
        jumps_dict[jump_start].append(rain_start)

    storm_candidates = {}
    for rain_start, jumps in storms_dict.items():

        def absolute_duration_difference(jump_start):
            """Get absolute difference in duration between storm and jump"""
            return -abs(duration_differences[(rain_start, jump_start)])

        candidates = sorted(jumps, key=absolute_duration_difference)
        storm_candidates[rain_start] = candidates

    jump_preferences = {}
    for jump_start, rains in jumps_dict.items():

        def time_offset(rain_start):
            """Compute the time offset between the start of storm and jump"""
            return jump_start - rain_start

        jump_preferences[jump_start] = {
            rain_start: -abs(time_offset(rain_start)) for rain_start in rains
        }
    # Mapping from jumps to storms
    jump_matches = find_stable_matching(storm_candidates, jump_preferences)
    unique_rain_intervals = []
    unique_jump_intervals = []
    for jump_start, rain_start in jump_matches.items():
        unique_rain_intervals.append((rain_start, storm_stops[rain_start]))
        unique_jump_intervals.append((jump_start, jump_stops[jump_start]))

    assert len(unique_rain_intervals) == len(unique_jump_intervals)
    return (unique_rain_intervals, unique_jump_intervals)


def find_stable_matching(storm_candidates, jump_preferences):
    """Find a stable matching between storms and jumps

    Uses a variant of the Gale-Shapley stable matching algorithm to
    find a matching between storms and jumps such that swapping two
    matches would not improve both matches.

    storms is a set of storms, each with the attribute "candidates"
    giving a sequence of candidate jumps, from worst to best.

    Each jump object has an attribute .preferences, which maps storms
    to a "match quality" for that storm (higher is better).

    Returns the stable matches as a dict mapping jumps to storms.

    """
    matchable_storms = {
        storm for storm, candidates in storm_candidates.items() if candidates
    }
    matches = dict()
    while matchable_storms:
        storm = matchable_storms.pop()
        # Inefficient but safe
        assert storm not in matches.items()
        storm_is_free = True
        assert storm_candidates[storm]
        jump = storm_candidates[storm].pop()
        if jump in matches:
            if jump_preferences[jump][storm] > jump_preferences[jump][matches[jump]]:
                assert matches[jump] not in matchable_storms
                #matchable_storms.add(matches[jump])
                #storm_candidates[matches[jump]]=[matches[jump]]
                matches[jump] = storm
                storm_is_free = False
        else:
            matches[jump] = storm
            storm_is_free = False
        if storm_is_free and storm_candidates[storm]:
            matchable_storms.add(storm)
    return matches


def check_for_uniform_time_steps(epoch):
    """Verify that time steps are uniform

    Checks that time steps in sorted array epoch are uniform,
    otherwise raising ValueError.

    Note that this may not work as expected for inexact datatypes
    (floats).

    """
    delta_t = np.diff(epoch)
    if delta_t.min() != delta_t.max():
        raise ValueError("Nonuniform time steps in {}".format(sorted(set(delta_t))))


def get_mystery_jump_mask(is_jump, is_raining, zeta_mm):
    """Flag everything from a mystery jump until the next rain

    A "mystery jump" is a jump in head with no rain.  The returned
    mask has the same length as is_jump and is_raining, and is marked
    True during intervals [jump with no rain, rain).  If it's raining,
    at the time of the jump, nothing happens; if not, all time points
    until (but not including) the next True value in is_raining are
    flagged.  All times with jumps and no rain are flagged; all times
    with rain and no jump are not flagged; and times with no rain and
    no jump will be flagged if and only if they occur between a
    "mystery jump" (jump with no rain) and a rain event.

    """

    """SA! added zeta_mm to stop mystery_jumps sooner. Instead of classifying as a mysteryjump untill is_raining=True
    now if it is not rising the mysterjump stops
    """

    assert_equal(len(is_jump), len(is_raining))
    mystery_jump_mask = np.zeros(len(is_jump), bool)
    in_mystery = True
    # pylint: disable=consider-using-enumerate
    for i in range(len(mystery_jump_mask)):
        if is_raining[i] or (zeta_mm[i] < zeta_mm[i-1]):
            in_mystery = False
        else:
            if is_jump[i]:
                in_mystery = True
        mystery_jump_mask[i] = in_mystery
    # Test assertions in docstring: mystery_jump_mask is true at least
    # wherever it's not raining and there's a jump, and false wherever
    # it's raining
    assert all(~mystery_jump_mask[np.nonzero(is_raining)[0]])
    assert all(mystery_jump_mask[np.nonzero(~is_raining & is_jump)[0]])
    return mystery_jump_mask


def get_true_interval_masks(boolean_vector):
    """Find a series of contigous intervals where boolean_vector is true

    Returns an iterator of masks for contiguous sets of indices where
    boolean_vector is true.  Each mask has the same length as boolean_vector.

    """
    if not boolean_vector.dtype == bool:
        raise ValueError("non-boolean input vector")
    int_vector = boolean_vector.astype(np.int64)
    is_start = np.concatenate(([0], int_vector[1:] - int_vector[:-1])) > 0
    # Indices increment for each block of True values in input vector
    indices = np.cumsum(is_start.astype(np.int64))
    indices[~boolean_vector] = 0
    assert (
        indices.astype(bool) == boolean_vector
    ).all(), "Non-zero indices correspond to True elements of input vector"
    unique_indices = sorted(set(indices))
    assert unique_indices[0] == 0
    del unique_indices[0]
    return (indices == index for index in unique_indices)


def assert_equal(a, b, message=None):  # pylint:disable=invalid-name
    """Convenience function for checking equality"""
    prefix = "{}: ".format(message) if message else ""
    assert a == b, "{}{} != {}".format(prefix, a, b)


def convert_epoch_to_datetime_text(epoch):
    """Convert UNIX epoch to ISO 8601 datetime text"""
    return datetime.datetime.fromtimestamp(epoch).strftime("%Y-%m-%d %H:%M:%S")
