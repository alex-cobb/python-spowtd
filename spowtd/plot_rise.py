#!/usr/bin/python3

"""Plot master rise curve

"""

import itertools
import matplotlib
import matplotlib.pyplot as plt

import numpy as np

import yaml

import spowtd.specific_yield as specific_yield_mod
import spowtd.simulate_rise as simulate_rise_mod

DEFAULT_COLORS = {
    'interstorm': '#CC79A760',
    'zeta_storm': '#00929260',
    'rain_storm': '#0000ff22',
    'rain': '#ff0000ff',
    'jump': '#0000ffff',
    'mystery_jump': '#ff00ffff',
}


def plot_rise(connection, parameters):
    """Plot master rise curve"""
    cursor = connection.cursor()

    fig = plt.figure(figsize=(8.96,6.72))
    axes = fig.add_subplot(1, 1, 1)
    axes.set_xlabel('Dynamic storage (cm)', fontsize=18)
    axes.set_ylabel('$\overline{z}_{WL}$ (cm)', fontsize=18)
    legendLines = [matplotlib.lines.Line2D([0], [0], color='#00929290', linewidth=2), matplotlib.lines.Line2D([0], [0], color='black', linewidth=2), matplotlib.lines.Line2D([0], [0], color='#490092', linewidth=2, linestyle='--')]
    axes.legend(legendLines, ['Individual rise events',  'Master rise curve', 'Fitted rise curve'],loc="lower right", fontsize=16)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.axhline(y=0, xmin=-100, xmax=100, color='grey', linestyle='--')

    #SA! 1 is for stretched plotting AND 0 is for non strecthed plottin
    IMERG = 0
    if IMERG == 1:
        cursor.execute(
            """
        SELECT interval_start_epoch AS rise_interval,
               rain_depth_offset_mm / 10 AS storage_cm,
               stretched_initial_zeta_mm / 10 AS zeta_cm
        FROM rising_curve_line_segment
        UNION ALL
        SELECT interval_start_epoch AS rise_interval,
               (rain_depth_offset_mm + stretched_rain_total_depth_mm) / 10 AS storage_cm,
               stretched_final_zeta_mm / 10 AS zeta_cm
        FROM rising_curve_line_segment
        ORDER BY rise_interval, storage_cm"""
        )
        for _, group in itertools.groupby(
            cursor.fetchall(), key=lambda row: row[0]
        ):
            (storage_cm, zeta_cm) = zip(*((t, z) for _, t, z in group))
            axes.plot(storage_cm, zeta_cm, '-', color='#00929260')

        cursor.execute(
            """
        SELECT mean_crossing_depth_mm / 10 AS storm_depth_cm,
               zeta_mm / 10 AS zeta_cm
        FROM average_rising_depth
        ORDER BY zeta_mm"""
        )
    else:
        cursor.execute(
            """
        SELECT interval_start_epoch AS rise_interval,
               rain_depth_offset_mm / 10 AS storage_cm,
               initial_zeta_mm / 10 AS zeta_cm
        FROM rising_curve_line_segment
        UNION ALL
        SELECT interval_start_epoch AS rise_interval,
               (rain_depth_offset_mm + rain_total_depth_mm) / 10 AS storage_cm,
               final_zeta_mm / 10 AS zeta_cm
        FROM rising_curve_line_segment
        ORDER BY rise_interval, storage_cm"""
        )
        for _, group in itertools.groupby(
                cursor.fetchall(), key=lambda row: row[0]
        ):
            (storage_cm, zeta_cm) = zip(*((t, z) for _, t, z in group))
            axes.plot(storage_cm, zeta_cm, '-', color='#00929260')

        cursor.execute(
            """
        SELECT mean_crossing_depth_mm / 10 AS storm_depth_cm,
               zeta_mm / 10 AS zeta_cm
        FROM average_rising_depth
        ORDER BY zeta_mm"""
        )
    (avg_storm_depth_cm, avg_zeta_cm) = zip(*cursor)

    axes.plot(avg_storm_depth_cm, avg_zeta_cm, '-k')

    if parameters is not None:
        plot_simulated_rise(
            axes, parameters, avg_zeta_cm, np.mean(avg_storm_depth_cm)
        )

    cursor.close()
    fname_long =('/data/leuven/324/vsc32460/AC/spowtd/FIG/' + 'rise_Itanga.png') #vsc-hard-mounts/leuven-data/
    plt.savefig(fname_long, dpi=300)
    plt.show()
    return 0


def plot_simulated_rise(axes, parameters, zeta_grid_cm, mean_storage_cm):
    """Plot simulated rise curve"""
    sy_parameters = yaml.safe_load(parameters)['specific_yield']
    W_mm = simulate_rise_mod.compute_rise_curve(
        specific_yield=specific_yield_mod.create_specific_yield_function(
            sy_parameters
        ),
        zeta_grid_mm=np.array(zeta_grid_cm, dtype=float) * 10,
        mean_storage_mm=mean_storage_cm * 10,
    )
    axes.plot(W_mm / 10, zeta_grid_cm, '--', color='#490092', linewidth=2.5)
