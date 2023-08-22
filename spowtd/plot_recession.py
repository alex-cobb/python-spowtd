#!/usr/bin/python3

"""Plot master recession curve

"""

import itertools

import matplotlib
import matplotlib.pyplot as plt

import numpy as np

import yaml

import spowtd.simulate_recession as simulate_recession_mod


def plot_recession(connection, parameters):
    """Plot master recession curve"""
    cursor = connection.cursor()

    fig = plt.figure(figsize=(8.96,6.72))
    axes = fig.add_subplot(1, 1, 1)
    axes.set_xlabel('Elapsed time (day)', fontsize=18)
    axes.set_ylabel('$\overline{z}_{WL}$ (cm)', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    legendLines = [matplotlib.lines.Line2D([0], [0], color='#CC79A7', linewidth=2), matplotlib.lines.Line2D([0], [0], color='black', linewidth=2), matplotlib.lines.Line2D([0], [0], color='#490092', linewidth=2, linestyle='--')]
    axes.legend(legendLines, [ 'Individual recession events', 'Master recession curve', 'Fitted recession curve'],loc="lower left", fontsize=16)
    plt.axhline(y=0, xmin=-100, xmax=100, color='grey', linestyle='--') 

    cursor.execute(
        """
    SELECT start_epoch AS zeta_interval,
           (epoch - start_epoch + time_offset_s) / (3600. * 24)
             AS elapsed_time_d,
           zeta_mm / 10 AS zeta_cm
    FROM zeta_interval AS zi
    JOIN recession_interval AS ri
      USING (start_epoch, interval_type)
    JOIN water_level AS wl
      ON wl.epoch >= zi.start_epoch
      AND wl.epoch <= zi.thru_epoch
    WHERE interval_type = 'interstorm'
    ORDER BY start_epoch, elapsed_time_d"""
    )
    for _, group in itertools.groupby(
        cursor.fetchall(), key=lambda row: row[0]
    ):
        (elapsed_time_d, zeta_cm) = zip(*((t, z) for _, t, z in group))
        axes.plot(elapsed_time_d, zeta_cm, '-', color='#CC79A7')

    cursor.execute(
        """
    SELECT CAST(elapsed_time_s AS double precision)
             / (3600 * 24) AS elapsed_time_d,
           zeta_mm / 10 AS zeta_cm
    FROM average_recession_time"""
    )
    (avg_elapsed_time_d, avg_zeta_cm) = zip(*cursor)
    cursor.close()

    axes.plot(avg_elapsed_time_d, avg_zeta_cm,'-', color='black')

    if parameters is not None:
        (_, _, elapsed_time_d) = simulate_recession_mod.simulate_recession(
            connection, parameters
        )
        axes.plot(elapsed_time_d, avg_zeta_cm, '--', color='#490092', linewidth=2.5)
    fname_long =('/data/leuven/324/vsc32460/AC/spowtd/FIG/' + 'recession_Itanga.png')
    plt.savefig(fname_long, dpi=300)
    plt.show()
    return 0
