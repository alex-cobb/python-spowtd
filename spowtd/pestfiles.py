"""Generate files for calibration with PEST

"""

import os

import yaml


def generate_rise_pestfiles(
    connection,
    parameter_file,
    outfile_type,
    configuration_file,
    outfile,
    precision=17,
):
    """Generate PEST files for calibration against rise curve"""
    assert outfile_type in ('tpl', 'ins', 'pst'), outfile_type
    parameters = yaml.safe_load(parameter_file)
    check_parameters(parameters)
    {
        'tpl': generate_rise_tpl_file,
        'ins': generate_rise_ins_file,
        'pst': generate_rise_pst_file,
    }[outfile_type](
        connection=connection,
        parameters=parameters,
        configuration=(
            {}
            if configuration_file is None
            else yaml.safe_load(configuration_file)
        ),
        outfile=outfile,
        precision=precision,
    )


def generate_curves_pestfiles(
    connection,
    parameter_file,
    outfile_type,
    configuration_file,
    outfile,
    precision=17,
):
    """Generate PEST files for calibration against master curves"""
    assert outfile_type in ('tpl', 'ins', 'pst'), outfile_type
    {
        'tpl': generate_curves_tpl_file,
        'ins': generate_curves_ins_file,
        'pst': generate_curves_pst_file,
    }[outfile_type](
        connection=connection,
        parameters=yaml.safe_load(parameter_file),
        configuration=(
            {}
            if configuration_file is None
            else yaml.safe_load(configuration_file)
        ),
        outfile=outfile,
        precision=precision,
    )


def generate_rise_tpl_file(
    connection, parameters, configuration, outfile, precision
):
    """Generate template file for calibration against rise curve"""
    lines = ['ptf @', 'specific_yield:']
    if parameters['specific_yield']['type'] == 'peatclsm':
        lines += [
            '  type: peatclsm',
            '  sd: @sd                      @',
            '  theta_s: @theta_s                 @',
            '  b: @b                       @',
            '  psi_s: @psi_s                   @',
        ]
    else:
        assert parameters['specific_yield']['type'] == 'spline'
        lines += ['  type: spline', '  zeta_knots_mm:']
        lines += [
            '    - {}'.format(value)
            for value in parameters['specific_yield']['zeta_knots_mm']
        ]
        lines += ['  sy_knots:  # Specific yield, dimensionless']
        lines += ['    - @sy_knot_{}@'.format(str(i).ljust(23))
                  for i in
                  range(1,
                        len(parameters['specific_yield']['sy_knots']) + 1)]
    lines += ['transmissivity:']
    if parameters['transmissivity']['type'] == 'peatclsm':
        lines += [
            '  type: peatclsm',
            '  Ksmacz0: {}  # m/s'.format(
                parameters['transmissivity']['Ksmacz0']
            ),
            '  alpha: {}  # dimensionless'.format(
                parameters['transmissivity']['alpha']
            ),
            '  zeta_max_cm: {}'.format(
                parameters['transmissivity']['zeta_max_cm']
            ),
        ]
    else:
        assert parameters['transmissivity']['type'] == 'spline'
        lines += ['  type: spline']
        lines += ['  zeta_knots_mm:']
        lines += [
            '    - {}'.format(value)
            for value in parameters['transmissivity']['zeta_knots_mm']
        ]
        lines += ['  K_knots_km_d:  # Conductivity, km /d']
        lines += [
            '    - {}'.format(value)
            for value in parameters['transmissivity']['K_knots_km_d']
        ]
        lines += [
            '  minimum_transmissivity_m2_d: {}  '
            '# Minimum transmissivity, m2 /d'.format(
                parameters['transmissivity']['minimum_transmissivity_m2_d']
            )
        ]
    outfile.write(os.linesep.join(lines))


def generate_rise_ins_file(
    connection, parameters, configuration, outfile, precision
):
    """Generate instruction file for calibration against rise curve"""
    cursor = connection.cursor()
    cursor.execute(
        """
    SELECT count(distinct zeta_number)
    FROM rising_interval_zeta"""
    )
    n_zeta = cursor.fetchone()[0]
    cursor.close()
    lines = ['pif @', '@# Rise curve simulation vector@']
    lines += ['l1 [e{}]3:24'.format(i + 1) for i in range(n_zeta)]
    outfile.write(os.linesep.join(lines))


def generate_rise_pst_file(
    connection, parameters, configuration, outfile, precision
):
    """Generate control file for calibration against rise curve"""
    # See Example 11.3 in pestman and Preface of addendum
    parameterization = parameters['specific_yield']['type']
    if parameterization not in ('peatclsm', 'spline'):
        raise ValueError(
            'Unrecognized parameterization "{}"'.format(parameterization)
        )
    if parameterization == 'spline':
        npar = len(parameters['specific_yield']['sy_knots'])
        npargp = 1  # One parameter group
    else:
        npar = 4
        npargp = 4  # A parameter group for each parameter
    cursor = connection.cursor()
    cursor.execute(
        """
    SELECT count(distinct zeta_number)
    FROM rising_interval_zeta"""
    )
    n_zeta = cursor.fetchone()[0]
    cursor.execute(
        """
    SELECT mean_crossing_depth_mm AS dynamic_storage_mm
    FROM average_rising_depth
    ORDER BY zeta_mm"""
    )
    avg_storage_mm = [row[0] for row in cursor.fetchall()]
    cursor.close()
    nobsgp = 1  # 1 observation group
    lines = [
        'pcf',
        '* control data',
        'restart  estimation',
        (
            '{}{}{}     0{}'.format(
                str(npar).rjust(5),
                str(n_zeta).rjust(6),
                str(npargp).rjust(6),
                str(nobsgp).rjust(6),
            )
        ),
        '    1     1 double point   1   0   0',
        '   5.0  2.0   0.3  0.03    10',
        '  3.0   3.0 0.001  0',
        '  0.1',
        '   30  0.01     4     3  0.01     3',
        '    1     1     1',
    ]
    if parameterization == 'spline':
        lines += [
            '* parameter groups',
            'sy_knot      relative 0.01  0.0  switch  2.0 parabolic',
            '* parameter data',
        ]
        lines += [
            'sy_knot_{}   none relative   NaN  0.01  1    sy_knot    1.0  0.0 1'.format(
                i + 1
            )
            for i in range(npar)
        ]
    else:
        lines += [
            '* parameter groups',
            'sd           relative 0.01  0.0  switch  2.0 parabolic',
            'theta_s      relative 0.01  0.0  switch  2.0 parabolic',
            'b            relative 0.01  0.0  switch  2.0 parabolic',
            'psi_s        relative 0.01  0.0  switch  2.0 parabolic',
        ]
        lines += [
            '* parameter data',
            'sd          none relative   NaN  0.0   2.0  sd         1.0  0.0 1',
            'theta_s     none relative   NaN  0.01  1    theta_s    1.0  0.0 1',
            'b           none relative   NaN  0.01  20.0 b          1.0  0.0 1',
            'psi_s       none relative   NaN  -1.0  -0.01  psi_s      1.0  0.0 1',
        ]
    lines += ['* observation groups', 'storageobs']
    lines += ['* observation data']
    lines += [
        (
            'e{{}}    {{:0.{}g}}    1.0   storageobs'.format(precision).format(
                i + 1, W
            )
        )
        for i, W in enumerate(avg_storage_mm)
    ]
    lines += [
        '* model command line',
        # XXX
        'bash simulate-rise.sh',
    ]
    lines += [
        '* model input/output',
        # XXX
        'rise_pars.yml.tpl  rise_pars.yml',
        'rise_observations.ins  rise_observations.yml',
    ]
    lines += ['* prior information']
    outfile.write(os.linesep.join(lines))


def generate_curves_tpl_file(
    connection, parameters, configuration, outfile, precision
):
    """Generate template file for calibration against master curves"""
    lines = ['ptf @', 'specific_yield:']
    if parameters['specific_yield']['type'] == 'peatclsm':
        lines += [
            '  type: peatclsm',
            '  sd: @sd                      @',
            '  theta_s: @theta_s                 @',
            '  b: @b                       @',
            '  psi_s: @psi_s                   @',
        ]
    else:
        assert parameters['specific_yield']['type'] == 'spline'
        lines += ['  type: spline', '  zeta_knots_mm:']
        lines += [
            '    - {}'.format(value)
            for value in parameters['specific_yield']['zeta_knots_mm']
        ]
        lines += ['  sy_knots:  # Specific yield, dimensionless']
        lines += ['    - @sy_knot_{}@'.format(str(i).ljust(23))
                  for i in
                  range(1,
                        len(parameters['specific_yield']['sy_knots']) + 1)]
    lines += ['transmissivity:']
    if parameters['transmissivity']['type'] == 'peatclsm':
        lines += [
            '  type: peatclsm',
            '  Ksmacz0: @Ksmacz0                 @  # m/s',
            '  alpha: @alpha                   @  # dimensionless',
            '  zeta_max_cm: {}'.format(
                parameters['transmissivity']['zeta_max_cm']
            ),
        ]
    else:
        assert parameters['transmissivity']['type'] == 'spline'
        lines += ['  type: spline']
        lines += ['  zeta_knots_mm:']
        lines += [
            '    - {}'.format(value)
            for value in parameters['transmissivity']['zeta_knots_mm']
        ]
        lines += ['  K_knots_km_d:  # Conductivity, km /d']
        lines += [
            '    - @K_knot_{}@'.format(str(i).ljust(17))
            for i in range(
                1, len(parameters['transmissivity']['K_knots_km_d']) + 1
            )
        ]
        lines += [
            '  minimum_transmissivity_m2_d: {}  '
            '# Minimum transmissivity, m2 /d'.format(
                '@T_min                   @'
            )
        ]
    outfile.write(os.linesep.join(lines))


def generate_curves_ins_file(
    connection, parameters, configuration, outfile, precision
):
    """Generate instruction file for calibration against master curves"""
    cursor = connection.cursor()
    cursor.execute(
        """
    SELECT count(distinct zeta_number)
    FROM rising_interval_zeta"""
    )
    n_rise_zeta = cursor.fetchone()[0]
    cursor.execute(
        """
    SELECT count(distinct zeta_number)
    FROM recession_interval_zeta"""
    )
    n_recession_zeta = cursor.fetchone()[0]
    cursor.close()
    lines = ['pif @', '@# Rise curve simulation vector@']
    lines += ['l1 [e{}]3:24'.format(i + 1) for i in range(n_rise_zeta)]
    lines += ['@# Recession curve simulation vector@']
    lines += [
        'l1 [e{}]3:24'.format(i + 1)
        for i in range(n_rise_zeta, n_rise_zeta + n_recession_zeta)
    ]
    outfile.write(os.linesep.join(lines))


def generate_curves_pst_file(
    connection, parameters, configuration, outfile, precision
):
    """Generate control file for calibration against master curves"""
    # See Example 11.3 in pestman and Preface of addendum
    parameterization = parameters['specific_yield']['type']
    if parameterization not in ('peatclsm', 'spline'):
        raise ValueError(
            'Unrecognized parameterization "{}"'.format(parameterization)
        )
    if parameterization == 'spline':
        n_Sy = len(parameters['specific_yield']['sy_knots'])
        n_T = len(parameters['transmissivity']['K_knots_km_d'])
        npar = n_Sy + n_T + 1  # For minimum transmissivity
        npargp = 3  # Three parameter groups
    else:
        assert parameterization == 'peatclsm'
        npar = 6
        npargp = 6  # A parameter group for each parameter
    cursor = connection.cursor()
    cursor.execute(
        """
    SELECT mean_crossing_depth_mm AS dynamic_storage_mm
    FROM average_rising_depth
    ORDER BY zeta_mm"""
    )
    avg_storage_mm = [row[0] for row in cursor.fetchall()]
    n_rise_zeta = len(avg_storage_mm)
    # Sort from highest to lowest water level, to match the way
    # recession curves are dumped and plotted.
    cursor.execute(
        """
    SELECT CAST(elapsed_time_s AS double precision)
             / (3600 * 24) AS elapsed_time_d
    FROM average_recession_time
    ORDER BY zeta_mm DESC"""
    )
    elapsed_time_d = [row[0] for row in cursor.fetchall()]
    n_recession_zeta = len(elapsed_time_d)
    cursor.close()
    nobsgp = 2  # 2 observation groups
    lines = [
        'pcf',
        '* control data',
        'restart  estimation',
        (
            '{}{}{}     0{}'.format(
                str(npar).rjust(5),
                str(n_rise_zeta + n_recession_zeta).rjust(6),
                str(npargp).rjust(6),
                str(nobsgp).rjust(6),
            )
        ),
        '    1     1 double point   1   0   0',
        '   5.0  2.0   0.3  0.03    10',
        '  3.0   3.0 0.001  0',
        '  0.1',
        '   30  0.01     4     3  0.01     3',
        '    1     1     1',
    ]
    if parameterization == 'spline':
        lines += [
            '* parameter groups',
            'sy_knot      relative 0.01  0.0  switch  2.0 parabolic',
            'k_knot       relative 0.01  0.0  switch  2.0 parabolic',
            'T_min        relative 0.01  0.0  switch  2.0 parabolic',
            '* parameter data',
        ]
        lines += [
            'sy_knot_{}  none relative  NaN  0.01     1       sy_knot  1.0  0.0  1'.format(
                i + 1
            )
            for i in range(n_Sy)
        ]
        lines += [
            'k_knot_{}   log  factor    NaN  1.0e-04  1.0e+5  k_knot   1.0  0.0  1'.format(
                i + 1
            )
            for i in range(n_T)
        ]
        lines += [
            'T_min      log  factor    NaN  1.0e-04  1.0e+5  T_min    1.0  0.0  1'
        ]
    else:
        lines += [
            '* parameter groups',
            'sd           relative 0.01  0.0  switch  2.0 parabolic',
            'theta_s      relative 0.01  0.0  switch  2.0 parabolic',
            'b            relative 0.01  0.0  switch  2.0 parabolic',
            'psi_s        relative 0.01  0.0  switch  2.0 parabolic',
            'Ksmacz0      relative 0.01  0.0  switch  2.0 parabolic',
            'alpha        relative 0.01  0.0  switch  2.0 parabolic',
        ]
        lines += [
            '* parameter data',
            'sd          none relative   NaN  0.0      2.0        sd         1.0  0.0  1',
            'theta_s     none relative   NaN  0.01     1          theta_s    1.0  0.0  1',
            'b           none relative   NaN  0.01     20.0       b          1.0  0.0  1',
            'psi_s       none relative   NaN  -1.0     -0.01      psi_s      1.0  0.0  1',
            'Ksmacz0     log  factor     NaN  1.0e-04  1.0e+5  Ksmacz0    1.0  0.0  1',
            'alpha       none relative   NaN  1        20.0       alpha      1.0  0.0  1',
        ]
    lines += ['* observation groups', 'storageobs', 'timeobs']
    lines += ['* observation data']
    lines += [
        (
            'e{{}}    {{:0.{}g}}    1.0   storageobs'.format(precision).format(
                i + 1, W
            )
        )
        for i, W in enumerate(avg_storage_mm)
    ]
    lines += [
        (
            'e{{}}    {{:0.{}g}}    1.0   timeobs'.format(precision).format(
                n_rise_zeta + i + 1, t
            )
        )
        for i, t in enumerate(elapsed_time_d)
    ]
    lines += [
        '* model command line',
        # XXX
        'bash simulate-curves.sh',
    ]
    lines += [
        '* model input/output',
        # XXX
        'curves_pars.yml.tpl  curves_pars.yml',
        'curves_observations.ins  curves_observations.yml',
    ]
    lines += ['* prior information']
    outfile.write(os.linesep.join(lines))


def check_parameters(parameters):
    """Check parameters for correctness"""
    if parameters['specific_yield']['type'] not in ('peatclsm', 'spline'):
        raise ValueError(
            'Unexpected specific yield type: {}'.format(
                parameters['specific_yield']['type']
            )
        )
    if parameters['transmissivity']['type'] not in ('peatclsm', 'spline'):
        raise ValueError(
            'Unexpected specific yield type: {}'.format(
                parameters['transmissivity']['type']
            )
        )
