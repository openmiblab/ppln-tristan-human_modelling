import os
import warnings

# Silence t-test warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

import pydmr
import numpy as np
import pingouin as pg
import pandas as pd


def describe(buildpath, drug):
    resultspath = os.path.join(buildpath, drug)
    tablespath = os.path.join(resultspath, 'Tables')
    os.makedirs(tablespath, exist_ok=True)
    file = os.path.join(resultspath, 'Results', 'all_results')
    dmr = pydmr.read(file, 'pandas') 

    summary = dmr['pars'].groupby("parameter")["value"].describe().reset_index()

    # Add other descriptors
    summary['IQR'] = summary['75%'].values - summary ['25%'].values
    summary['95%CI'] = 1.96 * summary['std'].values / np.sqrt(summary['count'].values)
    summary['95%CI/|mean|'] = summary['95%CI'].values / np.abs(summary ['mean'].values)
    summary['range'] = summary['max'].values - summary ['min'].values

    # Save results
    file = os.path.join(tablespath, f'describe.csv')
    summary.to_csv(file, index=False)


def ttest(buildpath, drug):
    resultspath = os.path.join(buildpath, drug)
    tablespath = os.path.join(resultspath, 'Tables')
    os.makedirs(tablespath, exist_ok=True)
    file = os.path.join(resultspath, 'Results', 'all_results')
    dmr = pydmr.read(file, 'pandas') 

    # Run pairwise (paired) tests across biomarker levels
    # - dv: dependent variable column name
    # - within: the factor whose levels you want to compare (biomarker)
    # - subject: column identifying repeated measures (subject_id)
    # - parametric=True -> t-tests (paired t)
    # - padjust -> multiple comparisons adjustment (e.g., 'fdr_bh', 'bonf')
    # - effsize -> compute effect size ('cohen' gives Cohen's d)
    res = pg.pairwise_tests(
        data=dmr['pars'],
        dv='value',
        within='parameter',
        subject='subject',
        parametric=True,
        padjust='none',   # or 'bonf', 'holm', 'none', etc.
        effsize='cohen'     # 'cohen', 'hedges', or None
    )

    # Save to CSV
    file = os.path.join(tablespath, f'ttest.csv')
    res.to_csv(file, index=False)


def ttest_translation(buildpath, drug):

    if drug == 'metformin':
        return
    if drug == 'rifampicin_clinical':
        return
    
    resultspath = os.path.join(buildpath, drug)
    tablespath = os.path.join(resultspath, 'Tables')
    os.makedirs(tablespath, exist_ok=True)
    file = os.path.join(resultspath, 'Results', 'all_results')
    dmr = pydmr.read(file, 'pandas') 

    # Rat reference values
    if drug == 'rifampicin':
        y_map = {
            'r_khe': -0.89,
            'r_kbh': -0.42,
        }
    elif drug == 'ciclosporin':
        y_map = {
            'r_khe': -0.96,
            'r_kbh': -0.58,
        }

    res = (
        dmr['pars']
        .loc[dmr['pars']['parameter'].isin(y_map)]
        .groupby('parameter')
        .apply(
            lambda x: pg.ttest(x['value'], y=y_map[x.name]), 
            include_groups=False
        )
    )

    # Save to CSV
    file = os.path.join(tablespath, f'ttest_translation.csv')
    res.to_csv(file, index=True)


def ttest_diurnal(buildpath, drug, scans):

    if scans == 1:
        return
    
    resultspath = os.path.join(buildpath, drug)
    tablespath = os.path.join(resultspath, 'Tables')
    os.makedirs(tablespath, exist_ok=True)
    file = os.path.join(resultspath, 'Results', 'all_results')
    dmr = pydmr.read(file, 'pandas') 

    # Test if mean parameter values are non-zero
    res = dmr['pars'].groupby('parameter').apply(
        lambda x: pg.ttest(x['value'], y=0),
        include_groups=False
    )

    # Save to CSV
    file = os.path.join(tablespath, f'ttest_diurnal.csv')
    res.to_csv(file, index=True)




def outcomes(buildpath, drug):

    resultspath = os.path.join(buildpath, drug)
    tablespath = os.path.join(resultspath, 'Tables')
    
    cols = ['Biomarker', 'Control', 'Treatment', 'Effect size', 'Effect size', 'p-value', 'range']
    units = ['', 'mL/min/100mL', 'mL/min/100mL', 'mL/min/100mL', '%', '', '%']

    file = os.path.join(tablespath, f'ttest.csv')
    ttest = pd.read_csv(file)
    file = os.path.join(tablespath, f'describe.csv')
    desc = pd.read_csv(file, index_col='parameter')

    rows = [
        units,
        [
            'khe',
            f"{np.round(6000 * desc.at['c_khe','mean'], 1)} +/- {np.round(6000 * desc.at['c_khe','95%CI'], 1)}", 
            f"{np.round(6000 * desc.at['d_khe','mean'], 1)} +/- {np.round(6000 * desc.at['d_khe','95%CI'], 1)}", 
            f"{np.round(6000 * desc.at['a_khe','mean'], 1)} +/- {np.round(6000 * desc.at['a_khe','95%CI'], 1)}", 
            f"{np.round(100 * desc.at['r_khe','mean'], 1)} +/- {np.round(100 * desc.at['r_khe','95%CI'], 1)}", 
            np.round(ttest[(ttest['A']=='c_khe') & (ttest['B']=='d_khe')]['p-unc'].values[0], 4),
            f"{np.round(100 * desc.at['r_khe','range'], 1)}",
        ],
        [   
            'kbh',
            f"{np.round(6000 * desc.at['c_kbh','mean'], 1)} +/- {np.round(6000 * desc.at['c_kbh','95%CI'], 1)}", 
            f"{np.round(6000 * desc.at['d_kbh','mean'], 1)} +/- {np.round(6000 * desc.at['d_kbh','95%CI'], 1)}", 
            f"{np.round(6000 * desc.at['a_kbh','mean'], 1)} +/- {np.round(6000 * desc.at['a_kbh','95%CI'], 1)}", 
            f"{np.round(100 * desc.at['r_kbh','mean'], 1)} +/- {np.round(100 * desc.at['r_kbh','95%CI'], 1)}", 
            np.round(ttest[(ttest['A']=='c_kbh') & (ttest['B']=='d_kbh')]['p-unc'].values[0], 4), 
            f"{np.round(100 * desc.at['r_kbh','range'], 1)}",
        ]
    ]

    file = os.path.join(tablespath, f'outcomes.csv')
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(file, index=False)



def outcomes_diurnal(buildpath, drug, scans):

    if scans == 1:
        return

    resultspath = os.path.join(buildpath, drug)
    tablespath = os.path.join(resultspath, 'Tables')
    
    cols = ['Biomarker', 'Visit', 'Start', 'End', 'Change', 'p-value', 'range']
    units = ['', '', 'mL/min/100mL', 'mL/min/100mL', '%', '', '%']

    file = os.path.join(tablespath, f'ttest_diurnal.csv')
    ttest = pd.read_csv(file)
    file = os.path.join(tablespath, f'describe.csv')
    desc = pd.read_csv(file, index_col='parameter')

    rows = [
        units,
        [
            'khe', 
            'Control',
            f"{np.round(6000 * desc.at['c_khe_i','mean'], 1)} +/- {np.round(6000 * desc.at['c_khe_i','95%CI'], 1)}", 
            f"{np.round(6000 * desc.at['c_khe_f','mean'], 1)} +/- {np.round(6000 * desc.at['c_khe_f','95%CI'], 1)}", 
            f"{np.round(100 * desc.at['c_dkhe','mean'], 1)} +/- {np.round(100 * desc.at['c_dkhe','95%CI'], 1)}",  
            np.round(ttest[ttest['parameter']=='c_dkhe']['p-val'].values[0], 4),
            f"{np.round(100 * desc.at['c_dkhe','range'], 1)}",
        ],
        [
            'khe', 
            'Treatment',
            f"{np.round(6000 * desc.at['d_khe_i','mean'], 1)} +/- {np.round(6000 * desc.at['d_khe_i','95%CI'], 1)}", 
            f"{np.round(6000 * desc.at['d_khe_f','mean'], 1)} +/- {np.round(6000 * desc.at['d_khe_f','95%CI'], 1)}", 
            f"{np.round(100 * desc.at['d_dkhe','mean'], 1)} +/- {np.round(100 * desc.at['d_dkhe','95%CI'], 1)}",  
            np.round(ttest[ttest['parameter']=='d_dkhe']['p-val'].values[0], 4),
            f"{np.round(100 * desc.at['d_dkhe','range'], 1)}",
        ],
        # [
        #     'kbh', 
        #     'Treatment',
        #     f"{np.round(6000 * desc.at['d_kbh_i','mean'], 1)} +/- {np.round(6000 * desc.at['d_kbh_i','95%CI'], 1)}", 
        #     f"{np.round(6000 * desc.at['d_kbh_f','mean'], 1)} +/- {np.round(6000 * desc.at['d_kbh_f','95%CI'], 1)}", 
        #     f"{np.round(100 * desc.at['d_dkbh','mean'], 1)} +/- {np.round(100 * desc.at['d_dkbh','95%CI'], 1)}",  
        #     np.round(ttest[ttest['parameter']=='d_dkbh']['p-val'].values[0], 4),
        #     f"{np.round(100 * desc.at['d_dkbh','range'], 1)}",
        # ],
    ]

    file = os.path.join(tablespath, f'outcomes_diurnal.csv')
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(file, index=False)


def constants(buildpath, drug):

    resultspath = os.path.join(buildpath, drug)
    tablespath = os.path.join(resultspath, 'Tables')

    file = os.path.join(tablespath, f'describe.csv')
    desc = pd.read_csv(file, index_col='parameter')

    cols = ['biomarker', 'units', 'value (95%CI)']


    rows = [
        [
            'Cardiac output',
            'L/min',
            f"{np.round(desc.at['CO','mean'] * 60 / 1000, 1)} +/- {np.round(desc.at['CO','95%CI'] * 60 / 1000, 1)}", 
        ],
        [
            'Heart-lung mean transit time',
            'sec',
            f"{np.round(desc.at['Thl','mean'], 0)} +/- {np.round(desc.at['Thl','95%CI'], 1)}", 
        ],
        [
            'Heart-lung dispersion',
            '%',
            f"{np.round(100 * desc.at['Dhl','mean'], 0)} +/- {np.round(100 * desc.at['Dhl','95%CI'], 0)}", 
        ],
        [
            'Organs blood mean transit time',
            'sec',
            f"{np.round(desc.at['To','mean'], 0)} +/- {np.round(desc.at['To','95%CI'], 0)}", 
        ],
        [
            'Organs extraction fraction',
            '%',
            f"{np.round(100 * desc.at['Eo','mean'], 0)} +/- {np.round(100 * desc.at['Eo','95%CI'], 0)}", 
        ],
        [
            'Organs extravascular mean transit time',
            'min',
            f"{np.round(desc.at['To_e','mean'] / 60, 1)} +/- {np.round(desc.at['To_e','95%CI'] / 60, 1)}", 
        ],
        [
            'Gut mean transit time',
            'sec',
            f"{np.round(desc.at['Tg','mean'], 0)} +/- {np.round(desc.at['Tg','95%CI'], 1)}", 
        ],
        [
            'Gut dispersion',
            '%',
            f"{np.round(100 * desc.at['Dg','mean'], 0)} +/- {np.round(100 * desc.at['Dg','95%CI'], 0)}", 
        ],
        [
            'Liver extracellular volume fraction',
            '%',
            f"{np.round(100 * desc.at['ve','mean'], 0)} +/- {np.round(100 * desc.at['ve','95%CI'], 0)}", 
        ],
    ]

    file = os.path.join(tablespath, f'constants.csv')
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(file, index=False)