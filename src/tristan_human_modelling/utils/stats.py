import os


import pydmr
import numpy as np
import pingouin as pg


def describe(file, tablespath):
    dmr = pydmr.read(file, 'pandas') 

    summary = dmr['pars'].groupby("parameter")["value"].describe().reset_index()

    # Add other descriptors
    summary['IQR'] = summary['75%'].values - summary ['25%'].values
    summary['95%CI'] = 1.96 * summary['std'].values / np.sqrt(summary['count'].values)
    summary['95%CI/|mean|'] = summary['95%CI'].values / np.abs(summary ['mean'].values)
    summary['range'] = summary['max'].values - summary ['min'].values

    # Save results
    file = os.path.join(tablespath, f'describe.csv')
    os.makedirs(tablespath, exist_ok=True)
    summary.to_csv(file, index=False)


def ttest(file, tablespath):
    # TODO: handle warnings - remove constants from list of variables
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
    os.makedirs(tablespath, exist_ok=True)
    res.to_csv(file, index=False)


def ttest_translation(file, tablespath, drug):
    
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
    os.makedirs(tablespath, exist_ok=True)
    file = os.path.join(tablespath, f'ttest_translation.csv')
    res.to_csv(file, index=True)


def ttest_diurnal(file, tablespath):
    # TODO: handle warnings - remove constants from list of variables
    dmr = pydmr.read(file, 'pandas') 

    # Test if mean parameter values are non-zero
    res = dmr['pars'].groupby('parameter').apply(
        lambda x: pg.ttest(x['value'], y=0),
        include_groups=False
    )
    # Save to CSV
    os.makedirs(tablespath, exist_ok=True)
    file = os.path.join(tablespath, f'ttest_diurnal.csv')
    res.to_csv(file, index=True)