import os

import numpy as np
import pandas as pd



def outcomes(statspath, tablespath):

    cols = ['Biomarker', 'Control', 'Treatment', 'Effect size', 'Effect size', 'p-value', 'range']
    units = ['', 'mL/min/100mL', 'mL/min/100mL', 'mL/min/100mL', '%', '', '%']

    file = os.path.join(statspath, f'ttest.csv')
    ttest = pd.read_csv(file)
    file = os.path.join(statspath, f'describe.csv')
    desc = pd.read_csv(file, index_col='parameter')

    rows = [
        units,
        [
            'khe',
            f"{np.round(6000 * desc.at['c_khe','mean'], 1)} +/- {np.round(6000 * desc.at['c_khe','95%CI'], 1)}", 
            f"{np.round(6000 * desc.at['d_khe','mean'], 1)} +/- {np.round(6000 * desc.at['d_khe','95%CI'], 1)}", 
            f"{np.round(6000 * desc.at['a_khe','mean'], 1)} +/- {np.round(6000 * desc.at['a_khe','95%CI'], 1)}", 
            f"{np.round(100 * desc.at['r_khe','mean'], 1)} +/- {np.round(100 * desc.at['r_khe','95%CI'], 1)}", 
            np.round(ttest[(ttest['A']=='c_khe') & (ttest['B']=='d_khe')]['p_unc'].values[0], 4),
            f"{np.round(100 * desc.at['r_khe','range'], 1)}",
        ],
        [   
            'kbh',
            f"{np.round(6000 * desc.at['c_kbh','mean'], 1)} +/- {np.round(6000 * desc.at['c_kbh','95%CI'], 1)}", 
            f"{np.round(6000 * desc.at['d_kbh','mean'], 1)} +/- {np.round(6000 * desc.at['d_kbh','95%CI'], 1)}", 
            f"{np.round(6000 * desc.at['a_kbh','mean'], 1)} +/- {np.round(6000 * desc.at['a_kbh','95%CI'], 1)}", 
            f"{np.round(100 * desc.at['r_kbh','mean'], 1)} +/- {np.round(100 * desc.at['r_kbh','95%CI'], 1)}", 
            np.round(ttest[(ttest['A']=='c_kbh') & (ttest['B']=='d_kbh')]['p_unc'].values[0], 4), 
            f"{np.round(100 * desc.at['r_kbh','range'], 1)}",
        ]
    ]

    os.makedirs(tablespath, exist_ok=True)
    file = os.path.join(tablespath, f'outcomes.csv')
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(file, index=False)



def outcomes_diurnal(statspath, tablespath, scans):

    if scans == 1:
        return
    
    cols = ['Biomarker', 'Visit', 'Start', 'End', 'Change', 'p-value', 'range']
    units = ['', '', 'mL/min/100mL', 'mL/min/100mL', '%', '', '%']

    file = os.path.join(statspath, f'ttest_diurnal.csv')
    ttest = pd.read_csv(file)
    file = os.path.join(statspath, f'describe.csv')
    desc = pd.read_csv(file, index_col='parameter')

    rows = [
        units,
        [
            'khe', 
            'Control',
            f"{np.round(6000 * desc.at['c_khe_i','mean'], 1)} +/- {np.round(6000 * desc.at['c_khe_i','95%CI'], 1)}", 
            f"{np.round(6000 * desc.at['c_khe_f','mean'], 1)} +/- {np.round(6000 * desc.at['c_khe_f','95%CI'], 1)}", 
            f"{np.round(100 * desc.at['c_dkhe','mean'], 1)} +/- {np.round(100 * desc.at['c_dkhe','95%CI'], 1)}",  
            np.round(ttest[ttest['parameter']=='c_dkhe']['p_val'].values[0], 4),
            f"{np.round(100 * desc.at['c_dkhe','range'], 1)}",
        ],
        [
            'khe', 
            'Treatment',
            f"{np.round(6000 * desc.at['d_khe_i','mean'], 1)} +/- {np.round(6000 * desc.at['d_khe_i','95%CI'], 1)}", 
            f"{np.round(6000 * desc.at['d_khe_f','mean'], 1)} +/- {np.round(6000 * desc.at['d_khe_f','95%CI'], 1)}", 
            f"{np.round(100 * desc.at['d_dkhe','mean'], 1)} +/- {np.round(100 * desc.at['d_dkhe','95%CI'], 1)}",  
            np.round(ttest[ttest['parameter']=='d_dkhe']['p_val'].values[0], 4),
            f"{np.round(100 * desc.at['d_dkhe','range'], 1)}",
        ],
        # [
        #     'kbh', 
        #     'Treatment',
        #     f"{np.round(6000 * desc.at['d_kbh_i','mean'], 1)} +/- {np.round(6000 * desc.at['d_kbh_i','95%CI'], 1)}", 
        #     f"{np.round(6000 * desc.at['d_kbh_f','mean'], 1)} +/- {np.round(6000 * desc.at['d_kbh_f','95%CI'], 1)}", 
        #     f"{np.round(100 * desc.at['d_dkbh','mean'], 1)} +/- {np.round(100 * desc.at['d_dkbh','95%CI'], 1)}",  
        #     np.round(ttest[ttest['parameter']=='d_dkbh']['p_val'].values[0], 4),
        #     f"{np.round(100 * desc.at['d_dkbh','range'], 1)}",
        # ],
    ]

    os.makedirs(tablespath, exist_ok=True)
    file = os.path.join(tablespath, f'outcomes_diurnal.csv')
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(file, index=False)


def constants(statspath, tablespath):

    file = os.path.join(statspath, f'describe.csv')
    desc = pd.read_csv(file, index_col='parameter')

    cols = ['biomarker', 'units', 'value (95%CI)']


    rows = [
        [
            'Cardiac output',
            'L/min',
            f"{np.round(desc.at['CO','mean'] * 60 / 1000, 1)} +/- {np.round(desc.at['CO','95%CI'] * 60 / 1000, 1)}", 
        ],
        [
            'Glomerular Filtration Rate',
            'mL/min',
            f"{np.round(desc.at['GFR','mean'] * 60, 1)} +/- {np.round(desc.at['GFR','95%CI'] * 60, 1)}", 
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

    os.makedirs(tablespath, exist_ok=True)
    file = os.path.join(tablespath, f'constants.csv')
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(file, index=False)