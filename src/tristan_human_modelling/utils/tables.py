import os

import numpy as np
import pandas as pd



def outcomes(statspath, tablespath):

    cols = ['Biomarker', 'Control', 'Treatment', 'Effect size', 'p-value', 'range']
    units = ['', 'mL/min/100mL', 'mL/min/100mL', '%', '', '%']

    file = os.path.join(statspath, f'ttest.csv')
    ttest = pd.read_csv(file)
    file = os.path.join(statspath, f'describe.csv')
    desc = pd.read_csv(file, index_col='parameter')

    rows = [
        units,
        [
            'k_e2h',
            f"{6000 * desc.at['k_1_e2h','mean']:.2g} +/- {6000 * desc.at['k_1_e2h','95%CI']:.2g}", 
            f"{6000 * desc.at['k_2_e2h','mean']:.2g} +/- {6000 * desc.at['k_2_e2h','95%CI']:.2g}", 
            f"{100 * desc.at['ke_e2h','mean']:.2g} +/- {100 * desc.at['ke_e2h','95%CI']:.2g}", 
            f"{ttest[(ttest['A']=='k_1_e2h') & (ttest['B']=='k_2_e2h')]['p_unc'].values[0]:.2g}",
            f"{100 * desc.at['ke_e2h','range']:.2g}",
        ],
        [   
            'k_h2b',
            f"{6000 * desc.at['k_1_h2b','mean']:.2g} +/- {6000 * desc.at['k_1_h2b','95%CI']:.2g}", 
            f"{6000 * desc.at['k_2_h2b','mean']:.2g} +/- {6000 * desc.at['k_2_h2b','95%CI']:.2g}", 
            f"{100 * desc.at['ke_h2b','mean']:.2g} +/- {100 * desc.at['ke_h2b','95%CI']:.2g}", 
            f"{ttest[(ttest['A']=='k_1_h2b') & (ttest['B']=='k_2_h2b')]['p_unc'].values[0]:.2g}", 
            f"{100 * desc.at['ke_h2b','range']:.2g}",
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
            'k_e2h', 
            'Control',
            f"{6000 * desc.at['ki_1_e2h','mean']:.2g} +/- {6000 * desc.at['ki_1_e2h','95%CI']:.2g}", 
            f"{6000 * desc.at['kf_1_e2h','mean']:.2g} +/- {6000 * desc.at['kf_1_e2h','95%CI']:.2g}", 
            f"{100 * desc.at['kc_1_e2h','mean']:.2g} +/- {100 * desc.at['kc_1_e2h','95%CI']:.2g}",  
            f"{ttest[ttest['parameter']=='kc_1_e2h']['p_val'].values[0]:.2g}",
            f"{100 * desc.at['kc_1_e2h','range']:.2g}",
        ],
        [
            'k_e2h', 
            'Treatment',
            f"{6000 * desc.at['ki_2_e2h','mean']:.2g} +/- {6000 * desc.at['ki_2_e2h','95%CI']:.2g}", 
            f"{6000 * desc.at['kf_2_e2h','mean']:.2g} +/- {6000 * desc.at['kf_2_e2h','95%CI']:.2g}", 
            f"{100 * desc.at['kc_2_e2h','mean']:.2g} +/- {100 * desc.at['kc_2_e2h','95%CI']:.2g}",  
            f"{ttest[ttest['parameter']=='kc_2_e2h']['p_val'].values[0]:.2g}",
            f"{100 * desc.at['kc_2_e2h','range']:.2g}",
        ],
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
            f"{desc.at['CO','mean'] * 60 / 1000:.2g} +/- {desc.at['CO','95%CI'] * 60 / 1000:.2g}", 
        ],
        [
            'Heart-lung mean transit time',
            'sec',
            f"{desc.at['T_hl','mean']:.2g} +/- {desc.at['T_hl','95%CI']:.2g}", 
        ],
        [
            'Heart-lung dispersion',
            '%',
            f"{100 * desc.at['D_hl','mean']:.2g} +/- {100 * desc.at['D_hl','95%CI']:.2g}", 
        ],
        [
            'Organs blood mean transit time',
            'sec',
            f"{desc.at['T_b_or','mean']:.2g} +/- {desc.at['T_b_or','95%CI']:.2g}", 
        ],
        [
            'Organs extraction fraction',
            '%',
            f"{100 * desc.at['E_or','mean']:.2g} +/- {100 * desc.at['E_or','95%CI']:.2g}", 
        ],
        [
            'Organs extravascular mean transit time',
            'min',
            f"{desc.at['T_e_or','mean'] / 60:.2g} +/- {desc.at['T_e_or','95%CI'] / 60:.2g}", 
        ],
        [
            'Organs venous return',
            '%',
            f"{100 * desc.at['vr_or','mean']:.2g} +/- {100 * desc.at['vr_or','95%CI']:.2g}", 
        ],
        [
            'Liver artery and gut mean transit time',
            'sec',
            f"{desc.at['T_lag','mean']:.2g} +/- {desc.at['T_lag','95%CI']:.2g}", #f"{x:.2g}"
        ],
        [
            'Liver fraction of cardiac output',
            '%',
            f"{100 * desc.at['fCO_li','mean']:.2g} +/- {100 * desc.at['fCO_li','95%CI']:.2g}", 
        ],
        # [
        #     'Gut dispersion',
        #     '%',
        #     f"{100 * desc.at['Dg','mean'], 0)} +/- {100 * desc.at['Dg','95%CI'], 0)}", 
        # ],
        [
            'Liver extracellular volume fraction',
            '%',
            f"{100 * desc.at['v_e_li','mean']:.2g} +/- {100 * desc.at['v_e_li','95%CI']:.2g}", 
        ],
    ]

    os.makedirs(tablespath, exist_ok=True)
    file = os.path.join(tablespath, f'constants.csv')
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(file, index=False)