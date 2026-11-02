import os
import logging

import numpy as np
import dcmri as dc

CONFIG = {
    'heartlung': 'chain',          # Smoother first pass data
    'organs': '2cxm',              # Long acquisition times reveal extravasation
    'lagut': 'comp',
    'liver': '1I-IC',              # Intracellular agent
    'baseline': 'measured',        # Data are available on precontrast T1
    'sequence': 'ZTE-3D-SPGR-SS',  # Echo times are assumed to be zero
    't2s_relaxation_ao': None,     # No T2* relaxation with TE=0
    't2s_relaxation_li': None,     # No T2* relaxation with TE=0
    'calibrate': True,             # Derive signal scaling factor from data
    'trigger': True,               # Signals are retrospectively triggered
    'estimate_bat': 'ttp',         # Use the dedicated BAT estimator
}

INIT = { 
    'dt': 0.1, 
    'dose_tolerance': 0.01,
    'H': 0.45,

    'fCO_li': 0.10,
    'CO': 100,
    'T_hl': 10,
    'D_hl': 0.5,
    'T_b_or': 20,
    'E_or': 0.15,
    'vr_or': 0.9,
    'T_e_or': 120,
    'T_lag': 30,
    #'D_lag': 0.5,
    'v_e_li': 0.25,
    'T_1_h': 4 * 60 * 60,
    'T_2_h': 4 * 60 * 60,
}

FREE = {
    'BAT_1': (-30, 30),
    'BAT_2': (-30, 30),

    'fCO_li': (0.01, 0.99),
    'CO': (30, 300),
    'T_hl': (0, 30),
    'D_hl': (0.2, 0.8),
    'T_b_or': (0, 60),
    'E_or': (0, 0.5),
    'vr_or': (0.5, 1.0),
    'T_e_or': (0, 800), 
    'T_lag': (0, 120),
    #'D_lag': (0, 1),
    'v_e_li': (0.01, 0.5),
    'T_1_h': (600, 6 * 60 * 60),
    'T_2_h': (600, 6 * 60 * 60),
}


def one_scan(subj, data, dir_results):
    config = CONFIG
    data |= INIT | { 
        'E_1_li': 0.1,
        'E_2_li': 0.1,
    }
    data['pfree'] = FREE | {
        'E_1_li': (0, 1),
        'E_2_li': (0, 1),
    }
    Model = dc.InverseAortaLiverDrug
    study='one scan - two visits'

    _fit_model(subj, data, dir_results, config, Model, study)



def two_scan(subj, data, dir_results):
    config = CONFIG | {
        'bolus': 'dual',               # Two independent injections performed
        'non_stationary': 'U',         # Uptake rate changes during the course of the experiment
    }
    data |= INIT | { 
        'Ei_1_li': 0.1,
        'Ei_2_li': 0.1,
        'Ef_1_li': 0.1,
        'Ef_2_li': 0.1,
    }
    data['pfree'] = FREE | {
        'BAT_3': (-30, 30),
        'BAT_4': (-30, 30),
        'Ei_1_li': (0, 1),
        'Ei_2_li': (0, 1),
        'Ef_1_li': (0, 1),
        'Ef_2_li': (0, 1),
    }
    study = 'two scans - two visits'
    Model = dc.InverseAortaLiverSplitDrug

    _fit_model(subj, data, dir_results, config, Model, study)



def _fit_model(subj, data, dir_results, config, Model, study):

    # Prepare output files
    pars_path = os.path.join(dir_results, 'Pars')
    plot_path = os.path.join(dir_results, 'Plots')
    pcov_path = os.path.join(dir_results, 'Cov')
    pvar_path = os.path.join(dir_results, 'Sdev')
    os.makedirs(pars_path, exist_ok=True)
    os.makedirs(plot_path, exist_ok=True)
    os.makedirs(pcov_path, exist_ok=True)
    os.makedirs(pvar_path, exist_ok=True)

    plot_file_all = os.path.join(plot_path, subj + '.png')
    plot_file_init = os.path.join(plot_path, subj + '_first_pass.png')
    plot_file_pcov = os.path.join(pcov_path, subj + '.png')
    plot_file_pvar = os.path.join(pvar_path, subj + '.png')
    pars_file = os.path.join(pars_path, subj + '.dmr')

    if os.path.exists(pars_file + '.zip'):
        return

    logging.info(f"Training subject {subj}")

    try:
        model = Model(**config)
        result = model(data, verbose=2, xtol=1e-3)
        logging.info(f"Success training subject {subj}")
    except:
        logging.exception(f"Failed training subject {subj}")
        return False

    model.save(result, pars_file, subj=subj, study=study)
    model.plot(data, result, fname=plot_file_all, show=False)
    model.plot(data, result, fname=plot_file_init, show=False, xlim=[0, 5 * 60])
    model.plot_pcov(result, plot_file_pcov, show=False)
    model.plot_pvar(result, plot_file_pvar, show=False, xlim=[-3, 3])
    
    return True