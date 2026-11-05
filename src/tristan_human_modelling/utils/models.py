import os
import logging

import numpy as np
import dcmri as dc

CONFIG = {
    'heartlung': 'chain',          # Smoother first pass data
    'organs': '2cxm',              # Long acquisition times reveal extravasation
    'lagut': 'comp',
    #'lagut': 'plucom',
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

    'fCO_li': 0.20,
    'CO': 100,
    'T_hl': 10,
    'D_hl': 0.05,
    'T_b_or': 5,
    'E_or': 0.25,
    'vr_or': 0.9,
    'T_e_or': 120,

    'T_lag': 30,
    #'D_lag': 0.5,
    # 'T_la': 1,
    # 'T_gu': 20,
    # 'ffa': 0.2,

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
    'D_hl': (0.01, 0.99),
    'T_b_or': (2, 60),
    'E_or': (0, 0.5),
    'vr_or': (0.5, 1.0),
    'T_e_or': (0, 800), 

    'T_lag': (0, 120),
    #'D_lag': (0, 1),
    # 'T_la': (0, 5),
    # 'T_gu': (0, 120),
    # 'ffa': (0, 1),

    'v_e_li': (0.01, 0.5),
    'T_1_h': (600, 6 * 60 * 60),
    'T_2_h': (600, 6 * 60 * 60),
}


def one_scan(subj, data, subjects, dir_results):
    config = CONFIG
    data |= INIT | { 
        'E_1_li': 0.1,
        'E_2_li': 0.1,
    }
    data['pfree'] = FREE | {
        'E_1_li': (0, 1),
        'E_2_li': (0, 1),
    }
    model = dc.InverseAortaLiverDrug(**config)
    study = subjects[subj]

    return _fit_model(subj, data, dir_results, model, study, scans=1)



def two_scan(subj, data, subjects, dir_results, dir_results_one_scan=None):
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
    study = subjects[subj]
    model = dc.InverseAortaLiverSplitDrug(**config)


    # Optional: Initialize with one-scan model results
    model_one_scan = dc.InverseAortaLiverDrug(**CONFIG)
    pars_one_scan = os.path.join(dir_results_one_scan, 'Pars', subj + '.dmr')
    results_one_scan, _, _ = model_one_scan.load(pars_one_scan)
    data |= results_one_scan['popt'] # BAT is incorrect but can be ignored
    data |= { 
        'Ei_1_li': data['E_1_li'],
        'Ei_2_li': data['E_2_li'],
        'Ef_1_li': data['E_1_li'],
        'Ef_2_li': data['E_2_li'],
    }

    # Fit
    return _fit_model(subj, data, dir_results, model, study, scans=2)



def _fit_model(subj, data, dir_results, model, study, scans=None):

    # Prepare output files
    pars_path = os.path.join(dir_results, 'Pars')
    pars_file = os.path.join(pars_path, subj + '.dmr')
    os.makedirs(pars_path, exist_ok=True)

    if os.path.exists(pars_file + '.zip'):
        return

    logging.info(f"Training subject {subj}")

    try:
        result = model(data, verbose=2)
        logging.info(f"Success training subject {subj}")
    except:
        logging.exception(f"Failed training subject {subj}")
        return False

    model.save(result, pars_file, subj=subj, study=study)

    _plot_model_fit(subj, data, dir_results, dir_results, model, scans=scans)

    return True


def plot_one_scan(subj, data, dir_results, dir_plots):
    config = CONFIG
    data |= INIT 
    model = dc.InverseAortaLiverDrug(**config)

    return _plot_model_fit(subj, data, dir_results, dir_plots, model, scans=1)


def plot_two_scan(subj, data, dir_results, dir_plots):
    config = CONFIG | {
        'bolus': 'dual',               # Two independent injections performed
        'non_stationary': 'U',         # Uptake rate changes during the course of the experiment
    }
    data |= INIT
    model = dc.InverseAortaLiverSplitDrug(**config)

    return _plot_model_fit(subj, data, dir_results, dir_plots, model,  scans=2)


def _plot_model_fit(subj, data, dir_results, dir_plots, model, scans=None):

    pars_path = os.path.join(dir_results, 'Pars')
    pars_file = os.path.join(pars_path, subj + '.dmr')
    if not os.path.exists(pars_file + '.zip'):
        return

    # Prepare output files
    plot_path = os.path.join(dir_plots,  'Plots')
    pcov_path = os.path.join(dir_plots, 'Cov')
    pvar_path = os.path.join(dir_plots, 'Sdev')
    os.makedirs(plot_path, exist_ok=True)
    os.makedirs(pcov_path, exist_ok=True)
    os.makedirs(pvar_path, exist_ok=True)

    plot_file_all = os.path.join(plot_path, subj + '.png')
    plot_file_init = os.path.join(plot_path, subj + '_first_pass.png')
    plot_file_2nd = os.path.join(plot_path, subj + '_second_pass.png')
    plot_file_pcov = os.path.join(pcov_path, subj + '.png')
    plot_file_pvar = os.path.join(pvar_path, subj + '.png')

    logging.info(f"Plotting subject {subj}")

    try:
        result, _, _ = model.load(pars_file)
    except:
        logging.exception(f"Failed plotting subject {subj}")
        return False

    model.plot(data, result, fname=plot_file_all, show=False)
    
    if scans==2:
        xlim = {1: [0, 5 * 60], 2: [0, 5 * 60]}
        model.plot(data, result, fname=plot_file_init, show=False, xlim=xlim)

        x1, x2 = data['tS_2_ao'][0], data['tS_4_ao'][0]
        xlim = {1: [x1, x1 + 15 * 60], 2: [x2, x2 + 15 * 60]}
        model.plot(data, result, fname=plot_file_2nd, show=False, xlim=xlim)
    else:
        model.plot(data, result, fname=plot_file_init, show=False, xlim=[0, 5 * 60])

    model.plot_pcov(result, plot_file_pcov, show=False)
    model.plot_pvar(result, plot_file_pvar, show=False, xlim=[-3, 3])

    logging.info(f"Success plotting subject {subj}")
    return True