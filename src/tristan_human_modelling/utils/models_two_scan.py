import os
import numpy as np
import pydmr

from .tables import export_covariance_matrix
from .figs import export_csv_to_heatmap
from dcmri import AortaLiverDynamicDrug

def train(state, datafile, pars_path, state_init=None, staged=0):
    subj = os.path.basename(state).split('.')[0]
    data = pydmr.read(datafile, format='nest')
    if state_init is not None:
        model = AortaLiverDynamicDrug().load(state_init)
    else:
        model = model_init(data, subj)
    vals, sdev, pcov = model_train(model, data, subj, verbose=2, staged=staged)
    model.save(state) 
    save_results(state, datafile, pars_path, sdev)

    pcov_csv = os.path.join(pars_path, subj + '_pcov.csv')
    pcov_png = os.path.join(pars_path, subj + '_pcov.png')
    export_covariance_matrix(pcov, list(vals.keys()), output_filename=pcov_csv)
    export_csv_to_heatmap(pcov_csv, pcov_png)


def save_plots(state, datafile, resultspath):   
    plotpath = os.path.join(resultspath, 'Plots')
    os.makedirs(plotpath, exist_ok=True)

    subj = os.path.basename(state).split('.')[0]
    name = subj
    file = os.path.join(plotpath, name)

    data = pydmr.read(datafile, format='nest')
    rois = data['rois'][subj]
    model = AortaLiverDynamicDrug().load(state)
    #model.set_params(tmax = [360 * 60, 360 * 60])
    xdata, ydata = model_data(rois)
    file = os.path.join(plotpath, subj + '.png')
    model.plot(xdata, ydata, clim=[1.0 / 1000, 0.2 / 1000], fname=file, show=False)


def save_results(state, datafile, dmrpath, sdev):
    subj = os.path.basename(state).split('.')[0]
    data = pydmr.read(datafile, format='nest')
    pars = data['pars'][subj]

    model = AortaLiverDynamicDrug().load(state)
    model._pars['c_dose_2'] = 0  # AUCs are computed on a single dose
    model._pars['d_dose_2'] = 0

    params = model.export_params(sdev, desc=True)
    params['c_T1_1'] = {'name': 'Control visit - liver T1-MOLLI at baseline', 'value': pars['control']['T1_liver_1'], 'unit': 'sec', 'sdev': 0}
    params['c_T1_2'] = {'name': 'Control visit - liver T1-MOLLI at 45min', 'value': pars['control']['T1_liver_2'], 'unit': 'sec', 'sdev': 0}
    params['c_T1_3'] = {'name': 'Control visit - liver T1-MOLLI at scan 2', 'value': pars['control']['T1_liver_3'], 'unit': 'sec', 'sdev': 0}
    params['d_T1_1'] = {'name': 'Drug visit - liver T1-MOLLI at baseline', 'value': pars['drug']['T1_liver_1'], 'unit': 'sec', 'sdev': 0}
    params['d_T1_2'] = {'name': 'Drug visit - liver T1-MOLLI at 45min', 'value': pars['drug']['T1_liver_2'], 'unit': 'sec', 'sdev': 0}
    params['d_T1_3'] = {'name': 'Drug visit - liver T1-MOLLI at scan 2', 'value': pars['drug']['T1_liver_3'], 'unit': 'sec', 'sdev': 0}
 
    # Exclude uninteresting
    excl = [
        'c_FA_1', 'c_FA_2', 'd_FA_1', 'd_FA_2', 'H', 'TR', 'TS', 
        'c_B1corr_1_a', 'd_B1corr_1_a', 'c_B1corr_1_l', 'd_B1corr_1_l',
        'c_B1corr_2_a', 'd_B1corr_2_a', 'c_B1corr_2_l', 'd_B1corr_2_l',
        'R20s_a', 'R20s_l', 
        'c_dose_1', 'c_dose_2', 'd_dose_1', 'd_dose_2', 
        'd_t_scan2', 'd_tmax', 'c_t_scan2', 'c_tmax',
        'dose_tolerance', 'dt', 'field_strength', 'rate'
    ]

    # to dmr
    file = os.path.join(dmrpath, subj + '.dmr')
    study = data['pars'][subj]['screening']['Study']
    dmr = {
        'data': {},
        'pars': {},   
        'sdev': {},
    }
    for key, val in params.items():
        if key in excl:
            continue
        if isinstance(val['value'], str):
            continue
        dmr['data'][key] = [val['name'], val['unit'], 'float']
        dmr['pars'][subj, study, key] = val['value']
        if val['sdev'] is not None:
            dmr['sdev'][subj, study, key] = val['sdev']
    pydmr.write(file, dmr)



def model_data(rois):

    # Get time arrays
    time_1_control = rois['control']['time_1'] - rois['control']['time_1'][0]
    time_2_control = rois['control']['time_2'] - rois['control']['time_1'][0]
    time_1_drug = rois['drug']['time_1'] - rois['drug']['time_1'][0]
    time_2_drug = rois['drug']['time_2'] - rois['drug']['time_1'][0]

    xdata = {
        ('ctrl', 'aorta', 1): time_1_control[rois['control']['aorta_1_accept']], 
        ('ctrl', 'aorta', 2): time_2_control[rois['control']['aorta_2_accept']], 
        ('ctrl', 'liver', 1): time_1_control[rois['control']['liver_1_accept']],
        ('ctrl', 'liver', 2): time_2_control[rois['control']['liver_2_accept']],
        ('drug', 'aorta', 1): time_1_drug[rois['drug']['aorta_1_accept']], 
        ('drug', 'aorta', 2): time_2_drug[rois['drug']['aorta_2_accept']], 
        ('drug', 'liver', 1): time_1_drug[rois['drug']['liver_1_accept']],
        ('drug', 'liver', 2): time_2_drug[rois['drug']['liver_2_accept']],
    }
    ydata = {
        ('ctrl', 'aorta', 1): rois['control']['aorta_1'][rois['control']['aorta_1_accept']],
        ('ctrl', 'aorta', 2): rois['control']['aorta_2'][rois['control']['aorta_2_accept']],
        ('ctrl', 'liver', 1): rois['control']['liver_1'][rois['control']['liver_1_accept']],
        ('ctrl', 'liver', 2): rois['control']['liver_2'][rois['control']['liver_2_accept']],
        ('drug', 'aorta', 1): rois['drug']['aorta_1'][rois['drug']['aorta_1_accept']],
        ('drug', 'aorta', 2): rois['drug']['aorta_2'][rois['drug']['aorta_2_accept']],
        ('drug', 'liver', 1): rois['drug']['liver_1'][rois['drug']['liver_1_accept']],
        ('drug', 'liver', 2): rois['drug']['liver_2'][rois['drug']['liver_2_accept']],
    }
    return xdata, ydata


def model_sigma(time):
    sigma = tuple([np.ones_like(t) for t in time.values()])

    # Prioritise fit to liver curve over aorta

    # weight x 8 for liver signals in the first scan
    # To improve the fit to the first uptake
    for i in [2, 6]:
        sigma[i][:] /= 8

    # weight x 4 for liver signals of the second scan
    for i in [3, 7]:
        sigma[i][:] /= 4

    return sigma


def model_init(data, subj):

    rois = data['rois'][subj]
    pars = data['pars'][subj]

    # Get time arrays
    time_1_control = rois['control']['time_1'] - rois['control']['time_1'][0]
    time_2_control = rois['control']['time_2'] - rois['control']['time_1'][0]
    time_1_drug = rois['drug']['time_1'] - rois['drug']['time_1'][0]
    time_2_drug = rois['drug']['time_2'] - rois['drug']['time_1'][0]

    if 'eGFR_abs' in pars['screening']:
        c_eGFR = pars['screening']['eGFR_abs']
        d_eGFR = pars['screening']['eGFR_abs']
    else:
        c_eGFR = pars['control']['eGFR_abs']
        d_eGFR = pars['drug']['eGFR_abs']
    eGFR = (c_eGFR + d_eGFR) / 2

    pars = { 

        # Simulation parameters
        'dt': 0.5, 
        #'dt': 0.1, 
        'c_tmax': max(time_2_control),
        #'d_tmax_1': max(time_1_drug), # not needed
        'd_tmax': max(time_2_drug),
        'dose_tolerance': 0.01,
        'field_strength': 3.0,

        # Injection parameters
        'weight': (pars['control']['weight'] + pars['drug']['weight']) / 2,
        'agent': 'gadoxetate',
        'c_dose_1': pars['control']['dose_1'], 
        'd_dose_1': pars['drug']['dose_1'],
        'rate': 1.0,
        'c_dose_2':pars['control']['dose_2'], 
        'd_dose_2': pars['drug']['dose_2'],
        # 'c_t_scan2': pars['control']['T1_time_3'] - rois['control']['time_1'][0] - 120, 
        # 'd_t_scan2': pars['drug']['T1_time_3'] - rois['drug']['time_1'][0] - 120,
        'c_t_scan2': time_2_control[0], 
        'd_t_scan2': time_2_drug[0],

        # Acquisition parameters
        'TR': pars['control']['TR'], 
        'c_FA_1': pars['control']['FA_1'],
        'c_FA_2': pars['control']['FA_2'],
        'd_FA_1': pars['drug']['FA_1'],
        'd_FA_2': pars['drug']['FA_2'],
        'TS': np.min(time_1_control[1:] - time_1_control[:-1]), 
        
        # Signal parameters
        'R10_a': np.mean([1/pars['control']['T1_aorta_1'], 1/pars['drug']['T1_aorta_1']]),
        'R10_l': np.mean([1/pars['control']['T1_liver_1'], 1/pars['drug']['T1_liver_1']]),

        # Tissue parameters
        'H': 0.45,
        'c_vol': pars['control']['liver_volume'], 
        'd_vol': pars['drug']['liver_volume'],

        # Fixed kinetics
        'FFl': 0.25,
        'GFR': eGFR / 60,
        'CO': 100,
        'Thl': 10,
        'Dhl': 0.2,
        'To': 20,
        'Eo': 0.15,
        'To_e': 120,
        'Tg': 30,
        've': 0.25,

        # Control kinetics
        'c_khe_i': 0.0025,
        'c_khe_f': 0.0025,
        'c_kbh': 0.00025,

        # Drug kinetics
        'd_khe_i': 0.0025,
        'd_khe_f': 0.0025,
        'd_kbh': 0.00025,

    }

    return AortaLiverDynamicDrug(**pars) 


def model_train(model: AortaLiverDynamicDrug, data, subj, verbose=0, staged=0):

    rois = data['rois'][subj]
    pars = data['pars'][subj]

    # Get data at valid (accepted) time points
    xdata, ydata = model_data(rois)

    n0 = [
        max([np.sum(xdata['ctrl', 'aorta', 1] < xdata['ctrl', 'aorta', 1][0] + pars['control']['t0']), 2]),
        max([np.sum(xdata['drug', 'aorta', 1] < xdata['drug', 'aorta', 1][0] + pars['drug']['t0']), 2]),
    ]

    # Train the model to the data
    return model.train(
        xdata, 
        ydata, 
        n0=n0,
        
        verbose=verbose, xtol=1e-6, n_runs=1, # max_nfev=1, # for curve_fit
        sigma=model_sigma(xdata),

        free = {
            'c_BAT_1': [-30, 30],
            'c_BAT_2': [-30, 30],
            'd_BAT_1': [-30, 30],
            'd_BAT_2': [-30, 30],

            'FFl': [0.01, 0.99],
            #'GFR': [0.5, 3],
            'CO': [30, 300],
            'Thl': [0, 30],
            'Dhl': [0.05, 0.95],
            'To': [0, 60],
            'Eo': [0, 0.5],
            'To_e': [0, 800], 
            'Tg': [15, 120],

            've': [0.01, 0.5],
            'c_khe_i': [0, 0.006],
            'c_khe_f': [0, 0.006],
            'c_kbh': [0, 0.0007],

            'd_khe_i': [0, 0.006],
            'd_khe_f': [0, 0.006],
            'd_kbh': [0, 0.0007],
        },
        
    )


