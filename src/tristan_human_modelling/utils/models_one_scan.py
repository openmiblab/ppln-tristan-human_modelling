import os
import numpy as np
import pydmr


#from dcmri import AortaLiverDrug
from .aorta_liver_drug import AortaLiverDrug


def train(state, datafile, state_init='', staged=None):
    subj = os.path.basename(state).split('.')[0]
    data = pydmr.read(datafile, format='nest')
    model = model_init(data, subj, state_init)
    model_train(model, data, subj, verbose=2)
    model.save(state)


def save_plots(state, datafile, resultspath):

    if not os.path.exists(resultspath):
        return    
    plotpath = os.path.join(resultspath, 'Plots')
    if not os.path.exists(plotpath):
        os.makedirs(plotpath)

    subj = os.path.basename(state).split('.')[0]
    name = subj
    file = os.path.join(plotpath, name)
    if os.path.exists(file + '.png'):
        return

    data = pydmr.read(datafile, format='nest')
    rois = data['rois'][subj]
    model = AortaLiverDrug().load(state)
    xdata, ydata = model_data(rois)
    file = os.path.join(plotpath, subj + '.png')
    model.plot(xdata, ydata, fname=file, show=False)


def save_results(state, datafile, dmrpath, drug):

    subj = os.path.basename(state).split('.')[0]
    file = os.path.join(dmrpath, subj + '.dmr')
    if os.path.exists(file):
        return

    data = pydmr.read(datafile, format='nest')
    # rois = data['rois'][subj]
    pars = data['pars'][subj]

    # xdata, ydata = model_data(rois)
    # tb, Sb, tl, Sl = xdata[0], ydata[0], xdata[1], ydata[1]
    
    model = AortaLiverDrug().load(state)
    
    # Take this out for now until AUCs redefined
    # model.set_params(dose2 = [0,0])
    #params = export_params(model, tb, Sb, tl, Sl)

    params = model.export_params()
    params['c_T1_1'] = {'name': 'Control visit - liver T1-MOLLI at baseline', 'value': pars['control']['T1_liver_1'], 'unit': 'sec', 'sdev': 0}
    params['c_T1_2'] = {'name': 'Control visit - liver T1-MOLLI at 45min', 'value': pars['control']['T1_liver_2'], 'unit': 'sec', 'sdev': 0}
    params['c_T1_3'] = {'name': 'Control visit - liver T1-MOLLI at scan 2', 'value': pars['control']['T1_liver_3'], 'unit': 'sec', 'sdev': 0}
    params['d_T1_1'] = {'name': 'Drug visit - liver T1-MOLLI at baseline', 'value': pars['drug']['T1_liver_1'], 'unit': 'sec', 'sdev': 0}
    params['d_T1_2'] = {'name': 'Drug visit - liver T1-MOLLI at 45min', 'value': pars['drug']['T1_liver_2'], 'unit': 'sec', 'sdev': 0}
    params['d_T1_3'] = {'name': 'Drug visit - liver T1-MOLLI at scan 2', 'value': pars['drug']['T1_liver_3'], 'unit': 'sec', 'sdev': 0}

    # Exclude uninteresting
    excl = [
        'c_FA', 'd_FA', 'H', 'TS', 
        'c_B1corr_a', 'd_B1corr_a', 'c_B1corr_l', 'd_B1corr_l',
        'c_R20s_a', 'c_R20s_l', 'd_R20s_a', 'd_R20s_l', 
        'c_dose', 'd_dose', 
        'd_tmax', 'c_tmax',
        'dose_tolerance', 'dt', 'field_strength', 'rate'
    ]
    params = {k:v for k, v in params.items() if k not in excl}
    to_dmr(file, subj, params, drug)
    return file



def export_params(model:AortaLiverDrug, tb, Sb, tl, Sl):

    # TODO: Redo AUC for control and drug visits separately

    # Compute AUC over 3hrs
    BAT = model.params()['c_BAT']
    t = model.time()
    C = model.conc()
    R1, R2s = model.relax()
    tAUCb = (BAT < t['ctrl', 'aorta']) & (t['ctrl', 'aorta'] < BAT + 180 * 60)
    tAUCl = (BAT < t['ctrl', 'liver']) & (t['ctrl', 'liver'] < BAT + 180 * 60)
    AUC_Cb = np.trapezoid(C['ctrl', 'aorta'][tAUCb], t['ctrl', 'aorta'][tAUCb]) 
    AUC_Cl = np.trapezoid(C['ctrl', 'liver'][tAUCl], t['ctrl', 'liver'][tAUCl])

    # Compute relative enhancement at 20mins
    tRE = BAT + 20*60
    R1b, tb = R1['ctrl', 'aorta'], t['ctrl', 'aorta']
    R1l, tl = R1['ctrl', 'liver'], t['ctrl', 'liver']
    RE_R1b = (R1b[tb < tRE][-1] - R1b[0])/R1b[0]
    RE_R1l = (R1l[tl < tRE][-1] - R1l[0])/R1l[0]
    S0b = np.mean(Sb[tb < BAT - 30])
    S0l = np.mean(Sl[tl < BAT - 30])
    RE_Sb = (Sb[tb < tRE][-1] - S0b)/S0b
    RE_Sl = (Sl[tl < tRE][-1] - S0l)/S0l

    # Compute AUC over 35min
    tAUCb = (BAT < t['ctrl', 'aorta']) & (t['ctrl', 'aorta'] < BAT + 35 * 60)
    tAUCl = (BAT < t['ctrl', 'liver']) & (t['ctrl', 'liver'] < BAT + 35 * 60)
    AUC35_Cb = np.trapezoid(C['ctrl', 'aorta'][tAUCb], t['ctrl', 'aorta'][tAUCb]) 
    AUC35_Cl = np.trapezoid(C['ctrl', 'liver'][tAUCl], t['ctrl', 'liver'][tAUCl])

    pars = model.export_params()
    pars['c_AUC_Cb']={'name': 'Control visit AUC for Cb (0-inf)', 'value': AUC_Cb, 'unit': 'M*sec', 'sdev': 0}
    pars['c_AUC_Cl']={'name': 'Control visit AUC for Cl (0-inf)', 'value': AUC_Cl, 'unit': 'M*sec', 'sdev': 0} 
    pars['c_AUC35_Cb']={'name': 'Control visit AUC for Cb (0-35min)', 'value': AUC35_Cb, 'unit': 'M*sec', 'sdev': 0}
    pars['c_AUC35_Cl']={'name': 'Control visit AUC for Cl (0-35min)', 'value': AUC35_Cl, 'unit': 'M*sec', 'sdev': 0} 
    pars['c_RE_R1b']={'name': 'Control visit RE for R1b at 20min', 'value': RE_R1b, 'unit': '', 'sdev': 0}
    pars['c_RE_R1l']={'name': 'Control visit RE for R1l at 20min', 'value': RE_R1l, 'unit': '', 'sdev': 0}
    pars['c_RE_Sb']={'name': 'Control visit RE for Sb at 20min', 'value': RE_Sb, 'unit': '', 'sdev': 0}
    pars['c_RE_Sl']={'name': 'Control visit RE for Sl at 20min', 'value': RE_Sl, 'unit': '', 'sdev': 0}    

    return pars 
 

def to_dmr(file, subj, pars, drug):

    # Build data dictionary
    dmr = {
        'data': {},
        'pars': {},   
        'sdev': {},
        'columns': ['group'],
        #'columns': ['group', 'label'],
    }
    study = drug
    for key, val in pars.items():
        if isinstance(val['value'], str):
            continue
        group = 'MRI - aorta' if key in AORTA_PARS else 'MRI - liver'
        dmr['data'][key] = [val['name'], val['unit'], 'float', group]
        dmr['pars'][subj, study, key] = val['value']
        dmr['sdev'][subj, study, key] = 0 if val['sdev'] is None else val['sdev']

    pydmr.write(file, dmr)


AORTA_PARS = [
    'c_RE_Sb', 'c_RE_R1b', 'c_S02_a', 'd_S02_a', 'c_BAT', 'd_BAT', 
    'GFR', 'CO', 'Thl', 'Dhl',
    'To', 'Eo', 'To_e', 'c_BAT2', 'd_BAT2', 'GFR',
    'c_AUC_R1b', 'c_AUC_Cb', 'c_AUC35_R1b', 'c_AUC35_Cb',
]


def model_data(rois):

    time_1_control = rois['control']['time_1'] - rois['control']['time_1'][0]
    time_1_drug = rois['drug']['time_1'] - rois['drug']['time_1'][0]

    xdata = (
        time_1_control[rois['control']['aorta_1_accept']], 
        time_1_control[rois['control']['liver_1_accept']],
        time_1_drug[rois['drug']['aorta_1_accept']], 
        time_1_drug[rois['drug']['liver_1_accept']],
    )
    ydata = (
        rois['control']['aorta_1'][rois['control']['aorta_1_accept']],
        rois['control']['liver_1'][rois['control']['liver_1_accept']],
        rois['drug']['aorta_1'][rois['drug']['aorta_1_accept']],
        rois['drug']['liver_1'][rois['drug']['liver_1_accept']],
    )
    return xdata, ydata


def model_init(data, subj, state_init):

    if os.path.exists(state_init):
        return AortaLiverDrug().load(state_init)

    rois = data['rois'][subj]
    pars = data['pars'][subj]

    # Get time arrays
    time_1_control = rois['control']['time_1'] - rois['control']['time_1'][0]
    time_1_drug = rois['drug']['time_1'] - rois['drug']['time_1'][0]

    if 'eGFR_abs' in pars['screening']:
        c_eGFR = pars['screening']['eGFR_abs']
        d_eGFR = pars['screening']['eGFR_abs']
    else:
        c_eGFR = pars['control']['eGFR_abs']
        d_eGFR = pars['drug']['eGFR_abs']
    eGFR = (c_eGFR + d_eGFR) / 2

    pars = { 

        # Simulation parameters
        'dt': 0.1, 
        'c_tmax': max(time_1_control),
        'd_tmax': max(time_1_drug),
        'dose_tolerance': 0.1, # no effect moving to 0.01
        'field_strength': 3.0,

        # Injection parameters
        'weight': (pars['control']['weight'] + pars['drug']['weight']) / 2,
        'agent': 'gadoxetate',
        'c_dose': pars['control']['dose_1'], 
        'd_dose': pars['drug']['dose_1'],
        'rate': 1.0,

        # Acquisition parameters
        'TR': pars['control']['TR'], 
        'c_FA': pars['control']['FA_1'],
        'd_FA': pars['drug']['FA_1'],
        'TS': np.min(time_1_control[1:] - time_1_control[:-1]), 
        
        # Signal parameters
        'c_R10_a': 1/pars['control']['T1_aorta_1'], 
        'd_R10_a': 1/pars['drug']['T1_aorta_1'],
        'c_R10_l': 1/pars['control']['T1_liver_1'], 
        'd_R10_l': 1/pars['drug']['T1_liver_1'],

        # Tissue parameters
        'H': 0.45,
        'c_vol': pars['control']['liver_volume'], 
        'd_vol': pars['drug']['liver_volume'],

        # Aorta parameters
        'c_BAT': 60, # irrelevant take out
        'd_BAT': 60, # irrelevant
        'CO': 100,
        'Thl': 10,
        'Dhl': 0.2,
        'To': 20,
        'Eo': 0.15,
        'To_e': 120,
        'GFR': eGFR / 60,

        # Liver parameters
        've': 0.15,
        'Tg': 30,
        'Dg': 0.85, # reducing makes no difference
        'c_khe': 0.0025,
        'd_khe': 0.0025,
        'c_kbh': 0.00025,
        'd_kbh': 0.00025,
    }

    return AortaLiverDrug(**pars) 


def model_train(model:AortaLiverDrug, data, subj, verbose=0):

    rois = data['rois'][subj]
    pars = data['pars'][subj]

    # Get data at valid (accepted) time points
    xdata, ydata = model_data(rois)

    # Train the model to the data
    model.train(
        xdata, 
        ydata, 
        n0=[
            max([np.sum(xdata[0] < pars['control']['t0']), 2]),
            max([np.sum(xdata[0] < pars['drug']['t0']), 2]),
        ], 
        xtol=1e-3, # no improvement moving this to 1e-6
        # max_nfev=1, # for debugging
        verbose=verbose,
        free = {
            'c_BAT': [-30, 30],
            'd_BAT': [-30, 30],

            # 'GFR': [0.5, 3],
            'CO': [0, 300],
            'Thl': [0, 30],
            'Dhl': [0.05, 0.95],
            'To': [0, 60],
            'Eo': [0, 0.5],
            'To_e': [0, 800], # 2 hit the ceiling
            # 'To_e': [0, 1200], # this keeps it away from the ceiling but little effect on outcomes
            've': [0.01, 0.3],

            # THIS IS CRITICAL!!!! Avoid hitting the ceiling.
            #'Tg': [15, 60], # Original setting
            'Tg': [15, 180], # Improved setting
            'Dg': [0, 1],

            'c_khe': [0, 0.006],
            'd_khe': [0, 0.006],
            'c_kbh': [0, 0.0006],
            'd_kbh': [0, 0.0006],
        },

    )


