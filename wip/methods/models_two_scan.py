import os
import numpy as np
import pydmr


from dcmri import AortaLiverDynamicDrug


def train(state, datafile, state_init):
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
    model = AortaLiverDynamicDrug().load(state)
    #model.set_params(tmax = [360 * 60, 360 * 60])
    xdata, ydata = two_scan_one_model_data(rois)
    file = os.path.join(plotpath, subj + '.png')
    model.plot(xdata, ydata, fname=file, show=False)


def save_results(state, datafile, resultspath):

    drug = os.path.basename(resultspath)
    dmrpath = os.path.join(resultspath, 'Results')
    subj = os.path.basename(state).split('.')[0]
    file = os.path.join(dmrpath, subj + '.dmr')
    if os.path.exists(file):
        return

    data = pydmr.read(datafile, format='nest')
    # rois = data['rois'][subj]
    pars = data['pars'][subj]

    # xdata, ydata = two_scan_one_model_data(rois)
    # tb, Sb, tl, Sl = xdata[0], ydata[0], xdata[2], ydata[2]
    
    model = AortaLiverDynamicDrug().load(state)
    
    # Take this out for now until AUCs redefined
    # model.set_params(dose_2 = [0,0])
    #params = export_params(model, tb, Sb, tl, Sl)

    params = model.export_params()
    params['c_T1_1'] = {'name': 'Control visit - liver T1-MOLLI at baseline', 'value': pars['control']['T1_liver_1'], 'unit': 'sec', 'sdev': 0}
    params['c_T1_2'] = {'name': 'Control visit - liver T1-MOLLI at 45min', 'value': pars['control']['T1_liver_2'], 'unit': 'sec', 'sdev': 0}
    params['c_T1_3'] = {'name': 'Control visit - liver T1-MOLLI at scan 2', 'value': pars['control']['T1_liver_3'], 'unit': 'sec', 'sdev': 0}
    params['d_T1_1'] = {'name': 'Drug visit - liver T1-MOLLI at baseline', 'value': pars['drug']['T1_liver_1'], 'unit': 'sec', 'sdev': 0}
    params['d_T1_2'] = {'name': 'Drug visit - liver T1-MOLLI at 45min', 'value': pars['drug']['T1_liver_2'], 'unit': 'sec', 'sdev': 0}
    params['d_T1_3'] = {'name': 'Drug visit - liver T1-MOLLI at scan 2', 'value': pars['drug']['T1_liver_3'], 'unit': 'sec', 'sdev': 0}
 
    to_dmr(file, subj, params, drug)
    return file



def export_params(model:AortaLiverDynamicDrug, tb, Sb, tl, Sl):

    # TODO: Redo AUC for control and drug visits separately

    # Compute AUC over 3hrs
    BAT = model.params('c_BAT_1')
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
    'RE_Sb', 'RE_R1b', 'c_S0_2_a', 'd_S0_2_a' 'c_BAT_1', 'd_BAT_1',
    # 'GFR', 
    'Eb',
    'CO','Thl','Dhl',
    'To', 'Eo', 'To_e', 'c_BAT_2', 'd_BAT_2', 
    'AUC_R1b','AUC_Cb', 'AUC35_R1b','AUC35_Cb',
]



def two_scan_one_model_data(rois):

    # Get time arrays
    time_1_control = rois['control']['time_1'] - rois['control']['time_1'][0]
    time_2_control = rois['control']['time_2'] - rois['control']['time_1'][0]
    time_1_drug = rois['drug']['time_1'] - rois['drug']['time_1'][0]
    time_2_drug = rois['drug']['time_2'] - rois['drug']['time_1'][0]

    xdata = (
        time_1_control[rois['control']['aorta_1_accept']], 
        time_2_control[rois['control']['aorta_2_accept']], 
        time_1_control[rois['control']['liver_1_accept']],
        time_2_control[rois['control']['liver_2_accept']],
        time_1_drug[rois['drug']['aorta_1_accept']], 
        time_2_drug[rois['drug']['aorta_2_accept']], 
        time_1_drug[rois['drug']['liver_1_accept']],
        time_2_drug[rois['drug']['liver_2_accept']],
    )
    ydata = (
        rois['control']['aorta_1'][rois['control']['aorta_1_accept']],
        rois['control']['aorta_2'][rois['control']['aorta_2_accept']],
        rois['control']['liver_1'][rois['control']['liver_1_accept']],
        rois['control']['liver_2'][rois['control']['liver_2_accept']],
        rois['drug']['aorta_1'][rois['drug']['aorta_1_accept']],
        rois['drug']['aorta_2'][rois['drug']['aorta_2_accept']],
        rois['drug']['liver_1'][rois['drug']['liver_1_accept']],
        rois['drug']['liver_2'][rois['drug']['liver_2_accept']],
    )
    return xdata, ydata


def two_scan_sigma(time, baseline=540, weight=0.01):
    sigma = tuple([np.ones_like(t) for t in time])
    for i in [1, 3, 5, 7]:
        sigma[i][time[i] < time[i][0] + baseline] = weight
    return sigma


def model_init(data, subj, state_init):

    if os.path.exists(state_init):
        return AortaLiverDynamicDrug().load(state_init)

    rois = data['rois'][subj]
    pars = data['pars'][subj]

    # Get time arrays
    time_1_control = rois['control']['time_1'] - rois['control']['time_1'][0]
    time_2_control = rois['control']['time_2'] - rois['control']['time_1'][0]
    time_1_drug = rois['drug']['time_1'] - rois['drug']['time_1'][0]
    time_2_drug = rois['drug']['time_2'] - rois['drug']['time_1'][0]

    pars = { 

        # Simulation parameters
        'dt': 0.1, 
        'c_tmax': max(time_2_control),
        'd_tmax_1': max(time_1_drug),
        'd_tmax': max(time_2_drug),
        'dose_tolerance': 0.1,
        'field_strength': 3.0,

        # Injection parameters
        'weight': (pars['control']['weight'] + pars['drug']['weight']) / 2,
        'agent': 'gadoxetate',
        'c_dose_1': pars['control']['dose_1'], 
        'd_dose_1': pars['drug']['dose_1'],
        'rate': 1.0,
        'c_dose_2':pars['control']['dose_2'], 
        'd_dose_2': pars['drug']['dose_2'],
        'c_t_scan2': pars['control']['T1_time_3'] - rois['control']['time_1'][0] - 120, 
        'd_t_scan2': pars['drug']['T1_time_3'] - rois['drug']['time_1'][0] - 120,

        # Acquisition parameters
        'TR': pars['control']['TR'], 
        'FA': pars['control']['FA_1'],
        'FA2': pars['control']['FA_2'],
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
        'c_BAT_2': 1200,
        'd_BAT_2': 1200,
        'c_BAT_1': 60,
        'd_BAT_1': 60,
        'CO': 100,
        'Thl': 10,
        'Dhl': 0.2,
        'To': 20,
        'Eo': 0.15,
        'To_e': 120,
        # 'GFR': pars['screening']['eGFR_abs'] / 60,
        'Eb': 0.05,

        # Liver parameters
        've': 0.15,
        'Tg': 30,
        'Dg': 0.85,
        'c_khe_i': 0.0025,
        'c_khe_f': 0.0025,
        'd_khe_i': 0.0025,
        'd_khe_f': 0.0025,
        'c_kbh': 0.00025,
        'd_kbh': 0.00025,
    }

    return AortaLiverDynamicDrug(**pars) 


def model_train(model: AortaLiverDynamicDrug, data, subj, verbose=0):

    rois = data['rois'][subj]
    pars = data['pars'][subj]

    # Get data at valid (accepted) time points
    xdata, ydata = two_scan_one_model_data(rois)
    # sigma = two_scan_sigma(xdata, baseline=540, weight=0.01)

    # Train the model to the data
    model.train(
        xdata, 
        ydata, 
        n0=[
            max([np.sum(xdata[0] < pars['control']['t0']), 2]),
            max([np.sum(xdata[0] < pars['drug']['t0']), 2]),
        ], 
        R102a=[
            1/pars['control']['T1_aorta_3'], 
            1/pars['drug']['T1_aorta_3'],
        ],
        R102l=[
            1/pars['control']['T1_liver_3'], 
            1/pars['drug']['T1_liver_3'],
        ],
        xtol=1e-3,
        # max_nfev=1,
        # sigma=sigma,
        verbose=verbose,
        free = {
            
            'c_S0_2_a': [0.01, 100],
            'd_S0_2_a': [0.01, 100],
            'c_S0_2_l': [0.01, 100],
            'd_S0_2_l': [0.01, 100],
            
            'c_BAT_2': [-30, 30],
            'd_BAT_2': [-30, 30],
            'c_BAT_1': [-30, 30],
            'd_BAT_1': [-30, 30],

            'Eb': [0.01, 0.15],
            'CO': [0, 300],
            'Thl': [0, 30],
            'Dhl': [0.05, 0.95],
            'To': [0, 60],
            'Eo': [0, 0.5],
            'To_e': [0, 800],
            've': [0.01, 0.3],
            'Tg': [15, 180], 
            'Dg': [0, 1],

            'c_khe_i': [0, 0.006],
            'c_khe_f': [0, 0.006],
            'd_khe_i': [0, 0.006],
            'd_khe_f': [0, 0.006],
            'c_kbh': [0, 0.0006],
            'd_kbh': [0, 0.0006],
        },
    )


