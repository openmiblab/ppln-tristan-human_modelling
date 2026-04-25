import numpy as np
import dcmri as dc


def two_scan_data(rois):
    xdata = (
        rois['time_1'][rois['aorta_1_accept']] - rois['time_1'][0], 
        rois['time_2'][rois['aorta_2_accept']] - rois['time_1'][0], 
        rois['time_1'][rois['liver_1_accept']] - rois['time_1'][0],
        rois['time_2'][rois['liver_2_accept']] - rois['time_1'][0],
    )
    ydata = (
        rois['aorta_1'][rois['aorta_1_accept']],
        rois['aorta_2'][rois['aorta_2_accept']],
        rois['liver_1'][rois['liver_1_accept']],
        rois['liver_2'][rois['liver_2_accept']],
    )
    return xdata, ydata


def two_scan(data, subj, visit, verbose=0):

    rois = data['rois'][subj][visit]
    pars = data['pars'][subj][visit]

    time_1 = rois['time_1'] - rois['time_1'][0]
    time_2 = rois['time_2'] - rois['time_1'][0]

    # Define default model
    model = dc.AortaLiver2scan(

        # Configuration
        kinetics='1I-IC-HFD', 
        non_stationary='UE', 
        sequence='SS', 

        # Parameters
        dt = 0.1,
        tmax = max(time_2),
        dose_tolerance = 0.1,
        field_strength = 3.0,

        # Injection parameters
        weight=pars['weight'],
        agent='gadoxetate',
        dose=pars['dose_1'],
        rate=1.0,
        dose2=pars['dose_2'],
        t_scan2 = pars['T1_time_3'] - rois['time_1'][0] - 120,

        # Acquisition parameters
        TR=pars['TR'], 
        FA=pars['FA_1'],
        FA2=pars['FA_2'],
        TS = np.min(time_1[1:] - time_1[:-1]),
        
        # Signal parameters
        R10a=1/pars['T1_aorta_1'],
        R10l=1/pars['T1_liver_1'],

        # Tissue parameters
        H=0.45,
        vol=pars['liver_volume'],

        # Aorta parameters
        BAT2=1200,
        BAT=60,
        CO=100,
        Thl=10,
        Dhl=0.2,
        To=20,
        Eo=0.15,
        Toe=120,
        Eb=0.05,

        # Liver parameters
        ve=0.15,
        Tg=30,
        Dg=0.85,
        khe_i=0.0025,
        khe_f=0.0025,
        Th_i=30 * 60,
        Th_f=30 * 60,
    )

    # Personalise model
    xdata, ydata = two_scan_data(rois)
    model.train(
        xdata, 
        ydata, 
        n0=max([np.sum(xdata[0] < pars['t0']), 2]), 
        R102a=1/pars['T1_aorta_3'],
        R102l=1/pars['T1_liver_3'],
        xtol=1e-3, 
        # max_nfev=1,
        joint=False,
        verbose=verbose,
        free = {
            'S02a': [0.01, 100],
            'S02l': [0.01, 100],

            'BAT2': [-30, 30],
            'BAT': [-30, 30],
            'CO': [0, 300],
            'Thl': [0, 30],
            'Dhl': [0.05, 0.95],
            'To': [0, 60],
            'Eo': [0, 0.5],
            'Toe': [0, 800],
            'Eb': [0.01, 0.15],

            've': [0.01, 0.3],
            'Tg': [15, 60], 
            'Dg': [0, 1],
            'khe_i': [0, 0.01], 
            'khe_f': [0, 0.01],
            'Th_i': [10 * 60, 600 * 60],
            'Th_f': [10 * 60, 600 * 60],
        }
    )
    return model




def one_scan_data(rois):

    xdata = (
        rois['time_1'][rois['aorta_1_accept']] - rois['time_1'][0], 
        rois['time_1'][rois['liver_1_accept']] - rois['time_1'][0],
    )
    ydata = (
        rois['aorta_1'][rois['aorta_1_accept']], 
        rois['liver_1'][rois['liver_1_accept']],
    )
    return xdata, ydata


def one_scan(data, subj, visit, verbose=0, tacq=None):

    rois = data['rois'][subj][visit]
    pars = data['pars'][subj][visit]

    time_1 = rois['time_1'] - rois['time_1'][0]

    # Define default model
    model = dc.AortaLiver(

        # Configuration
        kinetics='1I-IC-HFD', 
        non_stationary=None, 
        sequence='SS', 

        # Parameters
        dt = 0.5,
        tmax = max(time_1),
        dose_tolerance = 0.1,
        field_strength = 3.0,

        # Injection parameters
        weight=pars['weight'],
        agent='gadoxetate',
        dose=pars['dose_1'],
        rate=1.0,

        # Acquisition parameters
        TR=pars['TR'], 
        FA=pars['FA_1'],
        TS = np.min(time_1[1:] - time_1[:-1]),
        
        # Signal parameters
        R10a=1/pars['T1_aorta_1'],
        R10l=1/pars['T1_liver_1'],

        # Tissue parameters
        H=0.45,
        vol=pars['liver_volume'],

        # Aorta parameters
        BAT=60,
        CO=100,
        Thl=10,
        Dhl=0.2,
        To=20,
        Eo=0.15,
        Toe=120,
        Eb=0.05,

        # Liver parameters
        ve=0.3,
        Tg=30,
        Dg=0.65,
        khe_i=0.0025,
        khe_f=0.0025,
        Th_i=30*60,
        Th_f=30*60,
    )

    # Personalise model

    # Truncate data if requested
    xdata, ydata = one_scan_data(rois)
    if tacq is not None:
        idx0, idx1 = xdata[0]<tacq*60, xdata[1]<tacq*60
        xdata = (xdata[0][idx0], xdata[1][idx1])
        ydata = (ydata[0][idx0], ydata[1][idx1])

    # Train
    model.train(
        xdata, 
        ydata, 
        n0=max([np.sum(xdata[0] < pars['t0']), 2]), 
        xtol=1e-3, 
        joint=True,
        verbose=verbose,
        free = {
            'BAT': [-30, 30],
            'CO': [0, 300],
            'Thl': [0, 30],
            'Dhl': [0.05, 0.95],
            'To': [0, 60],
            'Eo': [0, 0.5],
            'Toe': [0, 800],
            'Eb': [0.01, 0.15],

            've': [0.01, 0.6],
            'Tg': [15, 60], 
            'Dg': [0, 1],
            'khe': [0, 0.005],
            'Th': [20*60, 10*60*60],
        }
    )
    return model