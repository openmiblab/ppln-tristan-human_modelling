import os
import time

import numpy as np
import dcmri as dc
import pydmr

from methods import tools, models



def compute(datafile, resultspath):
    train(datafile, resultspath)
    plot(datafile, resultspath)
    to_dmr(datafile, resultspath)



def train(datafile, resultspath):

    state_path = os.path.join(resultspath, 'State')
    if not os.path.exists(state_path):
        os.makedirs(state_path)

    start = time.time()

    data = pydmr.read(datafile, format='nest')
    for subj in data['rois'].keys():
        for visit in data['rois'][subj].keys():
            state = os.path.join(state_path, f"{subj}_{visit}.json")
            if os.path.exists(state):
                continue
            model = models.two_scan(data, subj, visit, verbose=2)
            model.save(state)

    print('Calculation time (mins): ', (time.time()-start)/60)



def plot(datafile, resultspath):
    state_path = os.path.join(resultspath, 'State')
    for f in os.listdir(state_path):
        state = os.path.join(state_path, f)
        save_plots(state, datafile, resultspath)



def to_dmr(datafile, resultspath):
    state_path = os.path.join(resultspath, 'State')
    results = []
    for f in os.listdir(state_path):
        state = os.path.join(state_path, f)
        file = save_results(state, datafile, resultspath)
        results.append(file)
    file = os.path.join(resultspath, 'all_results')
    pydmr.concat(results, file)




def save_plots(state, datafile, resultspath):

    subj, visit = os.path.basename(state).split('.')[0].split('_')

    plotpath = os.path.join(resultspath, 'Plots')
    if not os.path.exists(plotpath):
        os.makedirs(plotpath)
    name = subj + '_' + visit
    file = os.path.join(plotpath, name)
    if os.path.exists(file + '.png'):
        return

    data = pydmr.read(datafile, format='nest')
    rois = data['rois'][subj][visit]
    pars = data['pars'][subj][visit]

    xdata, ydata = models.two_scan_data(rois)
    t0 = rois['time_1'][0]

    t = [
        0, 
        #pars['T1_time_2']-t0,
        pars['T1_time_3']-t0,
    ]
    R1a = [
        1/pars['T1_aorta_1'], 
        #1/pars['T1_aorta_2'],
        1/pars['T1_aorta_3'],
    ]
    R1l = [
        1/pars['T1_liver_1'], 
        #1/pars['T1_liver_2'],
        1/pars['T1_liver_3'],
    ]

    model = dc.AortaLiver2scan().load(state)
    ya = [dc.signal_ss(model.params('S0a'), R1a[0], model.params('TR'), model.params('FA')),
          #dc.signal_ss(model.pars['S0a'], R1a[1], model.pars['TR'], model.pars['FA']),
          dc.signal_ss(model.params('S02a'), R1a[1], model.params('TR'), model.params('FA'))]
    yl = [dc.signal_ss(model.params('S0l'), R1l[0], model.params('TR'), model.params('FA')),
          #dc.signal_ss(model.pars['S0l'], R1l[1], model.pars['TR'], model.pars['FA']),
          dc.signal_ss(model.params('S02l'), R1l[1], model.params('TR'), model.params('FA'))]
    test = ((t,ya),(t,yl))
    model.plot(xdata, ydata, fname=file + '.png', ref=test, show=False)

    BAT = model.params('BAT')
    model.plot(xdata, ydata, xlim=[xdata[0][0], xdata[0][-1]], 
               fname=file + '_scan1_win1.png', ref=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+600], 
               fname=file + '_scan1_win2.png', ref=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+160], 
               fname=file + '_scan1_win3.png', ref=test, show=False)
    
    BAT = model.params('BAT2')
    model.plot(xdata, ydata, xlim=[xdata[3][0], xdata[3][-1]], 
               fname=file + '_scan2_win1.png', ref=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+600], 
               fname=file + '_scan2_win2.png', ref=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+160], 
               fname=file + '_scan2_win3.png', ref=test, show=False)


def save_results(state, datafile, resultspath):

    subj, visit = os.path.basename(state).split('.')[0].split('_')

    data = pydmr.read(datafile, format='nest')
    rois = data['rois'][subj][visit]
    pars = data['pars'][subj][visit]

    xdata, ydata = models.two_scan_data(rois)
    t0 = rois['time_1'][0]
    tb, Sb, tl, Sl = xdata[0], ydata[0], xdata[2], ydata[2]

    time = (xdata[0], xdata[1])

    model = dc.AortaLiver2scan().load(state)
    model.set_params(dose2 = 0)

    params = tools.export_params(model, tb, Sb, tl, Sl, pars)
    params['T1_3']=['Liver T1-MOLLI at scan 2', pars['T1_liver_3'], 'sec', 0]
    # params['t3_MOLLI']=['Time of T1-MOLLI at scan 2', pars['T1_time_3']/(60*60), 'hrs', 0]

    params['t0']=["Start time first acquisition", (t0+time[0][0])/(60*60), 
                'hrs',0]
    # params['t1']=["End time first acquisition", (t0+time[0][-1])/(60*60), 
    #             'hrs',0]
    # params['t2']=["Start time second acquisition", (t0+time[1][0])/(60*60), 
    #             'hrs',0]
    params['t3']=["End time second acquisition", (t0+time[1][-1])/(60*60), 
                'hrs',0]
    # params['dt1']=["Time step first acquisition", model.params('TS'), 'sec',0]
    # params['dt2']=["Time step second acquisition", model.params('TS'), 'sec',0]

    params = tools.to_tristan_units(params)

    dmrpath = os.path.join(resultspath, 'Results')
    return tools.to_dmr(dmrpath, subj, visit, params)


