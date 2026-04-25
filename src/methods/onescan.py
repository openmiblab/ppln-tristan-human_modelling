import os
import time

import dcmri as dc
import pydmr
import numpy as np

from methods import tools, models


def compute(datafile, resultspath):

    start = time.time()

    if not os.path.exists(resultspath):
        os.makedirs(resultspath)

    results = []
    data = pydmr.read(datafile, format='nest')
    for subj in data['rois'].keys():
        for visit in data['rois'][subj].keys():

            # Train model
            model = models.one_scan(data, subj, visit, verbose=2)

            # Save results
            save_plots(model, data, subj, visit, resultspath)
            file = save_results(model, data, subj, visit, resultspath)

            results.append(file)

    file = os.path.join(resultspath, 'all_results')
    pydmr.concat(results, file)

    print('Calculation time (mins): ', (time.time()-start)/60)


def compute_vart(datafile, resultspath, 
                 acq_times = [5,10,15,20,25,30,35,40]):

    start = time.time()

    if not os.path.exists(resultspath):
        os.makedirs(resultspath)

    results = []
    data = pydmr.read(datafile, format='nest')
    for subj in data['rois'].keys():
        for visit in data['rois'][subj].keys():
            for tacq in acq_times:
                
                # Train model
                model = models.one_scan(data, subj, visit, verbose=2, tacq=tacq)

                # Save results
                save_plots(model, data, subj, visit, resultspath, tacq=tacq)
                file = save_results(model, data, subj, visit, resultspath, tacq=tacq)

                results.append(file)
    file = os.path.join(resultspath, 'all_results')
    pydmr.concat(results, file)
    
    print('Calculation time (mins): ', (time.time()-start)/60)





def save_plots(model, data, subj, visit, path, tacq=None):

    rois = data['rois'][subj][visit]
    pars = data['pars'][subj][visit]

    study = visit if tacq is None else visit + '_' + str(tacq).zfill(2)
    name = subj + '_' + study
    xdata, ydata = models.one_scan_data(rois)
    t0 = rois['time_1'][0]
    path = os.path.join(path, 'Plots')
    file = os.path.join(path, name)
    
    t = [
        0, 
        #pars['T1_time_2']-t0,
    ]
    R1a = [
        1/pars['T1_aorta_1'], 
        #1/pars['T1_aorta_2'],
    ]
    R1l = [
        1/pars['T1_liver_1'], 
        #1/pars['T1_liver_2'],
    ]

    if not os.path.exists(path):
        os.makedirs(path)
    
    ya = [dc.signal_ss(model.params('S0a'), R1a[0], model.params('TR'), model.params('FA')),
          #dc.signal_ss(model.pars['S0a'], R1a[1], model.pars['TR'], model.pars['FA']),
          ]
    yl = [dc.signal_ss(model.params('S0l'), R1l[0], model.params('TR'), model.params('FA')),
          #dc.signal_ss(model.pars['S0l'], R1l[1], model.pars['TR'], model.pars['FA']),
          ]
    test=((t,ya),(t,yl))

    BAT = model.params('BAT')
    model.plot(xdata, ydata, 
               fname=file + '.png', ref=test, show=False)
    model.plot(xdata, ydata, xlim=[xdata[0][0], xdata[0][-1]], 
               fname=file + '_win1.png', ref=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+600], 
               fname=file + '_win2.png', ref=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+160], 
               fname=file + '_win3.png', ref=test, show=False) 


def save_results(model, data, subj, visit, path, tacq=None):

    study = visit if tacq is None else visit + '_' + str(tacq).zfill(2)

    # Save state
    state_path = os.path.join(path, 'State')
    os.makedirs(state_path, exist_ok=True)
    model.save(os.path.join(state_path, f"{subj}_{study}"))

    rois = data['rois'][subj][visit]
    pars = data['pars'][subj][visit]

    xdata, ydata = models.one_scan_data(rois)
    tb, Sb, tl, Sl = xdata[0], ydata[0], xdata[1], ydata[1]

    params = tools.export_params(model, tb, Sb, tl, Sl, pars)
    params = tools.to_tristan_units(params)
    
    dmrpath = os.path.join(path, 'Results')
    return tools.to_dmr(dmrpath, subj, study, params)