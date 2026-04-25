import os
import time

import miblab
import pydmr
from joblib import Parallel, delayed

# For two-scan protocol:
# import methods.models as models
# One-scan protocol:
import methods.models_one_scan
import methods.models_two_scan
import methods.preprocess
import methods.figs
import methods.tables

MODELS = {
    1: methods.models_one_scan,
    2: methods.models_two_scan,
}
YLIM = [40,40,4,4]
ylim = {
    'rifampicin': YLIM,
    'ciclosporin': YLIM,
    'metformin': YLIM,
}

def all(datapath, resultspath, drug, subject=None, scans=2, patients=False):

    # Define paths
    buildpath = os.path.join(resultspath, f'build_scans_{scans}')
    if patients:
        dataset = f'tristan_humans_patients_{drug}.dmr.zip'
        resultspath = os.path.join(buildpath, f'{drug}_clinical')
        data = os.path.join(resultspath, f'tristan_humans_patients_{drug}_deriv')
    else:
        dataset = f'tristan_humans_healthy_{drug}.dmr.zip'
        resultspath = os.path.join(buildpath, drug)
        data = os.path.join(resultspath, f'tristan_humans_healthy_{drug}_deriv')
    os.makedirs(resultspath, exist_ok=True)

    source = os.path.join(datapath, dataset)
    methods.preprocess.add_derived_variables(source, data)

    # Note
    # This is where the clinical data were saved for the first in human study:
    # datafile = miblab.zenodo_fetch(dataset, datapath, "15610541")
    # Obsolete
    # dataset = f'tristan_humans_healthy_{drug}.dmr.zip'
    # datafile = os.path.join(datapath, dataset)
    # When finished fetch from here
    # datafile = miblab.zenodo_fetch(dataset, datapath)
    
    # Compute
    train(data, resultspath, subject=subject, scans=scans)
    plot(data, resultspath, subject=subject, scans=scans)
    to_dmr(data, resultspath, subject=subject, scans=scans)
    methods.figs.effect_plot_combined(resultspath, ylim=[40,4])
    methods.figs.diurnal_k_combined(data, resultspath, ylim=ylim[drug])
    if patients:
        drug = f"{drug}_clinical" # bad hack
    methods.tables.describe(buildpath, drug)
    methods.tables.ttest(buildpath, drug)
    methods.tables.ttest_diurnal(buildpath, drug, scans)
    methods.tables.ttest_translation(buildpath, drug)
    methods.tables.outcomes(buildpath, drug)
    methods.tables.outcomes_diurnal(buildpath, drug, scans)
    methods.tables.constants(buildpath, drug)



def train(datafile, resultspath, subject=None, scans=2):
    start = time.time()
    data = pydmr.read(datafile, format='nest')
    # [train_subj(subj, datafile, resultspath, subject, scans) for subj in data['rois'].keys()]
    Parallel(n_jobs=-1)(delayed(train_subj)(subj, datafile, resultspath, subject, scans) for subj in data['rois'].keys())
    print('Calculation time (mins): ', (time.time()-start)/60)


def plot(datafile, resultspath, subject=None, scans=2):
    data = pydmr.read(datafile, format='nest')
    [plot_subj(subj, datafile, resultspath, subject, scans) for subj in data['rois'].keys()]
    #tasks = [dask.delayed(plot_subj)(subj, datafile, resultspath, subject) for subj in data['rois'].keys()]
    #dask.compute(*tasks)


def train_subj(subj, datafile, resultspath, subject=None, scans=2):
    if subject is not None:
        if subj != subject:
            return 
    state_path = os.path.join(resultspath, 'State')
    os.makedirs(state_path, exist_ok=True)
    state = os.path.join(state_path, subj)
    if os.path.exists(state):
        return
    # state_init_path = os.path.join(resultspath, 'StateInit')
    # state_init = os.path.join(state_init_path, subj")
    state_init = ''
    MODELS[scans].train(state, datafile, state_init)

    # Plot on the fly
    MODELS[scans].save_plots(state, datafile, resultspath)


def plot_subj(subj, datafile, resultspath, subject=None, scans=2):
    if subject is not None:
        if subj != subject:
            return
    state_path = os.path.join(resultspath, 'State')
    os.makedirs(state_path, exist_ok=True)
    state = os.path.join(state_path, subj)
    if not os.path.exists(state):
        return
    MODELS[scans].save_plots(state, datafile, resultspath)



def to_dmr(datafile, resultspath, subject=None, scans=2):

    state_path = os.path.join(resultspath, 'State')
    if not os.path.exists(state_path):
        return
    data = pydmr.read(datafile, format='nest')

    results = []
    for subj in data['rois'].keys():
        if subject is not None:
            if subj != subject:
                continue
        state = os.path.join(state_path, subj)
        if not os.path.exists(state):
            continue
        file = MODELS[scans].save_results(state, datafile, resultspath)
        results.append(file)

    file = os.path.join(resultspath, 'Results', 'all_results')
    pydmr.concat(results, file)
    wide_file = os.path.join(resultspath, 'Results', 'all_results_wide.csv')
    pydmr.pars_to_wide(file, wide_file)




