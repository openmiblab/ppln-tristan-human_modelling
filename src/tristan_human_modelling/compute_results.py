import os
import time

import miblab
import pydmr
from joblib import Parallel, delayed

# For two-scan protocol:
# import utils.models as models
# One-scan protocol:
import utils.models_one_scan
import utils.models_two_scan
import utils.preprocess
import utils.figs
import utils.tables

MODELS = {
    1: utils.models_one_scan,
    2: utils.models_two_scan,
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
    utils.preprocess.add_derived_variables(source, data)

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
    utils.figs.effect_plot_combined(resultspath, ylim=[40,4])
    utils.figs.diurnal_k_combined(data, resultspath, ylim=[40,40,4,4])
    if patients:
        drug = f"{drug}_clinical" # bad hack
    utils.tables.describe(buildpath, drug)
    utils.tables.ttest(buildpath, drug)
    utils.tables.ttest_diurnal(buildpath, drug, scans)
    utils.tables.ttest_translation(buildpath, drug)
    utils.tables.outcomes(buildpath, drug)
    utils.tables.outcomes_diurnal(buildpath, drug, scans)
    utils.tables.constants(buildpath, drug)



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




