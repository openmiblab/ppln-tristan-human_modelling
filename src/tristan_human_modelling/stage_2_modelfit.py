import os
import logging
import sys

from joblib import Parallel, delayed

import pydmr
from miblab.pipe import run_stage

import utils.models_one_scan
import utils.models_two_scan

MODELS = {
    1: utils.models_one_scan,
    2: utils.models_two_scan,
}


def run(build, logfile):
    dir_output = os.path.join(build, 'human_modelling', 'stage_2_modelfit')
    dir_data = os.path.join(build, 'human_modelling', 'stage_1_download')
    fit_data(dir_data, dir_output, logfile)
    # plot_data(dir_data, dir_output)


def fit_data(dir_data, dir_output, logfile):
    for scans in [1, 2]:
        dir_results = os.path.join(dir_output, f"scans_{scans}")
        file = os.path.join(dir_data, f'dataset.dmr.zip')
        data = pydmr.read(file, format='nest')
        Parallel(n_jobs=-1)(delayed(train_subject)(subj, file, dir_results, logfile, scans, staged=0) for subj in data['rois'].keys())
        #[(train_subject)(subj, file, dir_results, logfile, scans, staged=0) for subj in data['rois'].keys()]

def plot_data(dir_data, dir_output):
    for scans in [2, 1]:
        dir_results = os.path.join(dir_output, f"scans_{scans}")
        file = os.path.join(dir_data, f'dataset.dmr.zip')
        data = pydmr.read(file, format='nest')
        Parallel(n_jobs=-1)(delayed(plot_subject)(subj, file, dir_results, scans) for subj in data['rois'].keys())



def train_subject(subj, file, dir_results, logfile, scans=2, staged=0):
    # Re-initialize logging inside the child process
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(logfile),
            logging.StreamHandler(sys.stdout),
        ],
        force=True # Required to ensure the child process adopts this
    )
    logging.info(f"Training subject {subj} with {scans} scans")

    # Compute
    state_path = os.path.join(dir_results, 'State')
    os.makedirs(state_path, exist_ok=True)
    state = os.path.join(state_path, subj)
    pars_path = os.path.join(dir_results, 'Pars')
    os.makedirs(pars_path, exist_ok=True)

    MODELS[scans].train(state, file, pars_path, state_init=None, staged=staged) 
    MODELS[scans].save_plots(state, file, dir_results) # plot on the fly - optional


def plot_subject(subj, file, dir_results, scans=2):
    state_path = os.path.join(dir_results, 'State')
    state = os.path.join(state_path, subj)
    if os.path.exists(state):
        MODELS[scans].save_plots(state, file, dir_results)



if __name__ == '__main__':
    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, 'human_modelling', __file__)