import os
import logging

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
    fit_volunteers(dir_data, dir_output)
    fit_patients(dir_data, dir_output)
    # plot_volunteers(dir_data, dir_output)
    # plot_patients(dir_data, dir_output)


def fit_volunteers(dir_data, dir_output):
    for scans in [2, 1]:
        for drug in ['rifampicin', 'ciclosporin', 'metformin']:
            dir_results = os.path.join(dir_output, drug, f"scans_{scans}")
            file = os.path.join(dir_data, f'tristan_humans_healthy_{drug}.dmr.zip')
            data = pydmr.read(file, format='nest')
            Parallel(n_jobs=-1)(delayed(train_subject)(subj, file, dir_results, scans) for subj in data['rois'].keys())
            #[(train_subject)(subj, file, dir_results, scans) for subj in data['rois'].keys()]

def fit_patients(dir_data, dir_output):
    for scans in [2, 1]:
        dir_results = os.path.join(dir_output, 'patients_rifampicin', f"scans_{scans}")
        file = os.path.join(dir_data, f'tristan_humans_patients_rifampicin.dmr.zip')
        data = pydmr.read(file, format='nest')
        #staged=2 if scans==2 else 1
        Parallel(n_jobs=-1)(delayed(train_subject)(subj, file, dir_results, scans) for subj in data['rois'].keys())
        #[(train_subject)(subj, file, dir_results, scans, staged=staged) for subj in data['rois'].keys()]

def plot_volunteers(dir_data, dir_output):
    for scans in [1, 2]:
        for drug in ['rifampicin', 'ciclosporin', 'metformin']:
            dir_results = os.path.join(dir_output, drug, f"scans_{scans}")
            file = os.path.join(dir_data, f'tristan_humans_healthy_{drug}.dmr.zip')
            data = pydmr.read(file, format='nest')
            Parallel(n_jobs=-1)(delayed(plot_subject)(subj, file, dir_results, scans) for subj in data['rois'].keys())

def plot_patients(dir_data, dir_output):
    for scans in [1, 2]:
        dir_results = os.path.join(dir_output, 'patients_rifampicin', f"scans_{scans}")
        file = os.path.join(dir_data, f'tristan_humans_patients_rifampicin.dmr.zip')
        data = pydmr.read(file, format='nest')
        Parallel(n_jobs=-1)(delayed(plot_subject)(subj, file, dir_results, scans) for subj in data['rois'].keys())


def train_subject(subj, file, dir_results, scans=2, staged=0):
    state_path = os.path.join(dir_results, 'State')
    os.makedirs(state_path, exist_ok=True)
    state = os.path.join(state_path, subj)
    MODELS[scans].train(state, file, staged=staged)
    MODELS[scans].save_plots(state, file, dir_results) # plot on the fly - optional

def plot_subject(subj, file, dir_results, scans=2):
    state_path = os.path.join(dir_results, 'State')
    state = os.path.join(state_path, subj)
    if os.path.exists(state):
        MODELS[scans].save_plots(state, file, dir_results)



if __name__ == '__main__':
    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, 'human_modelling', __file__)