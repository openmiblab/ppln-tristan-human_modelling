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

    for scans in [1, 2]:
        for drug in ['rifampicin', 'ciclosporin', 'metformin']:
            dir_results = os.path.join(dir_output, drug, f"scans_{scans}")
            file = os.path.join(dir_data, f'tristan_humans_healthy_{drug}.dmr.zip')
            data = pydmr.read(file, format='nest')

            Parallel(n_jobs=-1)(delayed(train_subj)(subj, file, dir_results, scans) for subj in data['rois'].keys())
            
            # [(train_subj)(subj, file, dir_results, scans) for subj in data['rois'].keys()]


def train_subj(subj, file, dir_results, scans=2):
    state_path = os.path.join(dir_results, 'State')
    os.makedirs(state_path, exist_ok=True)
    state = os.path.join(state_path, subj)
    if os.path.exists(state):
        return

    MODELS[scans].train(state, file)

    # Optional: Plot on the fly
    MODELS[scans].save_plots(state, file, dir_results)



if __name__ == '__main__':
    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, 'human_modelling', __file__)