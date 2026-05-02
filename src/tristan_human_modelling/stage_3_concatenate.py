import os
import logging


import pydmr
from miblab.pipe import run_stage

import utils.models_one_scan
import utils.models_two_scan

MODELS = {
    1: utils.models_one_scan,
    2: utils.models_two_scan,
}


def run(build, logfile):
    dir_output = os.path.join(build, 'human_modelling', 'stage_3_concatenate')
    dir_data = os.path.join(build, 'human_modelling', 'stage_1_download')
    dir_fits = os.path.join(build, 'human_modelling', 'stage_2_modelfit')

    run_volunteers(dir_fits, dir_output, dir_data)
    run_patients(dir_fits, dir_output, dir_data)


def run_volunteers(dir_fits, dir_output, dir_data):
    for scans in [1, 2]:
        for drug in ['rifampicin', 'ciclosporin', 'metformin']:
            dir_state = os.path.join(dir_fits, drug, f"scans_{scans}", 'State')
            dir_results = os.path.join(dir_output, drug, f"scans_{scans}")
            datafile = os.path.join(dir_data, f'tristan_humans_healthy_{drug}.dmr.zip')
            to_dmr(datafile, dir_state, dir_results, scans, drug)


def run_patients(dir_fits, dir_output, dir_data):
    for scans in [1, 2]:
        dir_state = os.path.join(dir_fits, 'patients_rifampicin', f"scans_{scans}", 'State')
        dir_results = os.path.join(dir_output, 'patients_rifampicin', f"scans_{scans}")
        datafile = os.path.join(dir_data, f'tristan_humans_patients_rifampicin.dmr.zip')
        to_dmr(datafile, dir_state, dir_results, scans, 'rifampicin')


def to_dmr(datafile, state_path, dir_results, scans, drug):

    if not os.path.exists(state_path):
        return
    data = pydmr.read(datafile, format='nest')

    results = []
    for subj in data['rois'].keys():
        state = os.path.join(state_path, subj)
        if not os.path.exists(state):
            continue
        file = MODELS[scans].save_results(state, datafile, dir_results, drug)
        results.append(file)

    file = os.path.join(dir_results, 'all_results')
    pydmr.concat(results, file, cleanup=True)
    wide_file = os.path.join(dir_results, 'all_results_wide.csv')
    pydmr.pars_to_wide(file, wide_file)


if __name__ == '__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, 'human_modelling', __file__)