import os
import logging

from miblab.pipe import run_stage

import utils.models_one_scan
import utils.models_two_scan
from utils import figs

MODELS = {
    1: utils.models_one_scan,
    2: utils.models_two_scan,
}


def run(build, logfile):
    dir_output = os.path.join(build, 'human_modelling', 'stage_6_figs')
    dir_data = os.path.join(build, 'human_modelling', 'stage_1_download')
    dir_results = os.path.join(build, 'human_modelling', 'stage_3_concatenate')

    for scans in [1, 2]:
        for drug in ['rifampicin', 'ciclosporin', 'metformin']:
            file_data = os.path.join(dir_data, f'tristan_humans_healthy_{drug}.dmr.zip')
            file_output = os.path.join(dir_results, drug, f"scans_{scans}", 'all_results')
            dir_figs = os.path.join(dir_output, drug, f"scans_{scans}")
            
            figs.effect_plot_combined(file_output, dir_figs, drug)
            figs.diurnal_k_combined(file_data, file_output, dir_figs, drug)



if __name__ == '__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, 'human_modelling', __file__)