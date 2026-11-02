import os
import logging

from miblab.pipe import run_stage
from utils import figs

import utils.models

MODELS = {
    1: utils.models.one_scan,
    2: utils.models.two_scan,
}

PIPELINE = 'human_modelling'


def run(build, logfile):
    dir_output = os.path.join(build, 'human_modelling', 'stage_6_figs')
    dir_data = os.path.join(build, 'human_modelling', 'stage_1_download')
    dir_results = os.path.join(build, 'human_modelling', 'stage_3_concatenate')

    for scans in [1, 2]:
        for study in ['rifampicin', 'ciclosporin', 'metformin', 'patients_rifampicin']:
            file_data = os.path.join(dir_data, f'dataset.dmr.zip')
            file_output = os.path.join(dir_results, f"scans_{scans}", 'all_results')
            dir_figs = os.path.join(dir_output, study, f"scans_{scans}")
            
            figs.effect_plot_combined(file_output, dir_figs, study)
            figs.diurnal_k_combined(file_data, file_output, dir_figs, study)





if __name__ == '__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, 'human_modelling', __file__)