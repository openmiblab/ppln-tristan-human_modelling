import os
import logging


from miblab.pipe import run_stage

import utils.models_one_scan
import utils.models_two_scan
from utils import stats

MODELS = {
    1: utils.models_one_scan,
    2: utils.models_two_scan,
}


def run(build, logfile):
    dir_output = os.path.join(build, 'human_modelling', 'stage_4_stats')
    dir_data = os.path.join(build, 'human_modelling', 'stage_3_concatenate')

    for scans in [1, 2]:
        for drug in ['rifampicin', 'ciclosporin', 'metformin']:
            dir_results = os.path.join(dir_output, drug, f"scans_{scans}")
            file_data = os.path.join(dir_data, drug, f"scans_{scans}", 'all_results')

            stats.describe(file_data, dir_results)
            stats.ttest(file_data, dir_results)
            if drug in ['rifampicin', 'ciclosporin']:
                stats.ttest_translation(file_data, dir_results, drug)
            if scans==2:
                stats.ttest_diurnal(file_data, dir_results)


if __name__ == '__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, 'human_modelling', __file__)