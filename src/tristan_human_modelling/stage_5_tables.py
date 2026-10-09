import os
import logging

from miblab.pipe import run_stage

import utils.models_one_scan
import utils.models_two_scan
from utils import tables

MODELS = {
    1: utils.models_one_scan,
    2: utils.models_two_scan,
}


def run(build, logfile):
    dir_stats = os.path.join(build, 'human_modelling', 'stage_4_stats')
    dir_output = os.path.join(build, 'human_modelling', 'stage_5_tables')

    for scans in [1, 2]:
        for study in ['rifampicin', 'ciclosporin', 'metformin', 'patients_rifampicin']:
            dir_study_stats = os.path.join(dir_stats, study, f"scans_{scans}")
            dir_tables = os.path.join(dir_output, study, f"scans_{scans}")

            tables.outcomes(dir_study_stats, dir_tables)
            tables.outcomes_diurnal(dir_study_stats, dir_tables, scans)
            tables.constants(dir_study_stats, dir_tables)



if __name__ == '__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, 'human_modelling', __file__)