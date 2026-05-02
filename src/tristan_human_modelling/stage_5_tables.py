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

    run_volunteers(dir_stats, dir_output)
    run_patients(dir_stats, dir_output)


def run_volunteers(dir_input, dir_output):
    for scans in [1, 2]:
        for drug in ['rifampicin', 'ciclosporin', 'metformin']:
            dir_stats = os.path.join(dir_input, drug, f"scans_{scans}")
            dir_tables = os.path.join(dir_output, drug, f"scans_{scans}")

            tables.outcomes(dir_stats, dir_tables)
            tables.outcomes_diurnal(dir_stats, dir_tables, scans)
            tables.constants(dir_stats, dir_tables)


def run_patients(dir_input, dir_output):
    for scans in [1, 2]:
        dir_stats = os.path.join(dir_input, 'patients_rifampicin', f"scans_{scans}")
        dir_tables = os.path.join(dir_output, 'patients_rifampicin', f"scans_{scans}")

        tables.outcomes(dir_stats, dir_tables)
        tables.outcomes_diurnal(dir_stats, dir_tables, scans)
        tables.constants(dir_stats, dir_tables)



if __name__ == '__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, 'human_modelling', __file__)