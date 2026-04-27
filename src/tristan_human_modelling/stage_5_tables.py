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
    dir_output = os.path.join(build, 'human_modelling', 'stage_5_tables')
    for scans in [1, 2]:
        for drug in ['rifampicin', 'ciclosporin', 'metformin']:
            dir_stats = os.path.join(build, 'human_modelling', 'stage_4_stats', drug, f"scans_{scans}")
            dir_tables = os.path.join(dir_output, drug, f"scans_{scans}")

            tables.outcomes(dir_stats, dir_tables)
            tables.outcomes_diurnal(dir_stats, dir_tables, scans)
            tables.constants(dir_stats, dir_tables)



if __name__ == '__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, 'human_modelling', __file__)