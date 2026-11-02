import os
import logging

from miblab.pipe import run_stage

from utils import tables

PIPELINE = 'human_modelling'


def run(build, logfile):
    dir_stats = os.path.join(build, PIPELINE, 'stage_4_stats')
    dir_output = os.path.join(build, PIPELINE, 'stage_5_tables')

    for scans in [[1], [1,2]]:
        dir_scans = f"scans_{len(scans)}"
        for study in ['rifampicin', 'ciclosporin', 'metformin', 'patients_rifampicin']:
            dir_study_stats = os.path.join(dir_stats, study, dir_scans)
            dir_tables = os.path.join(dir_output, study, dir_scans)

            tables.outcomes(dir_study_stats, dir_tables)
            tables.outcomes_diurnal(dir_study_stats, dir_tables, scans)
            tables.constants(dir_study_stats, dir_tables)



if __name__ == '__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, PIPELINE, __file__)