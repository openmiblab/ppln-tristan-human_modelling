import os
import logging

from miblab.pipe import run_stage

from tristan_human_modelling.utils import tables

PIPELINE = 'human_modelling'


def run(build, logfile):
    dir_stats = os.path.join(build, PIPELINE, 'stage_4_stats')
    dir_output = os.path.join(build, PIPELINE, 'stage_5_tables')

    for scans in [[1], [1,2]]:
        dir_scans = f"scans_{len(scans)}"

        # Per study
        for study in ['healthy_rifampicin', 'healthy_ciclosporin', 'healthy_metformin', 'patients_rifampicin']:
            dir_study_stats = os.path.join(dir_stats, dir_scans, study)
            dir_tables = os.path.join(dir_output, dir_scans, study)

            tables.outcomes(dir_study_stats, dir_tables)
            tables.constants(dir_study_stats, dir_tables)
            if len(scans)==2:
                tables.outcomes_diurnal(dir_study_stats, dir_tables, scans)

        # Constants for all healthy volunteers
        dir_study_stats = os.path.join(dir_stats, dir_scans, 'healthy')
        dir_tables = os.path.join(dir_output, dir_scans, 'healthy')
        tables.constants(dir_study_stats, dir_tables)


if __name__ == '__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, PIPELINE, __file__)