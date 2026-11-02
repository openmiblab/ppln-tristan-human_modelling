import os
import logging


from miblab.pipe import run_stage

from utils import stats

PIPELINE = 'human_modelling'


def run(build, logfile):
    dir_input = os.path.join(build, PIPELINE, 'stage_2_concatenate')
    dir_output = os.path.join(build, PIPELINE, 'stage_3_stats')
    
    for scans in [[1], [1,2]]:
        dir_scans = f"scans_{len(scans)}"
        for study in ['rifampicin', 'ciclosporin', 'metformin', 'patients_rifampicin']:
            dir_results = os.path.join(dir_output, study, dir_scans)
            file_data = os.path.join(dir_input, dir_scans, 'all_results')

            stats.describe(file_data, dir_results, study)
            stats.ttest(file_data, dir_results, study)
            if len(scans)==2:
                stats.ttest_diurnal(file_data, dir_results, study)

    # run_volunteers(dir_output, dir_data)
    # run_patients(dir_output, dir_data)


def run_volunteers(dir_output, dir_data):

    for scans in [[1], [1,2]]:
        dir_scans = f"scans_{len(scans)}"
        for drug in ['rifampicin', 'ciclosporin', 'metformin']:
            dir_results = os.path.join(dir_output, drug, dir_scans)
            file_data = os.path.join(dir_data, drug, dir_scans, 'all_results')

            stats.describe(file_data, dir_results)
            stats.ttest(file_data, dir_results)
            if scans==2:
                stats.ttest_diurnal(file_data, dir_results)

def run_patients(dir_output, dir_data):

    for scans in [[1], [1,2]]:
        dir_scans = f"scans_{len(scans)}"
        dir_results = os.path.join(dir_output, 'patients_rifampicin', dir_scans)
        file_data = os.path.join(dir_data, 'patients_rifampicin', dir_scans, 'all_results')

        stats.describe(file_data, dir_results)
        stats.ttest(file_data, dir_results)
        if scans==2:
            stats.ttest_diurnal(file_data, dir_results)


if __name__ == '__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, 'human_modelling', __file__)