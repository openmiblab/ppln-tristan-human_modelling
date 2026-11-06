import os
import logging


from miblab.pipe import run_stage

from tristan_human_modelling.utils import stats

PIPELINE = 'human_modelling'


def run(build, logfile):
    dir_input = os.path.join(build, PIPELINE, 'stage_3_concatenate')
    dir_output = os.path.join(build, PIPELINE, 'stage_4_stats')
    
    for scans in [[1], [1,2]]:
        dir_scans = f"scans_{len(scans)}"
        file_data = os.path.join(dir_input, dir_scans, 'all_results')

        # Per study
        studies = ['healthy_rifampicin', 'healthy_ciclosporin', 'healthy_metformin', 'patients_rifampicin']
        for study in studies:
            dir_results = os.path.join(dir_output, dir_scans, study)
        
            stats.describe(file_data, dir_results, [f"tristan_humans_{study}"])
            stats.ttest(file_data, dir_results, study)
            if len(scans)==2:
                stats.ttest_diurnal(file_data, dir_results, study)

        # All healthy volunteers
        studies = ['healthy_rifampicin', 'healthy_ciclosporin', 'healthy_metformin']
        dir_results = os.path.join(dir_output, dir_scans, 'healthy')
        stats.describe(file_data, dir_results, [f"tristan_humans_{study}" for study in studies])


if __name__ == '__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, 'human_modelling', __file__)