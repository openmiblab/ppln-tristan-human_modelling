import os
import logging

import dcmri as dc

from miblab.pipe import run_stage
from tristan_human_modelling.utils import figs

PIPELINE = 'human_modelling'


def run(build, logfile):
    dir_output = os.path.join(build, 'human_modelling', 'stage_6_figs')
    dir_results = os.path.join(build, 'human_modelling', 'stage_3_concatenate')

    for scans in [[1], [1,2]]:
        dir_scans = f"scans_{len(scans)}"
        for study in ['healthy_rifampicin', 'healthy_ciclosporin', 'healthy_metformin', 'patients_rifampicin']:
            file_data = dc.fetch(f"tristan_humans_{study}")
            file_output = os.path.join(dir_results, dir_scans, 'all_results')
            dir_figs = os.path.join(dir_output, dir_scans, study)
            
            figs.effect_plot_combined(file_output, dir_figs, f"tristan_humans_{study}")
            figs.diurnal_k_combined(file_data, file_output, dir_figs, f"tristan_humans_{study}")


if __name__ == '__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, 'human_modelling', __file__)