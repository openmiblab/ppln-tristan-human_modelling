"""
Automatic download of data from Zenodo.
"""
import os

from miblab import pipe
from miblab_data.zenodo import download

from utils.preprocess import add_derived_variables

# Download from https://zenodo.org/records/19760095

def run(build, logfile):
    dir_output = os.path.join(build, 'human_modelling', 'stage_1_download')

    for dataset in [
        'tristan_humans_healthy_rifampicin.dmr.zip',
        'tristan_humans_healthy_metformin.dmr.zip',
        'tristan_humans_healthy_ciclosporin.dmr.zip',
        'tristan_humans_healthy_controls.dmr.zip',
        'tristan_humans_patients_rifampicin.dmr.zip',
    ]:
        file = download('19760095', dataset, dir_output)
        add_derived_variables(file, file)


if __name__ == '__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    pipe.run_stage(run, build, 'human_modelling', __file__)