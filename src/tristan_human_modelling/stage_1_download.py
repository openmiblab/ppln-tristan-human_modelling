"""
Automatic download of data from Zenodo.
"""
import os
import pydmr

from miblab import pipe
from miblab_data.zenodo import download

from utils.preprocess import add_derived_variables

STUDY = {
    'tristan_humans_healthy_rifampicin.dmr.zip': 'rifampicin',
    'tristan_humans_healthy_metformin.dmr.zip': 'metformin',
    'tristan_humans_healthy_ciclosporin.dmr.zip': 'ciclosporin',
    'tristan_humans_patients_rifampicin.dmr.zip': 'patients_rifampicin',
}

# Download from https://zenodo.org/records/20210596

def run(build, logfile):
    dir_output = os.path.join(build, 'human_modelling', 'stage_1_download')
    
    files = []
    for dataset, study in STUDY.items():
        file = download('20210596', dataset, dir_output)
        add_derived_variables(file, file, study)
        files.append(file)

    result = os.path.join(dir_output, 'dataset.dmr')
    pydmr.concat(files, result, cleanup=True)
    wide_file = os.path.join(dir_output, 'dataset_wide.csv')
    pydmr.pars_to_wide(result, wide_file)


if __name__ == '__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    pipe.run_stage(run, build, 'human_modelling', __file__)