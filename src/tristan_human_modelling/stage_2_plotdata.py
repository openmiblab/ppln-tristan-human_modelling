import os
import logging
import sys

from joblib import Parallel, delayed

import dcmri as dc
from miblab.pipe import run_stage

import tristan_human_modelling.utils.models as models

MODELS = {
    1: models.plot_one_scan,
    2: models.plot_two_scan,
}

PIPELINE = 'human_modelling'


def run(build, logfile):
    dir_input = os.path.join(build, PIPELINE, 'stage_1_modelfit')
    dir_output = os.path.join(build, PIPELINE, 'stage_2_plotdata')
    for scans in [[1], [1,2]]:
    #for scans in [ [1]]:
        dir_input_scans = os.path.join(dir_input, f"scans_{len(scans)}")
        dir_output_scans = os.path.join(dir_output, f"scans_{len(scans)}")
        data = dc.load_tristan_humans(scans=scans)
        success = Parallel(n_jobs=-1)(delayed(plot_subject)(subj, data, dir_input_scans, dir_output_scans, logfile, scans) for subj in data.keys())
        #success = [plot_subject(subj, data, dir_input_scans, dir_output_scans, logfile, scans) for subj in data.keys()]

        logging.info(f"--- {len(scans)} scans --- ")
        success = [s for s in success if s is not None]
        if success == []:
            logging.info(f"Found all results on disk. Nothing else to do here.")
        else:
            logging.info(f"Successfully plotted {sum(success)} / {len(success)} cases.")


def plot_subject(subj, data, dir_input, dir_output, logfile, scans):
    # Re-initialize logging inside the child process
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(logfile),
            logging.StreamHandler(sys.stdout),
        ],
        force=True # Required to ensure the child process adopts this
    )

    return MODELS[len(scans)](subj, data[subj], dir_input, dir_output) 



if __name__ == '__main__':
    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, PIPELINE, __file__)