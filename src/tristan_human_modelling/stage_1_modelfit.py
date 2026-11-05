import os
import logging
import sys

from joblib import Parallel, delayed

import dcmri as dc
from miblab.pipe import run_stage

import tristan_human_modelling.utils.models as models

PIPELINE = 'human_modelling'


#  TODO: try la and gut plucom for better first pass

def run(build, logfile):
    # run_one_scan(build, logfile)
    run_two_scans(build, logfile)


def run_one_scan(build, logfile):
    dir_output = os.path.join(build, PIPELINE, 'stage_1_modelfit')
    subjects = dc.load_tristan_human_subjects()

    scans = [1]

    dir_results = os.path.join(dir_output, f"scans_{len(scans)}")
    data = dc.load_tristan_humans(scans=scans)
    success = Parallel(n_jobs=-1)(delayed(train_subject_one_scan)(subj, data, subjects, dir_results, logfile) for subj in subjects.keys())
    # success = [train_subject_one_scan(subj, data, dir_results, logfile) for subj in data.keys()]

    logging.info(f"--- {len(scans)} scans --- ")
    success = [s for s in success if s is not None]
    if success == []:
        logging.info(f"Found all results on disk. Nothing else to do here.")
    else:
        logging.info(f"Successfully trained {sum(success)} / {len(success)} cases.")


def run_two_scans(build, logfile):
    dir_output = os.path.join(build, PIPELINE, 'stage_1_modelfit')
    subjects = dc.load_tristan_human_subjects()

    scans = [1, 2]

    dir_results = os.path.join(dir_output, f"scans_{len(scans)}")
    dir_results_one_scan = os.path.join(dir_output, f"scans_1")
    data = dc.load_tristan_humans(scans=scans)
    success = Parallel(n_jobs=-1)(delayed(train_subject_two_scans)(subj, data, subjects, dir_results, logfile, dir_results_one_scan) for subj in subjects.keys())
    # success = [train_subject_two_scans(subj, data, dir_results, logfile) for subj in data.keys()]

    logging.info(f"--- {len(scans)} scans --- ")
    success = [s for s in success if s is not None]
    if success == []:
        logging.info(f"Found all results on disk. Nothing else to do here.")
    else:
        logging.info(f"Successfully trained {sum(success)} / {len(success)} cases.")



def train_subject_one_scan(subj, data, subjects, dir_results, logfile):
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
    return models.one_scan(subj, data[subj], subjects, dir_results) 


def train_subject_two_scans(subj, data, subjects, dir_results, logfile, dir_results_one_scan=None):
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
    return models.two_scan(subj, data[subj], subjects, dir_results, dir_results_one_scan) 



if __name__ == '__main__':
    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, PIPELINE, __file__)