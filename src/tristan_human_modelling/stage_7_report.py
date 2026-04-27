import os
import logging

from miblab.pipe import run_stage

from utils import report


def run(build, logfile):
    dir_pipeline = os.path.join(build, 'human_modelling')
    dir_output = os.path.join(dir_pipeline, 'stage_7_report')
    report.two_scans(dir_pipeline, dir_output)
    report.one_scan(dir_pipeline, dir_output)


if __name__ == '__main__':
    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    run_stage(run, build, 'human_modelling', __file__)