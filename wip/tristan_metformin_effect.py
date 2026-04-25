import os
import argparse

import compute_results


def main(datapath, buildpath, scans=2):

    subject = None
    #subject = 'SHF-008'
    drug = 'metformin'

    compute_results.all(datapath, buildpath, drug, subject=subject, scans=scans)
    
    # report.all_results(
    #     results, 
    #     'report (complete)',
    #     title = 'Leeds pilot study',
    #     subtitle = f'{drug} (all results)',
    #     subject = 'D2.13 - Internal report',
    # )
    # report.key_results(
    #     results, 
    #     'report (summary)',
    #     title = 'Leeds pilot study',
    #     subtitle = f'{drug} (key results)',
    #     subject = 'D2.13 - Internal report',
    # )


if __name__ == '__main__':

    BASE = 'C:\\Users\\md1spsx\\Documents\\GitHub\\tristan-human-stage-2-modelling'
    DATA = os.path.join(BASE, 'data')
    BUILD = os.path.join(BASE, 'wip', 'build')

    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=str, default=DATA, help="Data folder")
    parser.add_argument("--build", type=str, default=BUILD, help="Build folder")
    parser.add_argument("--scans", type=int, default=1, help="Nr of scans")
    args = parser.parse_args()

    main(args.data, args.build, args.scans)

