import os
import miblab
import pydmr

from methods import models


DRUG = 'metformin'
POPULATION = 'healthy'
PARTICIPANT = 'SHF-014'
VISIT = 'control'


if __name__ == '__main__':

    # Define paths
    repo_dir = os.getcwd()
    data_folder = os.path.join(repo_dir, 'data')
    build_folder = os.path.join(repo_dir, 'build')
    os.makedirs(build_folder, exist_ok=True)

    # Download the data
    dataset = f'tristan_humans_{POPULATION}_{DRUG}.dmr.zip'
    datafile = miblab.zenodo_fetch(dataset, data_folder)

    # Read data
    data = pydmr.read(datafile, format='nest')
    pars = data['pars'][PARTICIPANT][VISIT]
    rois = data['rois'][PARTICIPANT][VISIT]

    # Train a model to predict the data
    model = models.two_scan(data, PARTICIPANT, VISIT)
    
    # Output results
    file = os.path.join(build_folder, 'case_study_control.png')
    x, y = models.two_scan_data(data['rois'][PARTICIPANT][VISIT])
    model.plot(x, y, fname=file)
    model.print_params(round_to=5)

    # Print rate constants
    khe = 6000 * model.params('khe')
    kbh = 6000 * model.params('kbh')

    # Print effect sizes
    print(f'khe: {khe:.1f} mL/min/100cm3')
    print(f'kbh: {kbh:.2f} mL/min/100cm3')


