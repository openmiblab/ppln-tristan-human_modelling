import os
import argparse

import pydmr


def bsa_dbdb(height, weight):
    # DuBois and DuBois formula

    # height in cm
    # weight in kg

    # returns BSA in m2

    BSA = 0.007184 * height**0.725 * weight**0.425

    return BSA


def egfr_ckd_epi(crea, age, sex, black):

    # crea in units of umol/L
    # age in yrs
    # sec: 1= Male, 0 = Female
    # black: True or False

    # returns eGFR in mL/min/1.73m2

    crea = 0.0113*crea # convert from umol/L to mg/dl
    if black:
        if sex==0:
            if crea<=0.7:
                return 166*(crea/0.7)**(-0.329) * 0.993**age
            else:
                return 166*(crea/0.7)**(-1.209) * 0.993**age
        else:
            if crea<=0.9:
                return 163*(crea/0.9)**(-0.411) * 0.993**age
            else:
                return 163*(crea/0.9)**(-1.209) * 0.993**age
    else:
        if sex==0:
            if crea<=0.7:
                return 144*(crea/0.7)**(-0.329) * 0.993**age
            else:
                return 144*(crea/0.7)**(-1.209) * 0.993**age
        else:
            if crea<=0.9:
                return 141*(crea/0.9)**(-0.411) * 0.993**age
            else:
                return 141*(crea/0.9)**(-1.209) * 0.993**age


def add_derived_variables(datafile, resultsfile):

    dmr = pydmr.read(datafile)
    # dmr['data'] = {k: v[:3] for k,v in dmr['data'].items()}

    # Extend data dictionary
    dmr['data']['eGFR'] = ['estimated Glomerular Filtration Rate (CKD-EPI)', 'mL/min/1.73m2', 'float']
    dmr['data']['BSA'] = ['Body surface area (duBois duBois)', 'm2', 'float']
    dmr['data']['eGFR_abs'] = ['estimated absolute Glomerular Filtration Rate (CKD-EPI)', 'mL/min', 'float']

    # Extend values
    pars = dmr['pars']
    subjects = {k[0] for k in pars}
    for subj in subjects:
        age = pars[(subj, 'screening', 'Age')]
        sex = pars[(subj, 'screening', 'Sex')]
        height = pars[(subj, 'screening', 'Height')]
        weight = pars[(subj, 'screening', 'Weight')]
        black = False
        bsa = bsa_dbdb(height, weight)
        dmr['pars'][(subj, 'screening', 'BSA')] = bsa

        for visit in ['control', 'drug', 'screening']:
            if (subj, visit, 'Crea') in pars:
                crea = pars[(subj, visit, 'Crea')]
                egfr = egfr_ckd_epi(crea, age, sex, black)
                dmr['pars'][(subj, visit, 'eGFR')] = egfr
                dmr['pars'][(subj, visit, 'eGFR_abs')] = egfr * (bsa / 1.73)


    pydmr.write(resultsfile, dmr)



if __name__ == '__main__':

    BASE = 'C:\\Users\\md1spsx\\Documents\\GitHub\\tristan-human-stage-2-modelling'
    DATA = os.path.join(BASE, 'data')
    BUILD = os.path.join(BASE, 'wip', 'build')

    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=str, default=DATA, help="Data folder")
    parser.add_argument("--build", type=str, default=BUILD, help="Build folder")
    args = parser.parse_args()

    add_derived_variables(args.data, args.build)