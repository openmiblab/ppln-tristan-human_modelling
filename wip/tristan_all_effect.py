"""
TRISTAN analysis human data - changes compared to first in human study

Bias correction:
- Relaxivity correction added (khe).
- Human khe relative to plasma rather than blood.

Parametrization:
- Joint fit of control and treatment data.
- Fixed values for parameters that are not affected by the drug. 
- No diurnal variation on kbh in the control visit
- Coupled aorta and liver concentrations by replacing Eb parameter 
  using GFR (constant) and CL. 
-> 25 free parameters per subject (8 curves)

Changes - implementation detail:
- Fitting normalized dimensionless parameters in [0,1]
- Internal time resolution to 0.1 sec instead of default 0.5. This is because with 
  a quarter dose injection the injection duration can be as short as 0.7 sec, needing 
  a higher temporal resolution to capture the input function shape accurately. 
  This comes at a cost of substantially increased computation time.
- Injected dose normalized to ensure exact dose is injected in the discrete time axis

Changes - data:
- Precontrast aorta T1 values in GE cases replaced by literature 
  values due to inflow effects in MOLLI

Changes - outcome measures:
- Outcome marker initial values rather than average over the experiment. 
  The drug administration is timed for maximum inhibition at the 
  start + the timing between 
  scans is variable and therefore also the mean k's.

"""

import os
import argparse

import tristan_rifampicin_effect
import tristan_rifampicin_clinical_effect
import tristan_metformin_effect
import tristan_ciclosporin_effect
import methods.report


if __name__ == '__main__':

    BASE = 'C:\\Users\\md1spsx\\Documents\\GitHub\\tristan-human-stage-2-modelling'
    DATA = os.path.join(BASE, 'data')
    BUILD = os.path.join(BASE, 'wip', 'build')

    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=str, default=DATA, help="Data folder")
    parser.add_argument("--build", type=str, default=BUILD, help="Build folder")
    args = parser.parse_args()

    for scans in [1,2]:
        tristan_rifampicin_clinical_effect.main(args.data, args.build, scans=scans)
        tristan_rifampicin_effect.main(args.data, args.build, scans=scans)
        tristan_ciclosporin_effect.main(args.data, args.build, scans=scans)
        tristan_metformin_effect.main(args.data, args.build, scans=scans)
        
    methods.report.all_drugs(args.build)
