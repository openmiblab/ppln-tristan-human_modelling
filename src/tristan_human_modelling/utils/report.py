import os

import miblab

DRUGS = ['rifampicin', 'ciclosporin', 'metformin']


def two_scans(        
        resultspath, 
        outputpath,
        title = 'Liver-mediated DDI study',
        subtitle = 'Two scans',
        subject = 'Internal report'
    ):

    print('Creating report..')

    doc = miblab.Report(
        outputpath,
        'two_scans',
        title = title,
        subtitle = subtitle,
        subject = subject,
    )

    doc.chapter('Numerical outcomes')

    doc.section('Rate constants')

    for drug in DRUGS + ['patients_rifampicin']:
        table = os.path.join(resultspath, 'stage_5_tables', drug, 'scans_2', 'outcomes.csv')
        caption = f"Primary outcomes for the drug {drug}."
        doc.table(table, caption=caption)

    doc.section('Changes during the visit', clearpage=True)

    for drug in DRUGS + ['patients_rifampicin']:
        table = os.path.join(resultspath, 'stage_5_tables', drug, 'scans_2', 'outcomes_diurnal.csv')
        caption = f"Secondary outcomes for the drug {drug}."
        doc.table(table, caption=caption) 

    doc.section('Constants', clearpage=True)

    for drug in DRUGS + ['patients_rifampicin']:
        table = os.path.join(resultspath, 'stage_5_tables', drug, 'scans_2', 'constants.csv')
        caption = f"Constants for the drug {drug}."
        doc.table(table, caption=caption) 

    doc.chapter('Plots of outcomes')

    doc.section('Effect plots', clearpage=True)

    for drug in DRUGS + ['patients_rifampicin']:
        fig = os.path.join(resultspath, 'stage_6_figs', drug, 'scans_2', '_effect_plot.png')
        caption = f"Effect of the drug {drug}"
        doc.figure(fig, width='6in', caption=caption)

    doc.section('Diurnal variations', clearpage=True)

    for drug in DRUGS + ['patients_rifampicin']:
        fig = os.path.join(resultspath, 'stage_6_figs', drug, 'scans_2', '_diurnal_function.png')
        caption = f"Effect of the drug {drug}"
        doc.figure(fig, width='7in', caption=caption)

    doc.chapter('Kinetics and model fits')

    for drug in DRUGS + ['patients_rifampicin']:
        doc.section(drug, clearpage=True)

        folder = os.path.join(resultspath, 'stage_2_modelfit', drug, 'scans_2', 'Plots')
        files = [f for f in os.listdir(folder)]
        for f in files:
            fig = os.path.join(folder, f)
            caption = f"Volunteer {f}, drug {drug}"
            doc.figure(fig, width='7in', caption=caption)

    doc.build()


def one_scan(        
        resultspath, 
        outputpath,
        title = 'Liver-mediated DDI study',
        subtitle = 'One scan',
        subject = 'Internal report'
    ):

    print('Creating report..')

    doc = miblab.Report(
        outputpath,
        'one_scan',
        title = title,
        subtitle = subtitle,
        subject = subject,
    )


    doc.chapter('Numerical outcomes')

    doc.section('Rate constants')

    for drug in DRUGS + ['patients_rifampicin']:
        table = os.path.join(resultspath, 'stage_5_tables', drug, 'scans_1', 'outcomes.csv')
        caption = f"Primary outcomes for the drug {drug}."
        doc.table(table, caption=caption) 

    doc.section('Constants', clearpage=True)

    for drug in DRUGS + ['patients_rifampicin']:
        table = os.path.join(resultspath, 'stage_5_tables', drug, 'scans_1', 'constants.csv')
        caption = f"Constants for the drug {drug}."
        doc.table(table, caption=caption) 

    doc.chapter('Plots of outcomes')

    doc.section('Effect plots', clearpage=True)

    for drug in DRUGS + ['patients_rifampicin']:
        fig = os.path.join(resultspath, 'stage_6_figs', drug, 'scans_1', '_effect_plot.png')
        caption = f"Effect of the drug {drug}"
        doc.figure(fig, width='6in', caption=caption)

    doc.chapter('Kinetics and model fits')

    for drug in DRUGS + ['patients_rifampicin']:
        doc.section(drug, clearpage=True)

        folder = os.path.join(resultspath, 'stage_2_modelfit', drug, 'scans_1', 'Plots')
        files = [f for f in os.listdir(folder)]
        for f in files:
            fig = os.path.join(folder, f)
            caption = f"Volunteer {f}, drug {drug}"
            doc.figure(fig, width='7in', caption=caption)

    doc.build()
