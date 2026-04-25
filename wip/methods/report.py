import os

import miblab

DRUGS = ['rifampicin', 'ciclosporin', 'metformin', 'rifampicin_clinical']


def all_drugs(        
        resultspath, 
        title = 'Liver-mediated DDI study',
        subtitle = 'Key results',
        subject = 'Internal report'
    ):

    print('Creating report..')

    doc = miblab.Report(
        resultspath,
        'all_drugs',
        title = title,
        subtitle = subtitle,
        subject = subject,
    )

    doc.chapter('Two scans')

    build = os.path.join(resultspath, 'build_scans_2')

    doc.section('Rates and effect sizes')

    for drug in DRUGS:

        table = os.path.join(build, drug, 'Tables', 'outcomes.csv')
        caption = (
            f"Primary outcomes for the drug {drug}."
        )
        doc.table(table, caption=caption) 

    doc.section('Changes during the visit', clearpage=True)

    for drug in DRUGS:

        table = os.path.join(build, drug, 'Tables', 'outcomes_diurnal.csv')
        caption = (
            f"Secondary outcomes for the drug {drug}."
        )
        doc.table(table, caption=caption) 

    doc.section('Constants', clearpage=True)

    for drug in DRUGS:
        
        table = os.path.join(build, drug, 'Tables', 'constants.csv')
        caption = (
            f"Constants for the drug {drug}."
        )
        doc.table(table, caption=caption) 

    doc.section('Effect plots', clearpage=True)

    for drug in DRUGS:
        fig = os.path.join(build, drug, 'Figures', '_effect_plot.png')
        caption = (
            f"Effect of the drug {drug}"
        )
        doc.figure(fig, width='6in', caption=caption)

    doc.section('Diurnal variations', clearpage=True)

    for drug in DRUGS:
        fig = os.path.join(build, drug, 'Figures', '_diurnal_function.png')
        caption = (
            f"Effect of the drug {drug}"
        )
        doc.figure(fig, width='7in', caption=caption)

    doc.section('Model fits', clearpage=True)

    for drug in DRUGS:

        doc.subsection(drug, clearpage=True)

        folder = os.path.join(build, drug, 'Plots')
        files = [f for f in os.listdir(folder)]
        for f in files:
            fig = os.path.join(folder, f)
            caption = f"Volunteer {f}, drug {drug}"
            doc.figure(fig, width='7in', caption=caption)



    doc.chapter('One scan')

    build = os.path.join(resultspath, 'build_scans_1')

    doc.section('Rates and effect sizes')

    for drug in DRUGS:
        table = os.path.join(build, drug, 'Tables', 'outcomes.csv')
        caption = (
            f"Primary outcomes for the drug {drug}."
        )
        doc.table(table, caption=caption) 

    doc.section('Constants', clearpage=True)

    for drug in DRUGS:
        
        table = os.path.join(build, drug, 'Tables', 'constants.csv')
        caption = (
            f"Constants for the drug {drug}."
        )
        doc.table(table, caption=caption) 

    doc.section('Effect plots', clearpage=True)

    for drug in DRUGS:
        fig = os.path.join(build, drug, 'Figures', '_effect_plot.png')
        caption = (
            f"Effect of the drug {drug}"
        )
        doc.figure(fig, width='6in', caption=caption)

    doc.section('Model fits', clearpage=True)

    for drug in DRUGS:

        doc.subsection(drug, clearpage=True)

        folder = os.path.join(build, drug, 'Plots')
        files = [f for f in os.listdir(folder)]
        for f in files:
            fig = os.path.join(folder, f)
            caption = f"Volunteer {f}, drug {drug}"
            doc.figure(fig, width='7in', caption=caption)

    doc.build()