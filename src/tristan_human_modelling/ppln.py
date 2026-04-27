from miblab import pipe

import tristan_human_modelling as ppln


def run(build, logfile):
    ppln.stage_1_download.run(build, logfile)
    ppln.stage_2_modelfit.run(build, logfile)
    ppln.stage_3_concatenate.run(build, logfile)
    ppln.stage_4_stats.run(build, logfile)
    ppln.stage_5_tables.run(build, logfile)
    ppln.stage_6_figs.run(build, logfile)
    ppln.stage_7_report.run(build, logfile)


if __name__=='__main__':

    build = r"C:\Users\md1spsx\Documents\Data\tristan"
    pipe.run_ppln(run, build, 'human_modelling')