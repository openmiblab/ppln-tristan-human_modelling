import os

from methods import calc
import methods.plot


def run(dataset, datapath, results,
        effect_range=([-100,200], [-100,500]),
        k_max = [30,5],
    ):

    if not os.path.exists(results):
        return

    # calc.build_results(path)
    # calc.descriptive_statistics(path)
    # calc.averages(path)
    # calc.pairwise_ttest(path)

    # plot.create_bar_chart(path)
    methods.plot.effect_plot_combined(results, ylim=k_max)
    methods.plot.diurnal_k(dataset, datapath, results, ylim=k_max)

    # tables.averages(path)
    # tables.pairwise_stats(path)
    # tables.cases(path)   