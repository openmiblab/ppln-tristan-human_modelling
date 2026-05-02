import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use the non-interactive backend for generating figures
import matplotlib.pyplot as plt
import pydmr


sym = {
    'control': 'b-',
    'drug': 'r-',
}
mark = {
    1: 'o',
    2: 'v',
    3: '^',
    4: '<',
    5: '>',
    6: 's',
    7: 'p',
    8: '*',
    9: 'x',
    10: 'd',
    11: 'X',
    12: 'o',
    13: 'v',
    14: '^',
    15: '<',
    16: '>',
    17: 's',
    18: 'p',
    19: '*',
    20: 'x',
    21: 'd',
    22: 'X',
}


def color(index: float) -> str:
    cmap = plt.get_cmap("tab20")  
    rgba_color = cmap(index)  # Get RGBA tuple
    hex_color = "#{:02x}{:02x}{:02x}".format(
        int(rgba_color[0] * 255),
        int(rgba_color[1] * 255),
        int(rgba_color[2] * 255)
    )  # Convert to hex

    return hex_color



def effect_plot_combined(file, figpath, drug):
    fig, (ax0, ax1, ax2) = plt.subplots(
        1, 3, width_ratios=[2, 4, 4], figsize=(8,3)
    )
    fig.subplots_adjust(wspace=0.5)

    def _effect_box_plots_combined():

        pars = pydmr.read(file)['pars']
        subjects = sorted(list({k[0] for k in pars}))
        all_data = [
            [100 * pars[(s, drug, 'r_khe')] for s in subjects],
            [100 * pars[(s, drug, 'r_kbh')] for s in subjects],
        ]

        pars = ['khe', 'kbh']

        linewidth = 1.0
        fontsize=10
        for axis in ['top','bottom','left','right']:
            ax0.spines[axis].set_linewidth(linewidth)

        # box plot
        boxprops = dict(linestyle='-', linewidth=linewidth, color='black')
        medianprops = dict(linestyle='-', linewidth=linewidth, color='black')
        whiskerprops = dict(linestyle='-', linewidth=linewidth, color='black')
        capprops = dict(linestyle='-', linewidth=linewidth, color='black')
        flierprops = dict(marker='o', markerfacecolor='white', markersize=6,
                    markeredgecolor='black', markeredgewidth=linewidth)
        bplot = ax0.boxplot(
            all_data,
            whis = [2.5,97.5],
            capprops=capprops,
            flierprops=flierprops,
            whiskerprops=whiskerprops,
            medianprops=medianprops,
            boxprops=boxprops,
            widths=0.3,
            vert=True,  # vertical box alignment
            patch_artist=True,  # fill with color
            labels=pars)  # will be used to label x-ticks
        ax0.set_xticklabels(labels=pars, fontsize=fontsize)
        #ax0.set_yticklabels(labels=ax0.get_yticklabels(), fontsize=fontsize)
        ax0.tick_params(axis='y', labelsize=fontsize)

        # fill with colors
        for patch in bplot['boxes']:
            patch.set_facecolor('lightsteelblue')

        # adding horizontal grid line
        ax0.yaxis.grid(True)
        ax0.set_ylabel('Relative effect size (%)', fontsize=fontsize)

    def _line_plots_combined():
        ylim = [40,4]
        pars = pydmr.read(file)['pars']
        subjects = sorted(list({k[0] for k in pars}))
        study = list({k[1] for k in pars})[0]

        fontsize=10
        markersize=6
        ax1.set_title('Hepatocellular uptake rate', fontsize=fontsize, pad=10)
        ax1.set_ylabel('khe (mL/min/100mL)', fontsize=fontsize)
        ax1.set_ylim(0, ylim[0])
        ax1.tick_params(axis='x', labelsize=fontsize)
        ax1.tick_params(axis='y', labelsize=fontsize)
        ax2.set_title('Biliary excretion rate', fontsize=fontsize, pad=10)
        ax2.set_ylabel('kbh (mL/min/100mL)', fontsize=fontsize)
        ax2.set_ylim(0, ylim[1])
        ax2.tick_params(axis='x', labelsize=fontsize)
        ax2.tick_params(axis='y', labelsize=fontsize)
        
        for i, s in enumerate(subjects):
            x = ['control', 'treatment']
            khe = [6000 * pars[(s, study, 'c_khe')], 6000 * pars[(s, study, 'd_khe')]]
            kbh = [6000 * pars[(s, study, 'c_kbh')], 6000 * pars[(s, study, 'd_kbh')]]
            si = i/len(subjects)  
            ax1.plot(x, khe, '-', label=s, marker=mark[int(i+1)], 
                    markersize=markersize, color=color(si))
            ax2.plot(x, kbh, '-', label=s, marker=mark[int(i+1)], 
                    markersize=markersize, color=color(si))
        ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    _effect_box_plots_combined()
    _line_plots_combined()

    file = os.path.join(figpath, '_effect_plot.png')
    os.makedirs(figpath, exist_ok=True)
    plt.savefig(fname=file, bbox_inches='tight')
    plt.close()




def diurnal_k_combined(datafile, output_file, figpath, drug, ylim=[40,40,4,4]):
    fontsize=10
    titlesize=12
    markersize=6
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(8,8))
    fig.subplots_adjust(
        left=0.1,
        right=0.9,
        bottom=0.1,
        top = 0.9, 
        wspace=0.3,
        #hspace=1,
        )
    ax = {
        'c_he': ax1,
        'd_he': ax2,
        'c_bh': ax3,
        'd_bh': ax4,
    }
    ax1.set_title('Control', fontsize=titlesize)
    ax1.set_xlabel('Time of day (hrs)', fontsize=fontsize)
    ax1.set_ylabel('khe (mL/min/100mL)', fontsize=fontsize)
    ax1.set_ylim(0, ylim[0])
    ax1.tick_params(axis='x', labelsize=fontsize)
    ax1.tick_params(axis='y', labelsize=fontsize)
    ax2.set_title('Treatment', fontsize=titlesize)
    ax2.set_xlabel('Time of day (hrs)', fontsize=fontsize)
    ax2.set_ylabel('khe (mL/min/100mL)', fontsize=fontsize)
    ax2.set_ylim(0, ylim[1])
    ax2.tick_params(axis='x', labelsize=fontsize)
    ax2.tick_params(axis='y', labelsize=fontsize)
    #ax3.set_title('Baseline', fontsize=titlesize)
    ax3.set_xlabel('Time of day (hrs)', fontsize=fontsize)
    ax3.set_ylabel('kbh (mL/min/100mL)', fontsize=fontsize)
    ax3.set_ylim(0, ylim[2])
    ax3.tick_params(axis='x', labelsize=fontsize)
    ax3.tick_params(axis='y', labelsize=fontsize)
    #ax4.set_title('Rifampicin', fontsize=titlesize)
    ax4.set_xlabel('Time of day (hrs)', fontsize=fontsize)
    ax4.set_ylabel('kbh (mL/min/100mL)', fontsize=fontsize)
    ax4.set_ylim(0, ylim[3])
    ax4.tick_params(axis='x', labelsize=fontsize)
    ax4.tick_params(axis='y', labelsize=fontsize)

    # Create line plots
    data = pydmr.read(datafile, format='nest')
    pars = pydmr.read(output_file)['pars']
    subjects = sorted(list({k[0] for k in pars}))
    for i, s in enumerate(subjects):

        # Simulated duration for the subject
        rois = data['rois'][s]
        time_2_control = rois['control']['time_2'] - rois['control']['time_1'][0]
        time_2_drug = rois['drug']['time_2'] - rois['drug']['time_1'][0]
        tmax = [max(time_2_control), max(time_2_drug)]

        for visit in ['c', 'd']:

            # Acquistion times for the visit
            t0 = rois['control']['time_1'][0] if visit=='c' else rois['drug']['time_1'][0]
            tm = tmax[0] if visit=='c' else tmax[1]
            t = np.arange(t0, t0 + tm, 1)

            # Skip if constant (one scan protocol)
            if (s, drug, f'{visit}_khe_i') not in pars:
                return
            
            # khe
            # Interpolate parameters on acquisition time
            khe_i = pars[(s, drug, f'{visit}_khe_i')]
            khe_f = pars[(s, drug, f'{visit}_khe_f')]
            khe = np.interp(t, [t[0], t[-1]], [khe_i, khe_f])              

            # Plot parameters
            si = i/len(subjects)
            ax[f"{visit}_he"].plot(
                t / 3600, khe * 6000, '-', 
                label=s, # marker=mark[int(i+1)], 
                markersize=markersize, color=color(si),
            )

            # kbh
            if visit=='c':
                kbh = t * 0 + pars[(s, drug, f'{visit}_kbh')]
            elif visit=='d':
                # vh = pars[(s, drug, f'vh')]
                # Th_i = vh / pars[(s, drug, f'{visit}_kbh_i')]
                # Th_f = vh / pars[(s, drug, f'{visit}_kbh_f')]
                # Th = np.interp(t, [t[0], t[-1]], [Th_i, Th_f])
                # kbh = vh / Th 
                kbh = t * 0 + pars[(s, drug, f'{visit}_kbh')]
            
            # Plot parameters
            si = i/len(subjects)
            ax[f"{visit}_bh"].plot(
                t / 3600, kbh * 6000, '-', 
                label=s, # marker=mark[int(i+1)], 
                markersize=markersize, color=color(si),
            )
                
    ax['d_he'].legend()

    plot_file = os.path.join(figpath, '_diurnal_function.png')
    os.makedirs(figpath, exist_ok=True)
    plt.savefig(fname=plot_file)
    plt.close()
