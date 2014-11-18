#! /usr/bin/env python2.7
import matplotlib as mpl


color_cn = (0.0, 0.5, 0.0)
color_mci = (0.8, 0.8, 0.0)
color_ad = (1.0, 0.0, 0.0)
progression_dict = {
    'red': ((0.0, color_cn[0], color_cn[0]), (0.5, color_mci[0], color_mci[0]), (1.0, color_ad[0], color_ad[0])),
    'green': ((0.0, color_cn[1], color_cn[1]), (0.5, color_mci[1], color_mci[1]), (1.0, color_ad[1], color_ad[1])),
    'blue': ((0.0, color_cn[2], color_cn[2]), (0.5, color_mci[2], color_mci[2]), (1.0, color_ad[2], color_ad[2]))
}
progression_cmap = mpl.colors.LinearSegmentedColormap('my_colormap', progression_dict)


def setup_axes(plt, ax, xgrid=True, ygrid=True, xspine=True, yspine=True):
    assert isinstance(ax, mpl.axes.Axes)

    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    if not xspine:
        ax.spines['bottom'].set_color('none')
    if not yspine:
        ax.spines['left'].set_color('none')

    if xgrid:
        ax.xaxis.grid(True, linestyle=':', which='major', color='grey', alpha=0.5)
    if ygrid:
        ax.yaxis.grid(True, linestyle=':', which='major', color='grey', alpha=0.5)
    ax.tick_params(axis='both', which='both', bottom='on', top='off', left='off', right='off')

    plt.setp(ax.get_yticklabels(), rotation=45, fontsize=11)
    plt.setp(ax.get_xticklabels(), fontsize=11)


def get_biomarker_string(biomarker):
    if biomarker == 'synth_cdrsb':
        return '$\mathcal{M}^{CDR-SB}$'
    elif biomarker == 'synth_mmse':
        return '$\mathcal{M}^{MMSE}$'
    elif biomarker == 'synth_hipp':
        return '$\mathcal{M}^{HV}$'
    elif biomarker == 'CDRSB':
        return 'CDR-SB'
    else:
        return biomarker


def set_boxplot_color(boxplot, index, color):
    box = boxplot['boxes'][index]
    box.set_facecolor(color + (0.2,))
    box.set_edgecolor(color)
    median = boxplot['medians'][index]
    median.set_color(color)
    cap = boxplot['caps'][2 * index]
    cap.set_color(color)
    cap = boxplot['caps'][2 * index + 1]
    cap.set_color(color)
    whisker = boxplot['whiskers'][2 * index]
    whisker.set_linestyle('-')
    whisker.set_color(color)
    whisker = boxplot['whiskers'][2 * index + 1]
    whisker.set_linestyle('-')
    whisker.set_color(color)
    flier = boxplot['fliers'][2 * index]
    flier.set_color(color)
    flier = boxplot['fliers'][2 * index + 1]
    flier.set_color(color)
