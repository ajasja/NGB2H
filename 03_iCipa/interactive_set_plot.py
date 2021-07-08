import os
import hepran
import hepran.bzipscore as bz
import hepran.bcipa as bc
import hepran.utils as u
from hepran.registers import *
import hepran.agadir as ag


import bokeh
import bokeh.plotting as bp
from bokeh.models import HoverTool

def get_formated_seq(row):
    row['seq1_disp'] = hepran.registers.html_format_amino_acids(row.seq1,'f')
    row['seq2_disp'] = hepran.registers.html_format_amino_acids(row.seq2,'f')
    return row

tooltips = [
    ('ID1, ID2', '@IDs'),
    ('Tm', '@Tm'),
    ('score', '@score'),
    ('seq1', '@seq1_disp{safe}'),
    ('seq2', '@seq2_disp{safe}'),
]



def draw_scatter_interactive(x, y, data, print_sum=False, title="",
                             y_range=None,
                             tooltips=tooltips,
                             save_to_file=False, width=400, height=400):
    """Draws a regression plot with r2 given in title"""

    on_target = data[data['on_target'] == True]
    off_t = data[data['on_target'] == False]

    p = bp.figure(title=title, tools="pan,wheel_zoom,box_zoom,reset,save",
                  y_range=y_range, width=width, height=height)
    p.toolbar.active_drag = None
    p.add_tools(HoverTool(tooltips=tooltips))

    #p.line(x1, res.fittedvalues, color="black", line_width=3)
    #p.line(x1, iv_u, color="gray", line_width=3, alpha=0.5)
    #p.line(x1, iv_l, color="gray", line_width=3, alpha=0.5)

    p.circle(x, y, fill_alpha=0.2, size=8, source=off_t)
    p.circle(x, y, color="red", fill_alpha=0.2, size=8, source=on_target)

    if y_range:
        p.line([y_range[0], y_range[1]], [y_range[0], y_range[1]],
               line_width=2, line_dash="dashed", line_color="black", line_alpha=0.5)


    p.xaxis.axis_label = x
    p.yaxis.axis_label = y

    #if save_to_file:
    #    print("SAVING")
    #    bp.output_file(title.replace(" ", "_") + '.html',  title)


    return p


def plot_orto_set_interactive(df, ids, title="", pairs=None, field="mean_RD", out_name=None, lower_is_better=False,
                  vmin=None, vmax=None, ax=None):
    if pairs is None:
        pairs = list(u.chunks(ids, 2))
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=[12, 8])

    ar = series_to_array(df[field], ids=ids)
    OI = bs.get_ortogonality_index(ar, aids=ids, pairs=pairs, lower_is_better=lower_is_better, normalized=False)
    GAP = bs.get_ortogonality_gap(ar, lower_is_better=lower_is_better, aids=ids, pairs=pairs, kind='minmax', normalized=False)
    NAM = bs.get_num_offtarget_above(ar, lower_is_better=lower_is_better, aids=ids, pairs=pairs, kind='min')

    atitle = "{}, {}\n(OI={:2.0f}, GAP={:2.1f}, NAM={:d})".format(title, field, OI, GAP, NAM)

    if lower_is_better:
        cmap = 'OrRd_r'
    else:
        cmap = 'OrRd'
    ax = bs.plot_matrix_mpl(mat=ar, aids=ids, cmap=cmap, title=atitle, ax=ax, vmin=vmin, vmax=vmax)
    # ax.figure.savefig("out/"+title)
    # fig.tight_layout()

    if not out_name is None:
        out_name = out_name.format(title, field)
        ax.figure.savefig(out_name)
    return ax.figure