import hepran
import hepran.bzipscore as bs
from hepran.draw_matrix import plot_matrix_mpl
import hepran.utils as u
import matplotlib as mpl
import matplotlib.pylab as plt
import numpy as np
import pandas as pd

def plot_pairs(df, ids, title="", pairs=None, field="mean_RD", out_name="exp_img/{}, {}.png", bcf_set=None):
    if pairs is None:
        pairs = list(u.chunks(ids,2))
    fig, axs = plt.subplots(1,2,figsize=[12,8])
    for n,condition in enumerate((0,6)):
        ar = series_to_array(df[field], ids=ids, condition=condition, bcf_set=bcf_set)

        OI = bs.get_ortogonality_index(ar, aids=ids, pairs=pairs)
        GAP = bs.get_ortogonality_gap(ar, lower_is_better=False, aids=ids, pairs=pairs, kind='minmax')
        NAM = bs.get_num_offtarget_above(ar, lower_is_better=False, aids=ids, pairs=pairs, kind='min')
        
        
        atitle="{}, {}, induction {} \n(OI={:2.0f}, GAP={:2.1f}, NAM={:d})".format(title, field, condition, OI, GAP, NAM)

        ax = bs.plot_matrix_mpl(mat=ar, aids=ids, cmap='OrRd', title=atitle, ax=axs[n])
        #ax.figure.savefig("out/"+title)
    #fig.tight_layout()    
    
    if not out_name is None:
        out_name = out_name.format(title, field)
        fig.savefig(out_name)
    return fig

def normalize_orto_set(mat, pairs, aids, lower_is_better=False, kind='med', normalize_zero=True, scale=100, q_range=(1,95)):
    """
    normalizes an orto set set to the average/median of on-taget and off-target sets.
    """
    
    on_target, off_target = bs.get_on_and_off_target_values(mat, pairs, aids)   
    
    if lower_is_better:
        mat = mat * -1
        
    if   kind.lower()=='med':
        max_val = np.nanmedian(on_target)
        min_val = np.nanmedian(off_target)
    elif (kind.lower() == 'avg') or (kind.lower() == 'mean'):
        max_val = np.nanmean(on_target)
        min_val = np.nanmean(off_target)
    elif kind.lower() == 'minmax':
        max_val = np.nanmax(on_target)
        min_val = np.nanmin(off_target)        
    elif kind.lower() == 'per':        
        max_val = np.percentile(mat, q_range[1])
        min_val = np.percentile(mat, q_range[0])       
    else:
        raise ValueError('Kind must be either avg or med')
    
    if normalize_zero:
        mat = mat - min_val
        mat = mat / (max_val-min_val)   
    else:
        mat = mat / max_val
    
    return mat*scale


def series_to_array_condition(df, ids=None, symetric=True, bcf_set=None, condition=6, method=np.sum):
    """From an indexd dataframe return an numpy array"""
    if ids is None:
        pass #implement
    res = np.zeros((len(ids), len(ids)))
    for ni, i in enumerate(ids):
        for nj, j in enumerate(ids):
            try:
                #print(i,j)
                #print(df.loc[i,j,bcf_set,condition])
                res[ni,nj]=df.loc[i,j,bcf_set,condition]
                if symetric:
                    res[nj,ni]=res[ni,nj]
            except (KeyError, IndexError) as e:
                if not symetric:
                    res[ni,nj]=np.nan
    return res    

def series_to_array(df, ids=None, symetric=True, extra_ids=()):
  """From an indexd dataframe return an numpy array"""
  if ids is None:
      pass #implement
  res = np.zeros((len(ids), len(ids)))
  for ni, i in enumerate(ids):
      for nj, j in enumerate(ids):
          try:
              #print(i,j)
              #print(df.loc[(i,j)+extra_ids])
              res[ni,nj]=df.loc[(i,j)+extra_ids]
              if symetric:
                  res[nj,ni]=res[ni,nj]
          #except (KeyError, IndexError) as e:
          except:
              if not symetric:
                  res[ni,nj]=np.nan
  return res    

def plot_orto_set(df, ids, title="", pairs=None, field="mean_RD", out_name=None, lower_is_better=False,
              vmin=None, vmax=None, ax=None, symetric=True, outline=True, make_symmetric=False):
    if pairs is None:
        pairs = list(u.chunks(ids,2))
    if ax is None:
        fig, ax = plt.subplots(1,1,figsize=[12,8])
    
    ar = series_to_array(df[field], ids=ids, symetric=symetric)
    if make_symmetric:
        #print("Making smmetric")
        ar =  (ar + ar.T )/2
    OI = bs.get_ortogonality_index(ar, aids=ids, pairs=pairs, lower_is_better=lower_is_better)
    GAP = bs.get_ortogonality_gap(ar, lower_is_better=lower_is_better, aids=ids, pairs=pairs, kind='minmax')
    NAM = bs.get_num_offtarget_above(ar, lower_is_better=lower_is_better, aids=ids, pairs=pairs, kind='min')


    atitle="{}, {}\n(OI={:2.0f}, GAP={:2.1f}, NAM={:d})".format(title, field,  OI, GAP, NAM)

    if lower_is_better:
        cmap = 'viridis'
    else:
        cmap= 'viridis_r'
    
    
    if outline:
        import matplotlib
        if ax is None:
            ax = plt.gca()
        for i, id1 in enumerate(ids):
            for j, id2 in enumerate(ids):
                if ((id1, id2) in pairs) or ((id2, id1) in pairs):
                    #rect= matplotlib.patches.Rectangle((i-0.5,j-0.5), 1, 1, edgecolor="white", fill=0, linewidth=3)
                    rect= matplotlib.patches.Circle((i,j), 0.08, color="black", fill=1, linewidth=3)
                    ax.add_patch(rect)

    fig = bs.plot_matrix_mpl(mat=ar, aids=ids, cmap=cmap, title=atitle, ax=ax, vmin=vmin, vmax=vmax, out_name=out_name)
    #ax.figure.savefig("out/"+title)
    #fig.tight_layout()    
    
   
    return fig
