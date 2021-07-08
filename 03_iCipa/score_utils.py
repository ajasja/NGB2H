from __future__ import division, absolute_import, print_function

import os

import hepran
import hepran.bzipscore as bz
import hepran.bcipa as bc
import hepran.utils as u
import hepran.registers as r
import hepran.agadir as ag

import pandas as pd
import numpy as np

import sklearn as sk
import sklearn
from sklearn import linear_model
from sklearn.externals import joblib

import matplotlib.pylab as plt
import seaborn as sns

from interactive_set_plot import *

#seq1 = 'DKNAALKAENAALEYEIAALEAEIAALEG'
#seq2 = 'DKNAALKAENAALEYEIAALEAEIAALEG'
#core_p, es_p = r.get_core_and_es_pairs(seq1, seq2, 'f', 'f')


def canonize(cp):
    """Returns the minimum lexial representation of string or reverse string
    ABC
    and
    CBA

    are equal, but this function returns always ABC
    """
    rev = cp[::-1]
    return min(rev, cp)


assert canonize("ABC") == "ABC"
assert canonize("CBA") == "ABC"


def get_crosspairs(core_p, es_p):
    N = len(core_p)
    res = []
    for i in range(N):
        res.append(
            canonize("".join(es_p[i]) + "".join(core_p[i]) + "".join(es_p[i + N])))
    return res


#get_crosspairs(core_p, es_p)

def get_CC_features(row):
    seq1, seq2 = r.trim_pair_seq(row.seq1, row.seq2)

    core_p, es_p = r.get_core_and_es_pairs(seq1, seq2, 'f', 'f')

    nterm = r.pairs_counter([core_p[0]])
    cterm = r.pairs_counter([core_p[-1]])
    # basically do a one hot nterm and cterm encoder
    for p in "NN IN II".split():
        row['nterm_c_' + p] = nterm.get(p, 0)

    for p in "NN IN II".split():
        row['cterm_c_' + p] = cterm.get(p, 0)

    core = r.pairs_counter(core_p)
    es = r.pairs_counter(es_p)
    row['c_NN'] = core['NN']
    row['c_IN'] = core['IN']
    row['c_II'] = core['II']
    #hacked because all are L currently
    row['cd_LL'] = len(seq1)//7
    row['es_EE'] = es['EE']
    row['es_EK'] = es['EK']
    row['es_KK'] = es['KK']


    #should the charges be trimmed?
    row['interface_charge1'] = round(r.get_interface_charge(row.seq1, 'f'), 1)
    row['interface_charge2'] = round(r.get_interface_charge(row.seq2, 'f'), 1)
    row['interface_repulsion'] = row['interface_charge1'] * row['interface_charge2']

    # vertical pairs
    # cv = core vertical
    cv = r.pairs_counter(r.get_vertical_pairs(core_p), sort=False)
    for p in "NN II NI IN".split():
        row['cv_' + p] = cv.get(p, 0)

    # esv = electrostatic vertical
    esv = r.pairs_counter(r.get_vertical_pairs(es_p), sort=False)
    for p in "KK EE EK KE".split():
        row['esv_' + p] = esv.get(p, 0)

    # cvt = core vertical triplets
    for p in "KK EE EK KE".split():
        row['esv_' + p] = esv.get(p, 0)

    cvt = r.pairs_counter(r.get_vertical_triplets(core_p), sort=False)
    for t in r.get_all_combinations("IN", 3):
        t = ''.join(t)
        row['cvt_' + t] = cvt.get(t, 0)

    esvt = r.pairs_counter(r.get_vertical_triplets(es_p), sort=False)
    for t in r.get_all_combinations("KE", 3):
        t = ''.join(t)
        row['esvt_' + t] = esvt.get(t, 0)

    cp = collections.Counter(get_crosspairs(core_p, es_p))
    for p in ['EEIIEE','EEINEE','EENNEE','EKIIEK','EKINEK','EKNIEK','EKNNEK','KKIIKK','KKINKK','KKNNKK']:
        row['cp_'+p] = cp.get(p, 0)

    return row


def get_CC_features_opt(row, seqs, pos=0):
    """Get just the essential features"""

    row['pos'] = pos

    seq1, seq2 = r.trim_pair_by_pos(seqs[row["ID1"]], seqs[row["ID2"]], pos)

    core_p, es_p = r.get_core_and_es_pairs(seq1, seq2, 'f', 'f')

    nterm = r.pairs_counter([core_p[0]])
    cterm = r.pairs_counter([core_p[-1]])
    # basically do a one hot nterm and cterm encoder
    for p in "NN IN II".split():
        row['nterm_c_' + p] = nterm.get(p, 0)

    for p in "NN IN II".split():
        row['cterm_c_' + p] = cterm.get(p, 0)

    core = r.pairs_counter(core_p)
    es = r.pairs_counter(es_p)
    row['c_NN'] = core['NN']
    row['c_IN'] = core['IN']
    row['c_II'] = core['II']
    #hacked because all are L currently
    row['cd_LL'] = len(seq1)//7
    row['es_EE'] = es['EE']
    row['es_EK'] = es['EK']
    row['es_KK'] = es['KK']


    #should the charges be trimmed?
    #row['interface_charge1'] = round(r.get_interface_charge(row.seq1, 'f'), 1)
    #row['interface_charge2'] = round(r.get_interface_charge(row.seq2, 'f'), 1)
    row['interface_repulsion'] =  round(r.get_interface_charge(seqs[row["ID1"]], 'f'), 1) * \
                                  round(r.get_interface_charge(seqs[row["ID2"]], 'f'), 1)

    # vertical pairs
    # cv = core vertical
    cv = r.pairs_counter(r.get_vertical_pairs(core_p), sort=False)
    for p in "NN II NI IN".split():
        row['cv_' + p] = cv.get(p, 0)

    # esv = electrostatic vertical
    esv = r.pairs_counter(r.get_vertical_pairs(es_p), sort=False)
    for p in "KK EE EK KE".split():
        row['esv_' + p] = esv.get(p, 0)

    return row

def score_opt(feat_dict, model, positons):
    """Given a linear model a list of positions performe the scoring"""

def make_model(target, fit_fields, Q, type="Ridge",):
    """Gets passed a dataset and makes a linear model, The dataset must also have a weights field"""

    Y = Q[target]
    X = Q[fit_fields]
    W = Q.weights

    if type== "LinearRegression":
        lm = linear_model.LinearRegression(normalize=False)
        lm.fit(X, Y, W)
    if type == "Ridge":
        lm = linear_model.Ridge(normalize=False, alpha=0.0003)
        lm.fit(X, Y, W)

    if type == "ElasticNet":
        lm = linear_model.ElasticNet(normalize=False, alpha=0.03, l1_ratio=.5)
        lm.fit(X, Y)
    if type == "SGDRegressor":
        lm = linear_model.SGDRegressor(alpha=0.03, l1_ratio=.15, max_iter=100)
        lm.fit(X, Y)
    if type == "BayesianRidge":
        #lm = linear_model.Ridge(normalize=False, alpha=0.03)
        lm = linear_model.BayesianRidge(normalize=False)
        lm.fit(X, Y)

    R2 = lm.score(X, Y, W)
    return lm, R2

def getRMSE(x, y, df):
    return np.sqrt(((df[x]-df[y])**2).mean())



from scipy.stats.stats import pearsonr

def getPearsonR(x, y, df):
    return pearsonr(df[x], df[y])[0]

def print_RMSE(df):
    print("RMSE all: ", getRMSE('Tm', 'score', df))
    print("RMSE >25: ", getRMSE('Tm', 'score', df.query('Tm > 25')))
    print("RMSE on_target: ", getRMSE('Tm', 'score', df.query('on_target')))

def get_model_features(lm, fit_fields, N_iter, fit_type):
    """Returns the model coeficients as pandas dataframe"""
    df = pd.DataFrame(list(zip(fit_fields, lm.coef_)), columns="feature coef".split())
    df['N_iter'] = N_iter
    df['fit_type'] = fit_type
    df.set_index('N_iter fit_type'.split(), inplace=True)
    return df


def get_metrics_df():
    return pd.DataFrame(
        columns="N_iter fit_type fit_class corrR R2_score RMSE med_abs_err explained_var Baysian_IC Akaike_IC N_samples N_feat".split())

def get_features_df():
    return pd.DataFrame(
        columns="N_iter fit_type feature coef".split())

# x=real, y=pred
def get_FIT_metrics(x, y, w, N_iter, fit_type, fit_class, N_feat, N_samples):
    r = {}
    r['N_iter'] = N_iter
    r['fit_type'] = fit_type
    r['fit_class'] = fit_class
    r['N_feat'] = N_feat
    r['N_samples'] = N_samples

    r['corrR'] = pearsonr(x, y)[0]
    r['R2_score'] = sklearn.metrics.r2_score(x, y)
    r['RMSE'] = np.sqrt(sklearn.metrics.mean_squared_error(x, y, w))
    r['med_abs_err'] = np.sqrt(sklearn.metrics.median_absolute_error(x, y))
    r['explained_var'] = sklearn.metrics.explained_variance_score(x, y, w)
    r['Baysian_IC'] = N_samples * np.log(r['RMSE'] ** 2) + N_feat * np.log(N_samples)
    r['Akaike_IC'] = 2 * N_feat + N_samples * np.log(N_samples * r['RMSE'] ** 2)

    return r


def mpl_plot_fit(fit_type, df):

    plt.plot([0, 80], [0, 80], 'k--', lw=4)
    plt.title('iCipa ' + fit_type)
    plt.ylim([0, 80])
    plt.xlim([0, 80])

    plt.gca().set_aspect('equal')

    sns.regplot(x='Tm', y='score', fit_reg=True, data=df)
    sns.regplot(x='Tm', y='score', fit_reg=True, data=df.query('on_target'), color='red')



def rescore_with_new_pos(df, lm, fit_fields, weight_string, set_weights, lm_type, extra_cols=[]):
    # make the alignemetns
    LEN = len('EIAALEAENAALEYKNAALKAKIAALKG')
    # reserve memory in the apply phase for the features. More than 100x speedup
    feature_cols = "nterm_c_NN          nterm_c_IN          nterm_c_II          cterm_c_NN          cterm_c_IN          cterm_c_II          c_NN                c_IN                c_II                es_EE               es_EK               es_KK               interface_charge1   interface_charge2   interface_repulsion cv_NN               cv_II               cv_NI               cv_IN               esv_KK              esv_EE              esv_EK              esv_KE              cvt_NNI             cvt_IIN             cvt_NNN             cvt_NIN             cvt_INN             cvt_III             cvt_INI             cvt_NII             esvt_EKK            esvt_EEK            esvt_EEE            esvt_EKE            esvt_KEE            esvt_KKK            esvt_KKE            esvt_KEK            cd_LL".split() + \
                   'cp_EEIIEE cp_EEINEE cp_EENNEE cp_EKIIEK cp_EKINEK cp_EKNIEK cp_EKNNEK cp_KKIIKK cp_KKINKK cp_KKNNKK avg_HP_12'.split()
    seqs_00 = df['seq1 seq2 Tm on_target pos score'.split() + feature_cols + extra_cols].reset_index()
    seqs_00.rename(columns={'pos': 'old_pos', 'score': 'old_score'}, inplace=True)
    seqs_00['pos'] = 0
    seqs_00['len'] = LEN

    seqs_p07 = df['seq1 seq2 Tm on_target pos score'.split() + feature_cols + extra_cols].reset_index()
    seqs_p07.rename(columns={'pos': 'old_pos', 'score': 'old_score'}, inplace=True)
    seqs_p07['pos'] = 7
    seqs_p07['seq1'] = seqs_p07.seq1.astype(str) + '-------'
    seqs_p07['seq2'] = '-------' + seqs_p07.seq2.astype(str)
    seqs_p07['len'] = LEN - 7

    seqs_m07 = df['seq1 seq2 Tm on_target pos score'.split() + feature_cols + extra_cols].reset_index()
    seqs_m07.rename(columns={'pos': 'old_pos', 'score': 'old_score'}, inplace=True)
    seqs_m07['pos'] = -7
    seqs_m07['seq1'] = '-------' + seqs_m07.seq1.astype(str)
    seqs_m07['seq2'] = seqs_m07.seq2.astype(str) + '-------'
    seqs_m07['len'] = LEN - 7

    # merge all to
    seqs = pd.concat([seqs_00, seqs_p07, seqs_m07], ignore_index=True)


    ###Get the new features
    #This takes a really long time
    seqs = seqs.apply(get_CC_features, axis=1)
    seqs['score'] = lm.predict(seqs[fit_fields])

    ####find best position
    #get index position
    idx = seqs.groupby("ID1 ID2".split())['score'].idxmax()
    seqs = seqs.loc[idx]
    seqs.set_index("ID1 ID2".split(), inplace=True)

    seqs = seqs.apply(get_formated_seq, axis=1)
    seqs['IDs'] = seqs.index
    seqs['old_func_new_pos'] = seqs.score


    #get new features
    seqs = seqs.apply(get_CC_features, axis=1)

    ###make new model
    set_weights(weight_string, seqs)
    lm, R2 = make_model('Tm', fit_fields, seqs, lm_type)
    seqs['score'] = lm.predict(seqs[fit_fields])
    seqs = seqs.apply(get_formated_seq, axis=1)
    seqs['IDs'] = seqs.index


    #color seqs for display
    return seqs, lm, R2


def get_FIT_dataframe(df, lm, N_iter, fit_type, N_feat, N_samples):
    r = get_metrics_df()

    Q = df
    if len(Q)>0:
        r = r.append(
        get_FIT_metrics(Q.Tm, Q.score, Q.weights, N_iter, fit_type, fit_class='all', N_feat=len(lm.coef_) + 1,
                        N_samples=len(Q)),
        ignore_index=True)

    Q = df.query('Tm>25')
    if len(Q) > 0:
        r = r.append(
        get_FIT_metrics(Q.Tm, Q.score, Q.weights, N_iter, fit_type, fit_class='Tm>25', N_feat=len(lm.coef_) + 1,
                        N_samples=len(Q)),
        ignore_index=True)


    Q = df.query('Tm>55')
    if len(Q) > 0:
        r = r.append(
        get_FIT_metrics(Q.Tm, Q.score, Q.weights, N_iter, fit_type, fit_class='Tm>55', N_feat=len(lm.coef_) + 1,
                        N_samples=len(Q)),
        ignore_index=True)

    Q = df.query('on_target == True')
    if len(Q) > 0:
        r = r.append(
        get_FIT_metrics(Q.Tm, Q.score, Q.weights, N_iter, fit_type, fit_class='on_target', N_feat=len(lm.coef_) + 1,
                        N_samples=len(Q)),
        ignore_index=True)

    r = r.set_index('N_iter fit_type fit_class'.split())
    return r


def ln_mean_RD_to_Tm(RD):
    return 19.94604942*RD - 57.672797381849335



basic = 'c_NN c_IN c_II es_EE es_EK es_KK'.split()
basicL = 'c_NN c_IN c_II es_EE es_EK es_KK cd_LL'.split()
basicCP = 'cp_EEIIEE cp_EEINEE cp_EENNEE cp_EKIIEK cp_EKINEK cp_EKNIEK cp_EKNNEK cp_KKIIKK cp_KKINKK cp_KKNNKK'.split()
nterm_core = 'nterm_c_NN nterm_c_IN nterm_c_II'.split()
cterm_core = 'cterm_c_NN cterm_c_IN cterm_c_II'.split()
core_ends = nterm_core + cterm_core
hel = ['avg_HP_12']

core_vertical = 'cv_NN cv_II cv_NI cv_IN'.split()
es_vertical = 'esv_KK esv_EE esv_EK esv_KE'.split()

core_t_vertical = 'cvt_INI cvt_IIN cvt_NII cvt_NNI cvt_NNN cvt_NIN cvt_INN cvt_III'.split()
es_t_vertical = 'esvt_EKK esvt_EEK esvt_EEE esvt_EKE esvt_KEE esvt_KKK esvt_KKE esvt_KEK'.split()
interface_rep = ['interface_repulsion']

fit_fields_dic = {
    "basic": basic,
    "basic-ends": basic + core_ends,
    "basic-rep": basic + interface_rep,
    "basic-core_vertical": basic + core_vertical,
    "basic-ends_rep": basic + core_ends + interface_rep,
    "basic-rep-core_vertical": basic + interface_rep + core_vertical,
    "basic-rep-nter_core-vertical_t_all": basic + interface_rep + nterm_core + core_t_vertical + es_t_vertical,
    "basic-rep-core_t_vertical": basic + interface_rep +  core_t_vertical,
    "basic-rep-nter_core-core_vertical": basic + interface_rep + nterm_core + core_vertical,
    "basic-rep-core_vertical": basic + interface_rep + core_vertical,
    "basic-rep-es_t_vertical": basic + interface_rep + es_t_vertical,
    "basic-rep-nter_core": basic + interface_rep + nterm_core,
    "basicL-rep-hel": basicL + interface_rep + hel, 
    "basicL": basicL,
    "basicL-ends": basicL + core_ends,
    "basicL-rep": basicL + interface_rep,
    "basicL-core_vertical": basicL + core_vertical,
    "basicL-ends_rep": basicL + core_ends + interface_rep,
    "basicL-rep-nter_core-vertical_t_all": basicL + interface_rep + nterm_core + core_t_vertical + es_t_vertical,
    "basicL-rep-core_t_vertical": basicL + interface_rep +  core_t_vertical,
    "basicL-rep-nter_core-core_vertical": basicL + interface_rep + nterm_core + core_vertical,
    "basicL-rep-es_t_vertical": basicL + interface_rep + es_t_vertical,
    "basicL-rep-core_vertical": basicL + interface_rep + core_vertical,
    "basicL-rep-es_vertical": basicL + interface_rep +  es_vertical,
    "basicL-rep-es_core_vertical": basicL + interface_rep +  es_vertical +core_vertical,
    "basicL-rep-nter_core": basicL + interface_rep + nterm_core,
    "basicL-rep-hel": basicL + interface_rep + hel,
    "basicL-nter_core": basicL + nterm_core,
    "basicCP": basicCP,
    "basicCP-ends": basicCP + core_ends,
    "basicCP-rep": basicCP + interface_rep,
    "basicCP-core_vertical": basicCP + core_vertical,
    "basicCP-ends_rep": basicCP + core_ends + interface_rep,
    "basicCP-rep-core_vertical": basicCP + interface_rep + core_vertical,
    "basicCP-rep-nter_core-vertical_t_all": basicCP + interface_rep + nterm_core + core_t_vertical + es_t_vertical,
    "basicCP-rep-core_t_vertical": basicCP + interface_rep +  core_t_vertical,
    "basicCP-rep-nter_core-core_vertical": basicCP + interface_rep + nterm_core + core_vertical,
    "basicCP-rep-core_vertical": basicCP + interface_rep + core_vertical,
    "basicCP-rep-es_t_vertical": basicCP + interface_rep + es_t_vertical,
    "basicCP-rep-nter_core": basicCP + interface_rep + nterm_core,
}


def score_overlap(df, lm, fit_fields, get_features, LEN=len('EIAALEAENAALEYKNAALKAKIAALKG')):
    # make the alignemetns

    # reserve memory in the apply phase for the features. More than 100x speedup
    seqs_00 = df.reset_index()
    seqs_00.rename(columns={'pos': 'old_pos', 'score': 'old_score'}, inplace=True)
    seqs_00['pos'] = 0
    seqs_00['len'] = LEN

    seqs_p07 = df.reset_index()
    seqs_p07.rename(columns={'pos': 'old_pos', 'score': 'old_score'}, inplace=True)
    seqs_p07['pos'] = 7
    seqs_p07['seq1'] = seqs_p07.seq1.astype(str) + '-------'
    seqs_p07['seq2'] = '-------' + seqs_p07.seq2.astype(str)
    seqs_p07['len'] = LEN - 7

    seqs_m07 = df.reset_index()
    seqs_m07.rename(columns={'pos': 'old_pos', 'score': 'old_score'}, inplace=True)
    seqs_m07['pos'] = -7
    seqs_m07['seq1'] = '-------' + seqs_m07.seq1.astype(str)
    seqs_m07['seq2'] = seqs_m07.seq2.astype(str) + '-------'
    seqs_m07['len'] = LEN - 7

    # merge all to
    seqs = pd.concat([seqs_00, seqs_p07, seqs_m07], ignore_index=True)

    ###Get the new features
    # This takes a really long time
    seqs = seqs.apply(get_features, axis=1)
    seqs['score'] = lm.predict(seqs[fit_fields])

    ####find best position
    # get index position
    idx = seqs.groupby("ID1 ID2".split())['score'].idxmax()
    seqs = seqs.loc[idx]
    seqs.set_index("ID1 ID2".split(), inplace=True)

    return seqs


def save_features_to_disk(dicts, i_part, out_dir, feature_function, positions, seqs):
    """Takes a list of dictionary with ID1 and ID2 and save a pandas dataframe to out_dir as ith part """
    #imports enable easier parallel distribution
    import hepran
    import hepran.bzipscore as bz
    import hepran.registers as r
    import pandas as pd
    import numpy as np
    from copy import copy

    res = []
    for i in np.arange(len(dicts)):
        for pos in positions:
            res.append(
                feature_function(copy(dicts[i]), seqs, pos)
            )

    # out_name = out_dir+"/"+str(i_part)+".feat.h5"
    out_name = "{out_dir}/{ipart:02d}.feat.h5".format(**locals())
    print(out_name)
    df = pd.DataFrame(data=res, columns=res[0].keys(), dtype=np.int16)
    df.to_hdf(out_name, 'feat', format='fixed', mode="w")

    return out_name


def score_from_features(feature_file, linear_model, fit_fields, fit_name, pos="ALL", out_dir='./scores',
                        lower_is_better=False, return_df=False):
    import pandas as pd
    import numpy as np
    import os

    out_path = out_dir + "/" + fit_name

    # create out directory if it does not exist
    try:
        os.makedirs(out_path);
    except:
        pass;

    # only take the number of the file
    base_name = os.path.basename(feature_file).split('.')[0]
    base_name = base_name + ".score.h5"
    out_name = out_path + "/" + base_name

    # load the file and do the scoring
    df = pd.read_hdf(feature_file)
    if pos == 'ALL':
        df['score'] = linear_model.predict(df[fit_fields])
        # sort by score and then drop duplicates.
        # MUCH MUCH faster than group by solution
        df.sort_values("ID1 ID2 score".split(), ascending=[True, True, lower_is_better], inplace=True)
        df.drop_duplicates("ID1 ID2".split(), keep='first', inplace=True)
        df.sort_index(inplace=True)
    else:
        # score only one position
        df = df.query('pos == @pos')
        df['score'] = linear_model.predict(df[fit_fields])

    df = df['ID1 ID2 score pos'.split()]
    # restrict types to save disk space.
    df.score = df.score.astype('float32');
    df.pos = df.pos.astype('int8');
    df.to_hdf(out_name, 'feat', format='fixed', mode="w")
    if return_df:
        return df
    else:
        return out_name


def assign_to_matrix(v, mat, pos_mat, ri):
    # pandas iloc and iat are EXTREMLY slow! 2.6 s vs 150 ms
    for n in np.arange(len(v)):
        i = ri[v[n, 0]]
        j = ri[v[n, 1]]
        # score is third column
        s = v[n, 2]
        p = v[n, 3]
        mat[i, j] = s
        mat[j, i] = s
        pos_mat[i, j] = p
        pos_mat[j, i] = p


# %time assign_to_matrix(v, mat, pos_mat, ri)


def score_dir_to_score_matrix(scores_dir, seqs_file, out_name=None):
    from glob import glob
    if out_name is None:
        out_name = scores_dir
    out_name_mat = scores_dir + ".bin"
    out_name_mat_pos = scores_dir + ".align.bin"

    scores_files = glob(scores_dir + '/*.score.h5')

    seqs = u.load_fasta(seqs_file)

    N = len(seqs)

    # allocat the matrices
    mat = np.zeros((N, N), dtype='float32')
    pos_mat = np.zeros((N, N), dtype='int8')

    # get lookuptable from ID to index
    ri = {seq_id: i for i, seq_id in enumerate(seqs.keys())}

    # assign to matrices
    for score_file in scores_files:
        v = pd.read_hdf(score_file).values
        assign_to_matrix(v, mat, pos_mat, ri)

        # save files
    mat.tofile(out_name_mat)
    pos_mat.tofile(out_name_mat_pos)

    return out_name_mat

def set_kernels_dir(CDIR):
    """Returns an instance to running cluster and sets the dir ot CDIR"""

    import ipyparallel as ipp
    rc = ipp.Client()
    dv = rc[:]
    lv = rc.load_balanced_view()

    #go to correct dirs
    dv.push({"CDIR": CDIR})
    dv.execute('import os; os.chdir(CDIR); os.getcwd()').get()
    return lv
    #%px import os
    #%px os.chdir(CDIR)