#!/usr/bin/env python3
import sys, os
from pandas.io.parsers import read_csv
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage
import itertools as it
import PyComplexHeatmap as pch
from scipy.stats import ttest_ind
from statsmodels.sandbox.stats.multicomp import multipletests

outdir = './pca_results'
os.system('mkdir -p '+outdir)

vst_df = read_csv('deseq2_vst_data.tsv', sep = '\t', index_col = 0)
norm_df = read_csv('deseq2_norm_data.tsv', sep = '\t', index_col = 0)
log_df = np.log(norm_df+1)

def get_pca(X, nc = 55):
    pca = PCA(n_components = nc).fit(X.T)
    rpca = []
    for i in range(20):
        RX = X.copy(); RX = np.array(RX)
        for x in RX:
            np.random.shuffle(x)
        RX = pd.DataFrame(RX)
        rpca.append(PCA(n_components = nc).fit(RX.T))
    rand_var_expl = pd.DataFrame([x.explained_variance_ratio_ for x in rpca])
    return pca, rand_var_expl

vst_pca, vst_rnd_pcavar = get_pca(vst_df, nc = 5)
norm_pca, norm_rnd_pcavar = get_pca(norm_df, nc = 5)
log_pca, log_rnd_pcavar = get_pca(log_df, nc = 5)

def plot_varianceratio(pca, rnd_pcavar, ax):
    ax.scatter(range(len(rnd_pcavar.columns)), pca.explained_variance_ratio_, label = 'signal')
    ax.scatter(range(len(rnd_pcavar.columns)), rnd_pcavar.mean(), label = 'randomized')
    ax.legend()
    return ax

fig, axs = plt.subplots(ncols = 3, figsize = (3*3*1.6, 3))
axs[0] = plot_varianceratio(vst_pca, vst_rnd_pcavar, axs[0])
axs[1] = plot_varianceratio(norm_pca, norm_rnd_pcavar, axs[1])
axs[2] = plot_varianceratio(log_pca, log_rnd_pcavar, axs[2])
axs[0].set_xlabel("PC"); axs[1].set_ylabel("PC");
axs[0].set_ylabel('explained variance ratio')
axs[0].set_title("VST normalized"); axs[1].set_title("Size factor normalized"); axs[2].set_title('Log-norm')
fig.savefig(outdir + '/PCA_variance.pdf', bbox_inches = 'tight')

df = pd.DataFrame(vst_pca.fit_transform(log_df.T), index = log_df.columns).T # pca components per sample
df.to_csv(outdir + '/PCA_coordinates.tsv', sep = '\t')

fig, ax = plt.subplots()
for label,c in zip(['ND','T2D'],['r','b']):
    cols = [c for c in df.columns if label in c]
    ax.scatter(df.loc[0,cols], df.loc[1,cols], c = c, label = label)
ax.set_xlabel('PC 1'); ax.set_ylabel('PC 2')
for c in df.columns:
    ax.text(df.loc[0,c], df.loc[1,c], '  '+c, va = 'center', ha = 'left', fontsize = 8)
ax.legend()
fig.savefig(outdir + '/PCA_map.pdf', bbox_inches = 'tight')

nddf = log_df[[c for c in log_df.columns if 'ND' in c]]
t2ddf = log_df[[c for c in log_df.columns if 'T2D' in c]]

def distance_point_line(xp, yp, a, b):
    """point coordinates are (xp,yp); straight line goes like y=a+b*x"""
    x = (xp+b*yp-a*b)/(1+b*b); y = a+b*x
    d2 = (x-xp)**2 + (y-yp)**2
    return np.sqrt(d2)

def coef_variation_df(ndf):
    """returns a dataframe with mu, SD, variance, coefficient of variance, and distance from Poissan behavior for each index"""
    cvdf = pd.DataFrame({'mu': ndf.mean(axis=1), 'SD': ndf.std(axis=1), 'var': ndf.var(axis=1)})
    cvdf['CV'] = cvdf.apply(lambda x: x['SD']/x['mu'], axis = 1)
    cvdf['distance'] = cvdf.apply(lambda x: distance_point_line(np.log10(x['mu']), np.log10(x['CV']), 0., -0.5), axis = 1)
    return cvdf

cv_nd = coef_variation_df(norm_df[nddf.columns])
cv_t2d = coef_variation_df(norm_df[t2ddf.columns])

fig, ax = plt.subplots()
ax.scatter(cv_nd['distance'], cv_t2d['distance'])
ax.scatter(cv_nd.loc['hsa-miR-127-5p','distance'], cv_t2d.loc['hsa-miR-127-5p','distance'])

fig, ax = plt.subplots()
ax.scatter(cv_nd['var'], cv_t2d['var'])

fig, ax = plt.subplots()
ax.scatter(cv_nd['var'], cv_t2d['var'])
ax.scatter(cv_nd.loc['hsa-miR-127-5p','var'], cv_t2d.loc['hsa-miR-127-5p','var'])
ax.set_xscale('log'); ax.set_yscale('log')

cvl_nd = coef_variation_df(log_df[nddf.columns])
cvl_t2d = coef_variation_df(log_df[t2ddf.columns])

fig, ax = plt.subplots()
ax.scatter(cvl_nd['var'], cvl_t2d['var'], label = 'var(log norm expr)')
ax.scatter(cvl_nd.loc['hsa-miR-127-5p','var'], cvl_t2d.loc['hsa-miR-127-5p','var'])
#ax.set_xscale('log'); ax.set_yscale('log')
xra = np.linspace(cvl_nd['var'].min(), cvl_nd['var'].max(), 100)
ax.plot(xra, xra, ls = '--', c = 'r')
ax.set_xlabel('ND'); ax.set_ylabel("T2D")

def bootstrap_variance(v):
    bootstrap_var = np.array([np.random.choice(v, size = len(v)).var() for i in range(100)])
    return pd.Series({'BSvar': bootstrap_var.mean(), 'BSvar_err': bootstrap_var.std()/np.sqrt(len(v))})

var_nndf = pd.DataFrame(log_df[nddf.columns].apply(lambda x: bootstrap_variance(x), axis = 1))
var_t2ddf = pd.DataFrame(log_df[t2ddf.columns].apply(lambda x: bootstrap_variance(x), axis = 1))
var_bs = var_nndf.merge(var_t2ddf, how = 'inner', left_index = True, right_index = True, suffixes = ['_ND','_T2D'])

fig, ax = plt.subplots()
ax.scatter(var_nndf['BSvar'], var_t2ddf['BSvar'], zorder = 2, alpha = 0.7, s = 5)
ax.errorbar(var_nndf['BSvar'], var_t2ddf['BSvar'], xerr = var_nndf['BSvar_err'], yerr = var_t2ddf['BSvar_err'], lw = 0.5, c = 'silver', zorder = 1, fmt = 'o', ms = 0)
#ax.errorbar(var_nndf.loc['hsa-miR-127-5p','BSvar'], var_t2ddf.loc['hsa-miR-127-5p','BSvar'], xerr = var_nndf.loc['hsa-miR-127-5p','BSvar_err'], yerr = var_t2ddf.loc['hsa-miR-127-5p','BSvar_err'], fmt = 'o', c = 'orange')
xra = np.linspace(0,1,100)
ax.plot(xra, xra, ls = '--', c = 'k', lw = 1)
ax.set_xlabel("Var(ND)"); ax.set_ylabel("Var(T2D)")
#ax.plot(xra, (xra+0.05)*1.1, ls = '--', c = 'k', lw = 1)
#ax.plot(xra, (xra-0.05)*0.9, ls = '--', c = 'k', lw = 1)
ax.fill_between(xra, (xra-0.05)*0.9,  (xra+0.05)*1.1, color = 'gray', alpha = 0.25, zorder = 0)

top_t2d_variable = var_bs[(var_bs['BSvar_T2D']-var_bs['BSvar_err_T2D'])>1.1*(var_bs['BSvar_ND']+0.05)]
#ax.scatter(var_nndf.loc[top_t2d_variable.index,'BSvar'], var_t2ddf.loc[top_t2d_variable.index,'BSvar'])
top_nd_variable = var_bs[(var_bs['BSvar_T2D']+var_bs['BSvar_err_T2D'])<0.9*(var_bs['BSvar_ND']-0.05)]
#ax.scatter(var_nndf.loc[top_nd_variable.index,'BSvar'], var_t2ddf.loc[top_nd_variable.index,'BSvar'])
fig.savefig(outdir + '/variance_comparison.pdf', bbox_inches = 'tight')

top_t2d_variable.to_csv(outdir + '/top_t2d_variable.tsv', sep = '\t')
