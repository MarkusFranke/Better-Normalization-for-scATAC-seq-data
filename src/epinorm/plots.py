import matplotlib.pyplot as plt
import numpy as np
import anndata as ad
from sklearn.metrics import RocCurveDisplay


def plot_count_distribution(adata: ad.AnnData, title: str):
    fig = plt.figure()
    ax = fig.add_subplot()
    fig.suptitle(title)
    ax.set_xlabel('fragment')
    ax.set_ylabel('count')

    data = adata.X.toarray().flatten()
    bins = np.arange(0, data.max() + 1.5) - 0.5
    ax.hist(data, log=True, bins=bins, edgecolor='black')


def plot_count_mean_var(adata: ad.AnnData, title: str):
    means = np.apply_along_axis(np.mean, 1, adata.X.toarray())
    variances = np.apply_along_axis(np.var, 1, adata.X.toarray())
    x_min = np.min(means)
    x_max = np.max(means)
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.grid()
    fig.suptitle(title)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('mean')
    ax.set_ylabel('variance')

    ax.scatter(means, variances, color='lightblue', edgecolors='white')
    ax.plot([x_min, x_max], [x_min, x_max], color='black', label='var = mean')
    ax.plot([x_min, x_max], [2 * x_min, 2 * x_max], color='black', linestyle='dotted', label='var = mean*2')
    ax.plot([x_min, x_max], [x_min / 2, x_max / 2], color='black', linestyle='dashed', label='var = mean/2')
    ax.legend()


def plot_roc_curve(fpr, tpr, auc, key, ax):
    display = RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=auc, estimator_name=key)
    display.plot(ax=ax)
