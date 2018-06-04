from __future__ import print_function, division, absolute_import

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


def get_colors(cfg):
    """Get colors from config file or set them with seaborn color palette.

    Args:
        cfg (dict): config settings.

    Returns:
        list: colors
    """
    try:
        colors = cfg['Plot']['colors']
        if not isinstance(colors, list):
            colors = [colors]
    except KeyError:
        colors = sns.color_palette('bright')
    return colors


def get_leg_args(cfg):
    """Get legend args from config file.

    Args:
        cfg (dict): config settings.

    Returns:
        dict: Keyword args to pass to legend.
    """
    leg_args = dict(scatterpoints=1, handletextpad=0.05, labelspacing=0.01,
                    borderpad=0.2, borderaxespad=0.5, loc=3)
    try:
        args = cfg['Legend']
        for k, v in args.items():
            try:
                v = float(v)
            except ValueError:
                continue
            leg_args[k] = v
    except KeyError:
        pass
    return leg_args


def get_path_collections(fig):
    """Get PathCollections (i.e., data points) from Axes instance.

    Args:
        fig: Axes instance.

    Returns:
        list: PathCollections.
    """
    children = fig.ax_joint.get_children()
    PathCollection = mpl.collections.PathCollection
    p = [item for item in children if isinstance(item, PathCollection)]
    return p


def joint_overplot(x, y, df, fig, color='r', marg_kws=None):
    """Overplot additional data on existing JointGrid instance.

    Args:
        x (str):
        y (str):
        df (DataFrame):
        fig: seaborn JointGrid instance.
        color (str): Color.
        marg_kws (dict): Keyword arguments to pass to plot_marginals().

    Returns:
        fig: seaborn JointGrid instance.
    """
    if marg_kws is None:
        marg_kws = dict(norm_hist=True,
                        hist_kws=dict(weights=df.Survivors.values))
    fig.x = df[x]
    fig.y = df[y]
    fig.plot_joint(plt.scatter, c=color)
    fig.plot_marginals(sns.distplot, color=color, kde=False, axlabel=False,
                       **marg_kws)
    return fig
