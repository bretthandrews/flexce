from __future__ import print_function, division, absolute_import
import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
from matplotlib.patheffects import withStroke
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
from general import pickle_read

class PlotSim:
    def __init__(self, suite_in):
        '''zones: list of zones to be analyzed'''
        #self.zone_names = [item.name for item in zn]
        self.suite = suite_in
        self.n_boxes = len(self.suite)
        self.setup()

    def setup(self):
        ''' Read in data from
        Reddy et al. (2003, 2006) [X/Fe]
        Ramirez et al. (2007) NLTE [O/Fe]
        APOGEE red clump [alpha/Fe]'''
        self.stem_pcaa_ce = os.path.expanduser('~') + '/pcaa_chemevol/'
        stem_data = self.stem_pcaa_ce + 'data/'
        stem_reddy = stem_data + 'reddy/'
        stem_ramirez = stem_data + 'ramirez13/'
        stem_apogee_redclump = stem_data + 'apogee_redclump/'
        self.reddy = pickle_read(stem_reddy + 'reddy.pck')
        self.reddy03 = pickle_read(stem_reddy + 'reddy03.pck')
        self.ramirez13 = pickle_read(stem_ramirez + 'ramirez13.pck')
        self.apogee_redclump = pickle_read(stem_apogee_redclump + 
                                           'apogee_redclump_v402.pck')

    def plot_xfe(self,
                 element,
                 runs='all',
                 time_pts=False,
                 line=True,
                 pts=False,
                 sampling_factor=1e8,
                 range=[0.17, 0.15, 0.75, 0.75],
                 xylim=[-3, 0.75, -0.2, 0.8],
                 reddypts=True,
                 apogeercpts=False,
                 ptsize=10,
                 pop=('tn', 'tntk', 'tk', 'tkh', 'h'),
                 datacolor=False,
                 time_pts_offset=[0.02, 0.02],
                 ofe_tick=[0.35, 0.45],
                 figlab=None,
                 legon=True,
                 legtitle=None,
                 leghlen=0.5,
                 legloc=3,
                 leglabsp=0.01,
                 legfont=20,
                 legtitlefont=20,
                 xmajortick=0.5,
                 xminortick=0.5,
                 ymajortick=0.2,
                 yminortick=0.1):
        '''Plot [X/Fe]--[Fe/H] where element = X.  time_pts=True will add
        points to the line to indicate the speed of chemical evolution.
        pts=True will plot as points instead of a line.  '''
        plt.rcParams['xtick.major.pad'] = 10
        plt.rcParams['ytick.major.pad'] = 8
        fig = plt.figure()
        ax = fig.add_axes(range)
        ax.set_axisbelow(True)
        ax.grid(color='gray', which='minor')
        ax.plot([-10, 10], [0, 0], c='k', ls='-', lw=1)
        ax.plot([0, 0], [-10, 10], c='k', ls='-', lw=1)
        if runs == 'all':
            runs = [item for item in self.suite.iterkeys()]
        t = [100, 500, 1000, 2000, 4000]
        t2 = ['100 Myr', '500 Myr', '1 Gyr', '2 Gyr', '4 Gyr']
        m = ['o', 's', '^', 'd', 'v']
        s = [60, 70, 80, 80, 80]
        myeffect = withStroke(foreground='w', linewidth=4)
        if reddypts:
            pdata, labdata = self.plot_reddy_ramirez(element, ax, pop,
                                                     datacolor, s=ptsize)
        if apogeercpts:
            self.plot_apogee_redclump(ax)
        p = []
        lab = []
        tlab = True
        xfe_max = -100.
        xfe_min = 100.
        #for i, (key, sim) in enumerate(self.suite.iteritems()):
        #    if key in runs:
        for i, key in enumerate(runs):
            sim = self.suite[key]
            ind = np.where(sim['ab'].elements == element)[0][0]
            ls = '-'
            if line:
                if 'ls' in sim:
                    ls = sim['ls']
                p.append(ax.plot(sim['ab'].feh, sim['ab'].xfe_all[ind],
                                 c=sim['c'], ls=ls, lw=5, zorder=9-i)[0])
            if pts:
                xs = []
                ys = []
                for tstep in range(1, sim['box'].n_steps):
                    ndots = np.ones(np.around(sim['box'].survivors[tstep] /
                                    sampling_factor).astype(int))
                    xs.append(ndots * sim['ab'].feh[tstep-1])
                    ys.append(ndots * sim['ab'].xfe_all[ind][tstep-1])
                np.random.seed(12345)
                xs = np.random.normal(np.concatenate(xs), 0.05)
                ys = np.random.normal(np.concatenate(ys), 0.02)
                kwargs = dict(facecolor=sim['c'], edgecolor='None', zorder=9-i)
                ax.scatter(xs, ys, s=3, **kwargs)
                p.append(ax.scatter([100], [100], s=40, **kwargs))
            lab.append(sim['name'])
            ax.plot(ofe_tick, np.ones(2) * sim['ab'].xfe_all[ind][-1],
                    c=sim['c'], ls=ls, lw=3, zorder=9-i)
            feht = interp1d(sim['box'].t[1:], sim['ab'].feh)(t)
            xfet = interp1d(sim['box'].t[1:], sim['ab'].xfe_all[ind])(t)
            if np.max(sim['ab'].xfe_all[ind]) > xfe_max:
                xfe_max = np.max(sim['ab'].xfe_all[ind])
            if np.min(sim['ab'].xfe_all[ind]) < xfe_min:
                xfe_min = np.min(sim['ab'].xfe_all[ind])
            if time_pts:
                for j, (x, y) in enumerate(zip(feht, xfet)):
                    ax.scatter(x, y, facecolor='w', edgecolor=sim['c'],
                               marker=m[j], lw=1.5, s=s[j], zorder=10)
                    if (key == 'sim0') or ('sim0' not in runs):
                        if tlab:
                            ax.text(x+time_pts_offset[0],
                                    y+time_pts_offset[1], '%s' % t2[j],
                                    color=sim['c'], fontsize=18,
                                    path_effects=[myeffect], zorder=10)
                if (key == 'sim0') or ('sim0' not in runs):
                    tlab = False
        if figlab is not None:
            ax.text(figlab[0], figlab[1], figlab[2], fontsize=28)
        ax.xaxis.set_major_locator(MultipleLocator(xmajortick))
        ax.xaxis.set_minor_locator(MultipleLocator(xminortick))
        ax.yaxis.set_major_locator(MultipleLocator(ymajortick))
        ax.yaxis.set_minor_locator(MultipleLocator(yminortick))
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(28)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(28)
        ax.set_xlabel('[Fe/H]', fontsize=28)
        ax.set_ylabel('[%s/Fe]' % element, fontsize=28)
        ax.axis(xylim)
        largs = dict(scatterpoints=1, handlelength=leghlen, handletextpad=0.5,
                     labelspacing=leglabsp, borderpad=0.5, borderaxespad=0.5)
        if legtitle is not None:
            largs['title'] = legtitle
        if legon:
            leg = ax.legend(p, lab, loc=legloc, **largs)
            plt.setp(leg.get_texts(), fontsize=legfont)
            plt.setp(leg.get_title(), fontsize=legtitlefont)
        if datacolor:
            if 'title' in largs:
                del largs['title']
            leg2 = ax.legend(pdata, labdata, loc=1, **largs)
            plt.setp(leg2.get_texts(), fontsize=12)
            plt.gca().add_artist(leg)


    def plot_xfe_mdf(self,
                     element,
                     runs='all',
                     ratio=None,
                     time_pts=False,
                     tlab=True,
                     line=True,
                     pts=False,
                     sampling_factor=1e8,
                     range=[0.15, 0.15, 0.63, 0.63],
                     mdf_range=[0.15, 0.8, 0.63, 0.15],
                     xfe_range=[0.8, 0.15, 0.15, 0.63],
                     xylim=[-3, 0.75, -0.2, 0.8],
                     mdf_xylim=[-2, 0.5001, 0, 10],
                     xfe_xylim=[0, 20, -0.1, 0.5],
                     mdf_args={'min': -3., 'max': 1., 'delta': 0.005},
                     xfe_args={'min': -0.1, 'max': 0.5, 'delta': 0.001},
                     mdf_distr_multiplelocator=5,
                     xfe_distr_multiplelocator=10,
                     composite_mdf=False,
                     composite_xfe=False,
                     reddypts=True,
                     apogeercpts=False,
                     ptsize=10,
                     pop=('tn', 'tntk', 'tk', 'tkh', 'h'),
                     datacolor=False,
                     time_pts_offset=[0.02, 0.02],
                     ofe_tick=[0.35, 0.45],
                     figlab=None,
                     legon=True,
                     legtitle=None,
                     leghlen=0.5,
                     legloc=3,
                     leglabsp=0.01,
                     legfont=20,
                     legtitlefont=20,
                     legborderpad=0.5,
                     legframeon=True,
                     legzorder=10,
                     xmajortick=0.5,
                     xminortick=0.5,
                     ymajortick=0.2,
                     yminortick=0.1,
                     ylabelpad=5):
        '''Plot [X/Fe]--[Fe/H] where element = X.  time_pts=True will add
        points to the line to indicate the speed of chemical evolution.
        pts=True will plot as points instead of a line.  '''        
        plt.rcParams['xtick.major.pad'] = 10
        plt.rcParams['ytick.major.pad'] = 8
        fig = plt.figure()
        ax = fig.add_axes(range)
        ax.set_axisbelow(True)
        ax.grid(color='gray', which='minor')
        ax.plot([-10, 10], [0, 0], c='k', ls='-', lw=1)
        ax.plot([0, 0], [-10, 10], c='k', ls='-', lw=1)
        if runs == 'all':
            runs = [item for item in self.suite.iterkeys()]
        frac_stars = np.ones(len(runs))
        if ratio is None:
            ratio = np.ones(len(runs))
        else:
            ind_ratio1 = np.where(np.array(ratio) == 1)[0]
            n_stars_ratio1 = np.sum(self.suite[runs[ind_ratio1]]['box'].survivors/1e6)
            for i, key in enumerate(runs):
                sim = self.suite[key]
                frac_stars[i] = (ratio[i] /
                                 (np.sum(sim['box'].survivors / 1e6) / n_stars_ratio1))
        t = [100, 500, 1000, 2000, 4000]
        t2 = ['100 Myr', '500 Myr', '1 Gyr', '2 Gyr', '4 Gyr']
        m = ['o', 's', '^', 'd', 'v']
        s = [60, 70, 80, 80, 80]
        myeffect = withStroke(foreground='w', linewidth=4)
        if reddypts:
            pdata, labdata = self.plot_reddy_ramirez(element, ax, pop,
                                                     datacolor, s=ptsize)
        if apogeercpts:
                self.plot_apogee_redclump(ax)
        p = []
        lab = []
        xfe_max = -100.
        xfe_min = 100.
        #for i, (key, sim) in enumerate(self.suite.iteritems()):
        #    if key in runs:
        for i, key in enumerate(runs):
            sim = self.suite[key]
            ind = np.where(sim['ab'].elements == element)[0][0]
            ls = '-'
            if line:
                if 'ls' in sim:
                    ls = sim['ls']
                p.append(ax.plot(sim['ab'].feh, sim['ab'].xfe_all[ind],
                                 c=sim['c'], ls=ls, lw=5, zorder=9-i)[0])
            if pts:
                xs = []
                ys = []
                for tstep in range(1, sim['box'].n_steps):
                    ndots = np.ones(np.around(sim['box'].survivors[tstep] / 
                                    sampling_factor * 
                                    frac_stars[i]).astype(int))
                    xs.append(ndots * sim['ab'].feh[tstep-1])
                    ys.append(ndots * sim['ab'].xfe_all[ind][tstep-1])
                np.random.seed(12345)
                xs = np.random.normal(np.concatenate(xs), 0.05)
                ys = np.random.normal(np.concatenate(ys), 0.02)
                kwargs = dict(facecolor=sim['c'], edgecolor='None', zorder=9-i)
                ax.scatter(xs, ys, s=3, **kwargs)
                p.append(ax.scatter([100], [100], s=40, **kwargs))
            lab.append(sim['name'])
            ax.plot(ofe_tick, np.ones(2) * sim['ab'].xfe_all[ind][-1],
                    c=sim['c'], ls=ls, lw=3, zorder=9-i)
            feht = interp1d(sim['box'].t[1:], sim['ab'].feh)(t)
            xfet = interp1d(sim['box'].t[1:], sim['ab'].xfe_all[ind])(t)
            if np.max(sim['ab'].xfe_all[ind]) > xfe_max:
                xfe_max = np.max(sim['ab'].xfe_all[ind])
            if np.min(sim['ab'].xfe_all[ind]) < xfe_min:
                xfe_min = np.min(sim['ab'].xfe_all[ind])
            if time_pts:
                for j, (x, y) in enumerate(zip(feht, xfet)):
                    ax.scatter(x, y, facecolor='w', edgecolor=sim['c'],
                               marker=m[j], lw=1.5, s=s[j], zorder=10)
                    if (key == 'sim0') or ('sim0' not in runs):
                        if tlab:
                            ax.text(x+time_pts_offset[0],
                                    y+time_pts_offset[1], '%s' % t2[j],
                                    color=sim['c'], fontsize=18,
                                    path_effects=[myeffect], zorder=10)
                if (key == 'sim0') or ('sim0' not in runs):
                    tlab = False
        if figlab is not None:
            ax.text(figlab[0], figlab[1], figlab[2], fontsize=28)
        ax.axis(xylim)
        ax.xaxis.set_major_locator(MultipleLocator(xmajortick))
        ax.xaxis.set_minor_locator(MultipleLocator(xminortick))
        ax.yaxis.set_major_locator(MultipleLocator(ymajortick))
        ax.yaxis.set_minor_locator(MultipleLocator(yminortick))
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(28)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(28)
        ax.set_xlabel('[Fe/H]', fontsize=28)
        ax.set_ylabel('[%s/Fe]' % element, fontsize=28, labelpad=ylabelpad)
        largs = dict(scatterpoints=1, handlelength=leghlen, handletextpad=0.5,
                     labelspacing=leglabsp, borderpad=legborderpad,
                     borderaxespad=0.5, frameon=legframeon)
        if legtitle is not None:
            largs['title'] = legtitle
        if legon:
            leg = ax.legend(p, lab, loc=legloc, **largs)
            plt.setp(leg.get_texts(), fontsize=legfont)
            plt.setp(leg.get_title(), fontsize=legtitlefont)
            leg.set_zorder(legzorder)
        if datacolor:
            if 'title' in largs:
                del largs['title']
            leg2 = ax.legend(pdata, labdata, loc=1, **largs)
            plt.setp(leg2.get_texts(), fontsize=12)
            plt.gca().add_artist(leg)
        # MDF
        ax = fig.add_axes(mdf_range)
        ax.set_axisbelow(True)
        ax.grid(color='gray')
        if composite_mdf:
            comp_stellar_pop = []
            comp_wt = []
        for i, key in enumerate(runs):
            sim = self.suite[key]
            kwargs = dict(c=sim['c'], lw=2, zorder=10-i)
            if key in runs:
                ind = np.where(sim['box'].survivors[1:] > 0.)[0]
                stellar_pop = sim['ab'].feh[ind]
                xs = np.arange(mdf_args['min'], mdf_args['max'], mdf_args['delta'])
                wt = (sim['box'].survivors[1:][ind] * frac_stars[i] /
                      1e6).astype(int)
                ys = self.weighted_kde(stellar_pop, wt, xs)
                if 'ls' in sim:
                    kwargs['ls'] = sim['ls']
                if 'alpha' in sim:
                    kwargs['alpha'] = sim['alpha']
                ax.plot(xs, ys, **kwargs)
                if composite_mdf:
                    comp_stellar_pop.append(stellar_pop)
                    comp_wt.append(wt)
        if composite_mdf:
            comp_stellar_pop = np.array(comp_stellar_pop)
            comp_wt = np.array(comp_wt)
            comp_stellar_pop = comp_stellar_pop.flatten()
            comp_wt = comp_wt.flatten()
            xs = np.arange(mdf_args['min'], mdf_args['max'], mdf_args['delta'])
            ys = self.weighted_kde(comp_stellar_pop, comp_wt, xs)
            if 'ls' in sim:
                kwargs['ls'] = sim['ls']
            kwargs['c'] = 'k'
            kwargs['zorder'] = 10
            ax.plot(xs, ys, **kwargs)
        ax.yaxis.set_major_locator(MultipleLocator(mdf_distr_multiplelocator))
        ax.axis(mdf_xylim)
        ax.set_ylabel(r'$P$ ([Fe/H])', fontsize=20)
        plt.setp(ax.get_xticklabels(), visible=False)
        # [X/Fe] distribution
        ax = fig.add_axes(xfe_range)
        ax.set_axisbelow(True)
        ax.grid(color='gray')
        if composite_xfe:
            comp_stellar_pop = []
            comp_wt = []
        for i, key in enumerate(runs):
            sim = self.suite[key]
            kwargs = dict(c=sim['c'], lw=2, zorder=10-i)
            if key in runs:
                ind = np.where(sim['box'].survivors[1:] > 0.)[0]
                wt = (sim['box'].survivors[1:][ind] * frac_stars[i] /
                      1e6).astype(int)
                xs = np.arange(xfe_args['min'], xfe_args['max'], xfe_args['delta'])
                ind_el = np.where(sim['ab'].elements_out == element)[0][0]
                stellar_pop = sim['ab'].xfe[ind_el, ind]
                ys = self.weighted_kde(stellar_pop, wt, xs)
                if 'ls' in sim:
                    kwargs['ls'] = sim['ls']
                if 'alpha' in sim:
                    kwargs['alpha'] = sim['alpha']
                ax.plot(ys, xs, **kwargs)
                if composite_xfe:
                    comp_stellar_pop.append(stellar_pop)
                    comp_wt.append(wt)
        if composite_xfe:
            comp_stellar_pop = np.array(comp_stellar_pop)
            comp_wt = np.array(comp_wt)
            comp_stellar_pop = comp_stellar_pop.flatten()
            comp_wt = comp_wt.flatten()
            xs = np.arange(xfe_args['min'], xfe_args['max'], xfe_args['delta'])
            ys = self.weighted_kde(comp_stellar_pop, comp_wt, xs)
            if 'ls' in sim:
                kwargs['ls'] = sim['ls']
            kwargs['c'] = 'k'
            kwargs['zorder'] = 10
            ax.plot(ys, xs, **kwargs)
        ax.xaxis.set_major_locator(MultipleLocator(xfe_distr_multiplelocator))
        ax.axis(xfe_xylim)
        ax.set_xlabel(r'$P$ ([%s/Fe])' % element, fontsize=20)
        plt.setp(ax.get_yticklabels(), visible=False)




    def plot_apogee_redclump(self, ax, contours=False, **kwargs):
        '''Plot the APOGEE red clump data in [alpha/Fe]--[Fe/H] space.'''
        feh = self.apogee_redclump['feh']
        afe = self.apogee_redclump['alphafe']
        ind_rc = self.apogee_redclump['ind_rc']
        ind_snr150 = self.apogee_redclump['ind_snr150']
        ax.scatter(feh[ind_snr150], afe[ind_snr150], s=3, c='k', zorder=1)



    def plot_reddy_ramirez(self, element, ax, pop, datacolor, s=10):
        '''Plot the data from Reddy et al. (2003, 2006) or Ramirez et
        al. (2013) (NLTE [O/Fe]) in [X/Fe]--[Fe/H] space.'''
        if element == 'O':
            d = self.ramirez13
            y = d['ofe']
        elif element in self.reddy['sym']:
            ind_r = np.where(self.reddy['sym'] == element)[0][0]
            d = self.reddy
            y = d['xfe'][ind_r]
        elif element in self.reddy03['sym']:
            ind_r03 = np.where(self.reddy03['sym'] == element)[0][0]
            d = self.reddy03
            y = d['xfe'][ind_r03]
        x = d['feh']
        colors = dict(tn='b', tntk='c', tk='lime', tkh='orange', h='r')
        label = dict(tn='thin', tntk='thin/thick', tk='thick',
                     tkh='thick/halo', h='halo')
        p = []
        lab = []
        for item in pop:
            try:
                ind = d['ind_' + item]
                kwargs = dict(facecolor='gray', edgecolor='none', s=s,
                              zorder=0)
                if datacolor:
                    kwargs['facecolor'] = colors[item]
                    lab.append(label[item])
                p.append(ax.scatter(x[ind], y[ind], **kwargs))
            except KeyError:
                print('%s not available' % item)
        arg = [p, lab]
        if not datacolor:
            arg = [[p[0]], ['all']]
        return arg


    def plot_reddy_ramirez_fex(self, element, ax, pop, datacolor):
        '''Plot the data from Reddy et al. (2003, 2006) or Ramirez et
        al. (2013) (NLTE [O/Fe]) in [Fe/X]--[X/H] space.'''
        if element == 'O':
            d = self.ramirez13
            y = d['ofe'] * -1.
        elif element in self.reddy['sym']:
            ind_r = np.where(self.reddy['sym'] == element)[0][0]
            d = self.reddy
            y = d['xfe'][ind_r] *1.
        elif element in self.reddy03['sym']:
            ind_r03 = np.where(self.reddy03['sym'] == element)[0][0]
            d = self.reddy03
            y = d['xfe'][ind_r03] * 1.
        x = y * -1. + d['feh']
        #y *= -1.
        colors = dict(tn='b', tntk='c', tk='lime', tkh='orange', h='r')
        label = dict(tn='thin', tntk='thin/thick', tk='thick',
                     tkh='thick/halo', h='halo')
        p = []
        lab = []
        for item in pop:
            try:
                ind = d['ind_' + item]
                kwargs = dict(facecolor='gray', edgecolor='none', s=10,
                              zorder=0)
                if datacolor:
                    kwargs['facecolor'] = colors[item]
                    lab.append(label[item])
                p.append(ax.scatter(x[ind], y[ind], **kwargs))
            except KeyError:
                print('%s not available' % item)
        arg = [p, lab]
        if not datacolor:
            arg = [[p[0]], ['all']]
        return arg


    def plot_feh_time(self, runs='all', pts=False,
                      range=[0.15, 0.15, 0.75, 0.75]):        
        '''Plot [Fe/H] as a function of time.  pts=True will plot as points
        instead of a line.  '''
        fig = plt.figure()
        ax = fig.add_axes(range)
        if runs == 'all':
            runs = [item for item in self.suite.iterkeys()]
        p = []
        lab = []
        for i, (key, sim) in enumerate(self.suite.iteritems()):
            if key in runs:
                if not pts:
                    p.append(ax.plot(sim['box'].t[1:], sim['ab'].feh,
                                     c=sim['c'], lw=2, zorder=9-i)[0])
                else:
                    p.append(ax.scatter(sim['box'].t[1:], sim['ab'].feh,
                                        facecolor=sim['c'], edgecolor='None',
                                        zorder=9-i))
                lab.append(sim['name'])
        ax.plot([-5000, 50000], [0, 0], 'k:')
        ax.set_xlabel('Time [Myr]')
        ax.set_ylabel('[Fe/H]')
        ax.axis([-200, 12000, -3., 0.5])
        leg = ax.legend(p, lab, scatterpoints=1, handlelength=1.,
                        handletextpad=0.05, labelspacing=0.01, borderpad=0.2,
                        borderaxespad=0.2, loc=4)
        plt.setp(leg.get_texts(), fontsize=12)


    def plot_xfe_time(self, element, runs='all', pts=False,
                      range=[0.15, 0.15, 0.75, 0.75]):        
        '''Plot [X/Fe] as a function of time.  pts=True will plot as points
        instead of a line.  '''
        fig = plt.figure()
        ax = fig.add_axes(range)
        if runs == 'all':
            runs = [item for item in self.suite.iterkeys()]
        p = []
        lab = []
        for i, (key, sim) in enumerate(self.suite.iteritems()):
            ind = np.where(sim['ab'].elements == element)[0][0]
            if key in runs:
                if not pts:
                    p.append(ax.plot(sim['box'].t[1:], sim['ab'].xfe_all[ind],
                                     c=sim['c'], lw=2, zorder=9-i)[0])
                else:
                    p.append(ax.scatter(sim['box'].t[1:],
                                        sim['ab'].xfe_all[ind],
                                        facecolor=sim['c'], edgecolor='None',
                                        zorder=9-i))
                lab.append(sim['name'])
        ax.plot([-5000, 50000], [0, 0], 'k:')
        ax.set_xlabel('Time [Myr]')
        ax.set_ylabel('[%s/Fe]' % element)
        rng = ax.axis()
        ax.axis([-200, 12000, rng[2], rng[3]])
        leg = ax.legend(p, lab, scatterpoints=1, handlelength=1.,
                        handletextpad=0.05, labelspacing=0.01, borderpad=0.2,
                        borderaxespad=0.2, loc=1)
        plt.setp(leg.get_texts(), fontsize=12)


    def plot_mdf_ofe(self,
                     runs='all',
                     feh_args={'min': -3., 'max': 1., 'delta': 0.005},
                     ofe_args={'min': -0.1, 'max': 0.5, 'delta': 0.001},
                     mdf_range=[-1, 0.20001, 0, 20],
                     ofe_range=[-0.1, 0.50001, 0., 45],
                     legtitle=None,
                     mdf_xtextloc=-0.4,
                     ofe_xtextloc1=-0.09,
                     ofe_xtextloc2=0.25,
                     figsize=(8, 6)):
        '''Plot MDF (metallicity distribution function; i.e., histogram of
        [Fe/H]). Also plot [O/Fe] distribution for low metallicity and high
        metallicity stars.  Note: np.histogram fails if weights are too large
        '''
        if runs == 'all':
            runs = [item for item in self.suite.iterkeys()]
        fig = plt.figure(figsize=figsize)
        fig.subplots_adjust(left=0.125, right=0.95, bottom=0.125, top=0.95,
                            hspace=0.3)
        xlabels = ['[Fe/H]', '[O/Fe]']
        axis_ranges = [mdf_range, ofe_range]
        for j in range(2):
            ax = fig.add_subplot(2, 1, j+1)
            ax.set_axisbelow(True)
            ax.grid(color='gray')
            p = []
            lab = []
            for i, key in enumerate(runs):
                sim = self.suite[key]
                kwargs = dict(c=sim['c'], lw=2, zorder=10-i)
                if key in runs:
                    if j == 0:
                        ind = np.where(sim['box'].survivors[1:] > 0.)[0]
                        stellar_pop = sim['ab'].feh[ind]
                        xs = np.arange(feh_args['min'], feh_args['max'], feh_args['delta'])
                        wt = (sim['box'].survivors[1:][ind]/1e6).astype(int)
                        ys = self.weighted_kde(stellar_pop, wt, xs)
                        p.append(ax.plot(xs, ys, **kwargs)[0])
                        lab.append(sim['name'])
                        fr = '%0.1e' % (ys[np.where(xs<=-1)[0]].sum()/ys.sum())
                        fr2 = 'e-'.join(fr.split('e-0'))
                        yt = axis_ranges[j][3] * (0.7 - 0.1 * i)
                        ax.text(mdf_xtextloc, yt, fr2, color=sim['c'], fontsize=16)
                    elif j == 1:
                        xs = np.arange(ofe_args['min'], ofe_args['max'], ofe_args['delta'])
                        for k in range(2):
                            #ind = np.where(sim['box'].survivors[1:] > 0.)[0]
                            if k == 0:
                                indt = np.where(sim['ab'].feh <= -0.5)[0]
                                kwargs['ls'] = '--'
                                kwargs['dashes'] = (30, 2, 1, 2)
                            elif k == 1:
                                indt = np.where(sim['ab'].feh > -0.5)[0]
                                kwargs['ls'] = '-'
                                kwargs['dashes'] = (None, None)
                            ind2 = np.intersect1d(ind, indt)
                            stellar_pop = sim['ab'].xfe[2, ind2]
                            wt = (sim['box'].survivors[1:][ind2]/1e6).astype(int)
                            ys = self.weighted_kde(stellar_pop, wt, xs)
                            p.append(ax.plot(xs, ys, **kwargs)[0])
            if j == 0:
                ax.text(mdf_xtextloc, axis_ranges[j][3] * 0.8,
                    'fraction of [Fe/H]$<-$1', fontsize=16)
                largs = dict(scatterpoints=1, handlelength=1.,
                             handletextpad=0.5, labelspacing=0.01,
                             borderpad=0.25, borderaxespad=0.8, loc=2)
                if legtitle is not None:
                    largs['title'] = legtitle
                l1 = ax.legend(p, lab, **largs)
                plt.gca().add_artist(l1)
                plt.setp(l1.get_texts(), fontsize=16)
                plt.setp(l1.get_title(), fontsize=16)
            elif j == 1:
                yt = axis_ranges[j][3] * 0.75
                ax.text(ofe_xtextloc1, yt, r'[Fe/H]$\geq -0.5$')
                ax.text(ofe_xtextloc2, yt, r'[Fe/H]$<-$0.5')
            ax.set_xlabel(xlabels[j], fontsize=20)
            ylab = r'$P$ (%s)' % xlabels[j]
            #if j == 0:
            #    ylab = 'log(' + ylab + ')'
            ax.set_ylabel(ylab, fontsize=20)
            ax.axis(axis_ranges[j])
            if j == 1:
                ax.yaxis.set_major_locator(MultipleLocator(10))



    def weighted_kde(self, data, weights, xs):
        data_wt = np.array([data[i] for i in range(len(data))
                            for w in range(weights[i])])
        density = gaussian_kde(data_wt)
        ys = density(xs)
        return ys


    def plot_histogram(self):
        '''Generalized plotting function to create histograms for an MDF or an
        [O/Fe] distribution.'''



    def plot_mdf(self, runs='all', feh_min=-5., feh_max=1., dfeh=0.1,
                 step=True, range=[0.15, 0.15, 0.75, 0.75], show_gcs=False):
        '''Plot MDF (metallicity distribution function; i.e., histogram of
        [Fe/H]).'''
        fig = plt.figure()
        ax = fig.add_axes(range)
        if runs == 'all':
            runs = [item for item in self.suite.iterkeys()]
        p = []
        lab = []
        for i, (key, sim) in enumerate(self.suite.iteritems()):
            if key in runs:
                bins = np.arange(feh_min, feh_max, dfeh)
                # np.histogram fails if weights are too large
                N, bins0 = np.histogram(sim['ab'].feh, bins=bins, normed=True,
                                        weights=sim['box'].survivors[1:]/1000.)
                kwargs = dict(c=sim['c'], lw=2, zorder=10-i)
                if step:
                    kwargs['drawstyle'] = 'steps-mid'
                p.append(ax.plot(bins[1:], N / N.sum(), **kwargs)[0])
                lab.append(sim['name'])
        if show_gcs:
            gcs = np.genfromtxt(
                self.stem_pcaa_ce + 'data/casagrande11/catalog.dat',
                dtype=[('id', '|S4'),('HIP', int),('Name', '|S13'),
                       ('m_name', '|S4'),('plx', float),('e_plx', float),
                       ('logg', float),('Teff', float),('e_Teff', float),
                       ('Fbol', float),('feh', float),('meh', float),
                       ('afe', float)],
                delimiter=(4,7,13,4,7,7,6,6,7,13,6,6,6))
            gcs['feh'][np.where(gcs['feh'] < -5)] = np.nan
            gcs['meh'][np.where(gcs['meh'] < -5)] = np.nan
            gcs['afe'][np.where(gcs['afe'] < -5)] = np.nan
            bins = np.arange(feh_min, feh_max, dfeh)
            N1, bins1 = np.histogram(gcs['feh'], bins=bins, normed=True)
            kwargs = dict(c='k', lw=4)
            if step:
                kwargs['drawstyle'] = 'steps-mid'
            p.append(ax.plot(bins[1:], N1 / N1.sum(), **kwargs)[0])
            lab.append('GCS')
        ax.set_xlabel('[Fe/H]')
        ax.set_ylabel(r'$N$')
        ax.axis([-4., 0.75, 0., 0.35])
        leg = ax.legend(p, lab, scatterpoints=1, handlelength=1.,
                        handletextpad=0.5, labelspacing=0.01, borderpad=0.5,
                        borderaxespad=0.5, loc=2)
        plt.setp(leg.get_texts(), fontsize=16)


    def plot_xfe_distr(self, element='O', runs='all', xfe_min=-0.25,
                       xfe_max= 1., dxfe=0.05, step=True,
                       range=[0.15, 0.15, 0.75, 0.75]):
        '''Plot [X/Fe] ditribution function (i.e., histogram of [X/Fe]).'''
        fig = plt.figure()
        ax = fig.add_axes(range)
        if runs == 'all':
            runs = [item for item in self.suite.iterkeys()]
        p = []
        lab = []
        for i, (key, sim) in enumerate(self.suite.iteritems()):
            ind = np.where(sim['ab'].elements == element)[0][0]
            if key in runs:
                bins = np.arange(xfe_min, xfe_max, dxfe)
                N, bins0 = np.histogram(sim['ab'].xfe_all[ind], bins=bins,
                                        normed=True,
                                        weights=sim['box'].survivors[1:]/100.)
                kwargs = dict(c=sim['c'], lw=2, zorder=10-i)
                if step:
                    kwargs['drawstyle'] = 'steps-mid'
                p.append(ax.plot(bins[1:], N / N.sum(), **kwargs)[0])
                lab.append(sim['name'])
        ax.set_xlabel('[%s/Fe]' % element)
        ax.set_ylabel(r'N/N$_\mathrm{tot}$')
        ax.axis([-0.25, 1., 0., 0.5])
        leg = ax.legend(p, lab, scatterpoints=1, handlelength=1.,
                        handletextpad=0.05, labelspacing=0.01, borderpad=0.2,
                        borderaxespad=0.2, loc=1)
        plt.setp(leg.get_texts(), fontsize=12)



    def plot_mstar_time(self, runs='all', pts=False,
                       range=[0.15, 0.15, 0.75, 0.75]):        
        '''Plot stellar mass as a function of time.  pts=True will plot as
        points instead of a line.  '''
        fig = plt.figure()
        ax = fig.add_axes(range)
        if runs == 'all':
            runs = [item for item in self.suite.iterkeys()]
        p = []
        lab = []
        for i, (key, sim) in enumerate(self.suite.iteritems()):
            if key in runs:
                if not pts:
                    p.append(ax.plot(sim['box'].t[1:],
                                     np.log10(np.cumsum(
                        np.sum(sim['box'].mstar_left[1:], axis=1))),
                                     c=sim['c'], lw=2, zorder=9-i)[0])
                else:
                    p.append(ax.scatter(sim['box'].t[1:],
                                        np.log10(np.cumsum(
                        np.sum(sim['box'].mstar_left[1:], axis=1))),
                                        facecolor=sim['c'], edgecolor='None',
                                        zorder=9-i))
                lab.append(sim['name'])
        ax.set_xlabel('Time [Myr]')
        ax.set_ylabel(r'log(M$_\star$ [M$_\odot$])')
        leg = ax.legend(p, lab, scatterpoints=1, handlelength=1.,
                        handletextpad=0.05, labelspacing=0.01, borderpad=0.2,
                        borderaxespad=0.2, loc=4)
        plt.setp(leg.get_texts(), fontsize=12)


    def plot_mgas_time(self, runs='all', pts=False,
                       range=[0.15, 0.15, 0.75, 0.75]):        
        '''Plot gas mass as a function of time.  pts=True will plot as points
        instead of a line.  '''
        fig = plt.figure()
        ax = fig.add_axes(range)
        if runs == 'all':
            runs = [item for item in self.suite.iterkeys()]
        p = []
        lab = []
        for i, (key, sim) in enumerate(self.suite.iteritems()):
            if key in runs:
                if not pts:
                    p.append(ax.plot(sim['box'].t[1:],
                             np.log10(np.sum(sim['box'].mgas_iso[1:], axis=1)),
                                     c=sim['c'], lw=2, zorder=9-i)[0])
                else:
                    p.append(ax.scatter(sim['box'].t[1:],
                                        np.log10(np.cumsum(
                        np.sum(sim['box'].mgas_iso[1:], axis=1))),
                                        facecolor=sim['c'], edgecolor='None',
                                        zorder=9-i))
                lab.append(sim['name'])
        ax.set_xlabel('Time [Myr]')
        ax.set_ylabel(r'log(M$_\mathrm{gas}$ [M$_\odot$])')
        leg = ax.legend(p, lab, scatterpoints=1, handlelength=1.,
                        handletextpad=0.05, labelspacing=0.01, borderpad=0.2,
                        borderaxespad=0.2, loc=1)
        plt.setp(leg.get_texts(), fontsize=12)


    def plot_mass_time(self, runs='all', range=[0.15, 0.15, 0.75, 0.75]):
        '''Plot stellar and gas mass as a function of time.'''        
        fig = plt.figure()
        ax = fig.add_axes(range)
        if runs == 'all':
            runs = [item for item in self.suite.iterkeys()]
        pstar = []
        pgas = []
        pwarm = []
        lab = []
        for i, (key, sim) in enumerate(self.suite.iteritems()):
            if key in runs:
                kwargs = dict(c=sim['c'], lw=2, zorder=9-i)
                pstar.append(ax.plot(sim['box'].t[1:], np.log10(np.cumsum(
                    np.sum(sim['box'].mstar_left[1:], axis=1))), **kwargs)[0])
                pgas.append(ax.plot(sim['box'].t[1:], np.log10(np.sum(
                    sim['box'].mgas_iso[1:], axis=1)), ls='--', **kwargs)[0])
                pwarm.append(ax.plot(sim['box'].t[1:], np.log10(np.sum(
                    sim['box'].mwarmgas_iso[1:], axis=1)), ls=':',
                                     **kwargs)[0])
                lab.append(sim['name'])
        ax.set_xlabel('Time [Myr]')
        ax.set_ylabel(r'log(Mass [M$_\odot$])')
        leg = ax.legend(pstar, lab, scatterpoints=1, handlelength=1.,
                        handletextpad=0.05, labelspacing=0.01, borderpad=0.2,
                        borderaxespad=0.2, loc=4)
        plt.setp(leg.get_texts(), fontsize=12)
        leg2 = ax.legend([pstar[0], pgas[0], pwarm[0]],
                         ['stars', 'cold gas', 'warm gas'],
                         scatterpoints=1, handlelength=1.5, handletextpad=0.05,
                         labelspacing=0.01, borderpad=0.1, borderaxespad=0.2,
                         loc=3)
        plt.setp(leg2.get_texts(), fontsize=12)
        plt.gca().add_artist(leg)


    def plot_gas_flow_rate(self, runs='all', log=True,
                           range=[0.15, 0.15, 0.75, 0.75]):
        '''Plot SFR, statistical SFR (before Poisson sampling), inflow rate,
        and outflow rate as a function of time.'''
        fig = plt.figure()
        ax = fig.add_axes(range)
        if runs == 'all':
            runs = [item for item in self.suite.iterkeys()]
        psfr = []
        pstat = []
        pin = []
        pout = []
        pcool = []
        pej = []
        lab = []
        for i, (key, sim) in enumerate(self.suite.iteritems()):
            if key in runs:
                dt = sim['box'].dt * 1e6
                sfr = np.sum(sim['box'].mstar[1:], axis=1) / dt
                sfr_stat = sim['box'].sfr[1:]
                inflow = np.sum(sim['box'].inflow[1:], axis=1) / dt
                outflow = np.sum(sim['box'].outflow[1:], axis=1) / dt
                gas_cooling = np.sum(sim['box'].gas_cooling[1:], axis=1) / dt
                stellar_ej = ((np.sum(sim['box'].snii[1:], axis=1) + 
                              np.sum(sim['box'].snia[1:], axis=1) +
                              np.sum(sim['box'].agb[1:], axis=1)) / dt)
                val = dict(sfr=sfr, sfr_stat=sfr_stat, inflow=inflow,
                           outflow=outflow, gas_cooling=gas_cooling,
                           stellar_ej=stellar_ej)
                if log:
                    for k, v in val.iteritems():
                        val[k] = np.log10(v)
                    ax.set_ylabel(r'log($\dot{\mathrm{M}}$' +
                                  ' [M$_\odot$ yr$^{-1}$])')
                else:
                    ax.set_ylabel(r'$\dot{\mathrm{M}}$ [M$_\odot$ yr$^{-1}$]')
                kwargs = dict(c=sim['c'], lw=2, zorder=9-i)
                t = sim['box'].t[1:]
                psfr.append(ax.plot(t, val['sfr'], **kwargs)[0])
                pstat.append(ax.plot(t, val['sfr_stat'], ls='-.', **kwargs)[0])
                pin.append(ax.plot(t, val['inflow'], ls='--', **kwargs)[0])
                pout.append(ax.plot(t, val['outflow'], ls = ':', **kwargs)[0])
                pcool.append(ax.plot(t, val['gas_cooling'], ls = '--',
                                     dashes=(15,5), **kwargs)[0])
                pej.append(ax.plot(t, val['stellar_ej'], ls = '--',
                                   dashes=(25,5), **kwargs)[0])
                lab.append(sim['name'])
        ax.set_xlabel('Time [Myr]')
        rtmp = list(ax.axis())
        rtmp[1] *= 1.5
        ax.axis(rtmp)
        leg = ax.legend(psfr, lab, scatterpoints=1, handlelength=1.,
                        handletextpad=0.05, labelspacing=0.01, borderpad=0.2,
                        borderaxespad=0.2, loc=1)
        plt.setp(leg.get_texts(), fontsize=12)
        leg2 = ax.legend([psfr[0], pstat[0], pin[0], pout[0], pcool[0],
                          pej[0]], ['SFR', 'statistical SFR', 'inflow rate',
                                    'outflow rate', 'gas cooling rate',
                                    'stellar gas return'], scatterpoints=1,
                         handlelength=1.5, handletextpad=0.05,
                         labelspacing=0.01, borderpad=0.1, borderaxespad=0.2,
                         loc=4)
        plt.setp(leg2.get_texts(), fontsize=12)
        plt.gca().add_artist(leg)


    def plot_snia_time(self, runs='all', pts=False,
                       range=[0.15, 0.15, 0.75, 0.75]):        
        '''Plot number of SNIa as a function of time.  pts=True will plot as
        points instead of a line.  '''        
        fig = plt.figure()
        ax = fig.add_axes(range)
        if runs == 'all':
            runs = [item for item in self.suite.iterkeys()]
        p = []
        lab = []
        for i, (key, sim) in enumerate(self.suite.iteritems()):
            if key in runs:
                dt = sim['box'].dt
                if not pts:
                    p.append(ax.plot(sim['box'].t[1:], sim['box'].NIa[1:] / dt,
                                     c=sim['c'], lw=2, zorder=9-i)[0])
                else:
                    p.append(ax.scatter(
                        sim['box'].t[1:], sim['box'].NIa[1:] / dt,
                        facecolor=sim['c'], edgecolor='None', zorder=9-i)[0])
                lab.append(sim['name'])
        ax.set_xlabel('Time [Myr]')
        ax.set_ylabel(r'$N_\mathrm{SNIa}$ / Myr')
        leg = ax.legend(p, lab, scatterpoints=1, handlelength=1.,
                        handletextpad=0.05, labelspacing=0.01, borderpad=0.2,
                        borderaxespad=0.2, loc=1)
        plt.setp(leg.get_texts(), fontsize=12)



    def plot_fex(self, element, runs='all', time_pts=False, pts=False,
                 range=[0.17, 0.15, 0.75, 0.75], xylim=[-3, 0.75, -0.8, 0.2],
                 pop=('tn', 'tntk', 'tk', 'tkh', 'h'), datacolor=True,
                 figlab=None, legtitle=None, legloc=2, leglabsp=0.01):
        '''Plot [Fe/X]--[X/H] where element = X.  time_pts=True will add
        points to the line to indicate the speed of chemical evolution.
        pts=True will plot as points instead of a line.  '''        
        plt.rcParams['xtick.major.pad'] = 10
        plt.rcParams['ytick.major.pad'] = 8
        fig = plt.figure()
        ax = fig.add_axes(range)
        if runs == 'all':
            runs = [item for item in self.suite.iterkeys()]
        t = [100, 500, 1000, 2000, 4000]
        t2 = ['100 Myr', '500 Myr', '1 Gyr', '2 Gyr', '4 Gyr']
        m = ['o', 's', '^', 'd', 'v']
        s = [60, 70, 80, 80, 80]
        myeffect = withStroke(foreground='w', linewidth=4)
        pdata, labdata = self.plot_reddy_ramirez_fex(element, ax, pop,
                                                     datacolor)
        p = []
        lab = []
        tlab = True
        fex_max = -100.
        fex_min = 100.
        #for i, (key, sim) in enumerate(self.suite.iteritems()):
        #    if key in runs:
        for i, key in enumerate(runs):
            sim = self.suite[key]
            ind = np.where(sim['ab'].elements == element)[0][0]
            if not pts:
                p.append(ax.plot(sim['ab'].xh_all[ind],
                                 sim['ab'].xfe_all[ind] * -1.,
                                 c=sim['c'], lw=5, zorder=9-i)[0])
            else:
                p.append(ax.scatter(sim['ab'].xh_all[ind],
                                    sim['ab'].xfe_all[ind] * -1.,
                                    facecolor=sim['c'], edgecolor='None',
                                    zorder=9-i))
            lab.append(sim['name'])
            xht = interp1d(sim['box'].t[1:], sim['ab'].xh_all[ind])(t)
            fext = interp1d(sim['box'].t[1:], sim['ab'].xfe_all[ind] * -1.)(t)
            if np.max(sim['ab'].xfe_all[ind] * -1.) > fex_max:
                fex_max = np.max(sim['ab'].xfe_all[ind] * -1.)
            if np.min(sim['ab'].xfe_all[ind] * -1.) < fex_min:
                fex_min = np.min(sim['ab'].xfe_all[ind] * -1.)
            if time_pts:
                for j, (x, y) in enumerate(zip(xht, fext)):
                    ax.scatter(x, y, facecolor='w', edgecolor=sim['c'],
                               marker=m[j], lw=1.5, s=s[j], zorder=10)
                    if (key == 'sim0') or ('sim0' not in runs):
                        if tlab:
                            ax.text(x+0.02, y+0.02, '%s' % t2[j],
                                    color=sim['c'], fontsize=18,
                                    path_effects=[myeffect], zorder=10)
                if (key == 'sim0') or ('sim0' not in runs):
                    tlab = False
        ax.plot([-10, 10], [0, 0], 'k:')
        ax.plot([0, 0], [-10, 10], 'k:')
        if figlab is not None:
            ax.text(figlab[0], figlab[1], figlab[2], fontsize=28)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(28)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(28)
        ax.set_xlabel('[%s/H]' % element, fontsize=28)
        ax.set_ylabel('[Fe/%s]' % element, fontsize=28)
        ax.axis(xylim)
        largs = dict(scatterpoints=1, handlelength=0.5, handletextpad=0.5,
                     labelspacing=leglabsp, borderpad=0.5, borderaxespad=0.5)
        if legtitle is not None:
            largs['title'] = legtitle
        leg = ax.legend(p, lab, loc=legloc, **largs)
        plt.setp(leg.get_texts(), fontsize=20)
        if datacolor:
            if 'title' in largs:
                del largs['title']
            leg2 = ax.legend(pdata, labdata, loc=1, **largs)
            plt.setp(leg2.get_texts(), fontsize=12)
        plt.gca().add_artist(leg)


    def find_confidence_interval(x, pdf, confidence_level):
        return pdf[pdf > x].sum() - confidence_level

    def plot_rc_alpha(ax, snr150=False):
        nbins_x = 30
        nbins_y = 40
        ind = indrc
        if snr150:
            ind = ind_snr150
        H, xedges, yedges = np.histogram2d(metals[ind], alphafe[ind],
                                           bins=(nbins_x, nbins_y), normed=True)
        x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1, nbins_x))
        y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y, 1))
        pdf = (H * np.multiply(x_bin_sizes, y_bin_sizes).T)
        sig05 = optimize.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.38))
        sig1 = optimize.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
        sig2 = optimize.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
        sig3 = optimize.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
        levels = [sig1, sig05, 1.]
        X = 0.5 * (xedges[1:] + xedges[:-1])
        Y = 0.5 * (yedges[1:] + yedges[:-1])
        Z = pdf.T
        ax.contourf(X, Y, Z, levels=levels, origin='lower',
                    colors=('darkgray', 'gray'), zorder=9)
        ax.contour(X, Y, Z, levels=[levels[0], levels[1]], origin='lower',
                   colors='k', zorder=10)
        for i in range(nbins_x):
            for j in range(nbins_y):
                if Z[j, i] <= sig1 * 1.2:
                    ind_tmp0 = np.where(
                        (metals >= xedges[i]) & (metals < xedges[i+1]) &
                        (alphafe >= yedges[j]) & (alphafe < yedges[j+1]))[0]
                    ind_tmp = np.intersect1d(ind, ind_tmp0)
                    if len(ind_tmp > 0):
                        ax.scatter(metals[ind_tmp], alphafe[ind_tmp], c='k', s=5)

