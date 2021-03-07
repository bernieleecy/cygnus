import MDAnalysis as mda
from MDAnalysis.analysis import distances, align, diffusionmap
from MDAnalysis.analysis.base import analysis_class

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import seaborn as sns
from typing import List, Dict, Tuple

sns.set_style("ticks")
sns.set_palette("colorblind")


'''
---------------
 FUNCTIONS
---------------
'''

def traj_dist_calc(u, ref, sel):
    '''
    For running calc_dist on all frames of the supplied trajectory
    Designed to work with a single reference atom
    Returns a (25005,1,X) matrix, where X is the number of atoms in the
    selection
    '''
    dist = _Distances(u.trajectory, u, ref, sel).run()
    return np.squeeze(dist.results)


@analysis_class
def _Distances(u, ref, sel):
    '''
    Calculates distance between all atoms of the reference and selection
    Designed to work with a single reference atom
    With the decorator, it becomes an analysis class
    '''
    dist = distances.distance_array(ref.positions, sel.positions,
                                    box=u.dimensions)
    return dist


def compare_rmsf(rmsf_objs, ax=None):
    '''
    Compare RMSF
    '''
    if ax is None:
        ax = plt.gca()
    fig = ax.get_figure()

    for obj in rmsf_objs:
        ax.plot(obj.df.iloc[:,0], label=obj.label)
        print(obj.label)

    ax.set(xlabel='Residue number',
           ylabel='RMSF (nm)',
           ylim=(0,1.0),
           )

    ax.legend(frameon=False)
    sns.despine()
    fig.tight_layout()


class MDData:
    def __init__(self, gro_file, traj_files):
        self.gro_file = gro_file
        self.traj_files = traj_files
        self.u = mda.Universe(gro_file, traj_files)

    def process_data(self):
        raise NotImplementedError

    def md_plot(self):
        raise NotImplementedError


class PairwiseRMSD(MDData):
    def __init__(self, gro_file, traj_files):
        self.gro_file = gro_file
        self.traj_files = traj_files

    def process_data(self, stride, selection, timestep):
        '''
        Process MD data
        '''
        self.u = mda.Universe(self.gro_file, self.traj_files,
                              in_memory=True, in_memory_step=stride)
        self.traj_len = len(self.u.trajectory)
        self.n_runs = len(self.traj_files)

        # calculate total simulation time and convert from ps to ns
        # use floor division to handle the extra frame in the traj reader
        self.sim_time = int((self.traj_len * stride * timestep) // 1000)

        # do the pairwise RMSD calculation
        align.AlignTraj(self.u, self.u, select=selection, in_memory=True).run()
        self.matrix = diffusionmap.DistanceMatrix(self.u, select=selection,
                                                  verbose=True).run()


    def md_plot(self, ax=None, **kwargs):
        '''
        Plotting pairwise RMSD data
        '''
        if ax is None:
            ax = plt.gca()
        fig = ax.get_figure()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)

        # to mark out divisions between runs
        frames_per_run = self.traj_len // self.n_runs
        time_per_run = self.sim_time // self.n_runs

        # set the ticks to match the number of frames used
        ticks = np.arange(0, self.traj_len, frames_per_run)

        # set the times to match the actual run lengths
        # add 1 to self.sim_time to generate the last time
        times = np.arange(0, self.sim_time + 1, time_per_run)

        im = ax.imshow(self.matrix.dist_matrix, cmap='inferno', origin='lower')
        fig.colorbar(im, cax=cax, orientation='vertical',
                     label=r'RMSD ($\rm \AA$)')

        # lines to mark out divisions between runs
        for i in range(frames_per_run, self.traj_len, frames_per_run):
            ax.axvline(x=i, color='white', linestyle='--')
            ax.axhline(y=i, color='white', linestyle='--')

        ax.set(xlabel='Time (ns)',
               ylabel='Time (ns)',
               xticks=ticks,
               yticks=ticks,
               xticklabels=times,
               yticklabels=times,
               title='Pairwise RMSD',
               **kwargs
               )

        fig.tight_layout()
        return fig, ax


class InterAtomDist(MDData):
    '''
    For calculating distances between sets of atoms
    '''
    def __init__(self, gro_file, traj_files, label):
        self.gro_file = gro_file
        self.traj_files = traj_files
        self.label = label
        self.u = mda.Universe(gro_file, traj_files)

    def process_data(self, dist_dict):
        self.distances = []
        print(dist_dict)

        for key, entry in dist_dict.items():
            if isinstance(entry, list):
                print(key, entry)
                for i in range(len(entry)):
                    dist = self._process_dict_dist(key, entry[i])
                    self.distances.append(dist)
        else:
            print(key, entry)
            dist = self._process_dict_dist(key, entry)
            self.distances.append(dist)

        self.df = pd.DataFrame(self.distances).T
        self.df.columns = [i+1 for i in range(len(self.df.columns))]

    def md_plot(self, ax=None,
                scale: str = 'area',
                **kwargs):
        '''
        For plotting a single violin plot of distances
        '''
        if ax is None:
            ax = plt.gca()
        fig = ax.get_figure()

        sns.violinplot(data=self.df, ax=ax, scale=scale)

        ax.set(xlabel='Interaction',
               ylabel=r'Distance ($\rm \AA$)',
               **kwargs
               )

        sns.despine()
        fig.tight_layout()
        return fig, ax

    def hued_violins(self, other, ax=None,
                     xlabel: str = 'Interaction',
                     ylabel: str = r'Distance ($\rm \AA$)',
                     **kwargs):
        '''
        For making violin plots
        Accepts 2 instances of HBondsData and their respective labels
        '''
        if len(self.df.columns) != len(other.df.columns):
            raise ValueError('Dataframes must have the same number of columns')

        self.df['Category'] = self.label
        other.df['Category'] = other.label

        combined_df = pd.concat([self.df, other.df])
        long_df = combined_df.melt(id_vars='Category')

        if ax is None:
            ax = plt.gca()
        fig = ax.get_figure()

        sns.violinplot(x='variable', y='value', hue='Category',
                       data=long_df, inner=None, scale='area',
                       scale_hue=True, split=True)

        ax.set(xlabel=xlabel,
               ylabel=ylabel,
               **kwargs
               )
        ax.legend(frameon=False)

        sns.despine()
        fig.tight_layout()
        return fig, ax

    def _process_dict_dist(self, key, entry):
        '''
        Takes the universe + the dictionary key and entry as arguments.
        Dictionary key is the reference, entry/entries as the selection
        Designed to work with a SINGLE reference atom

        traj_dist_calc is then run to calculate the interatomic distance
        between  the reference and selection over all frames.
        For an entry with 1 reference atom and 1 selection atom, the distance
        array is returned as a (25005,1,1) matrix ((25000,) after np.squeeze)

        With more than one selection atom, the distance array is returned as a
        (25005,1,x) matrix (where x = number of atoms in the selection, will be
        (25005,x) after np.squeeze)
        This is due to residues like aspartic/glutamic acid, where both
        carboxylate atoms get selected (e.g. name OE1 OE2)

        To average out the ref-OE1 and ref-OE2 distance, the dimensions of the
        distance array are checked with dist.ndim A (25005,2) matrix has
        dist.ndim = 2, so it will get averaged across axis=1 (i.e. average the
        two columns)
        '''
        ref = self.u.select_atoms(key)
        sel = self.u.select_atoms(entry)
        dist = traj_dist_calc(self.u, ref, sel)

        if dist.ndim > 1:
            print(np.mean(dist, axis=1).shape)
            return(np.mean(dist, axis=1))
        else:
            return(dist)


class RMSFData(MDData):
    '''
    GROMACS is preferred for calculating RMSD and RMSF
    Processing of the xvg file is done by default
    '''
    def __init__(self, xvg_file, label):
        '''
        Reads the xvg_file and label, stores the label
        Then processes the xvg_file
        '''
        self.xvg_file = xvg_file
        self.label = label
        with open(self.xvg_file) as infile:
            res_no = []
            val = []

            for line in infile:
                if not (line.startswith("#") or line.startswith("@")):
                    # two spaces between each entry
                    data = line.strip().split('  ')
                    res_no.append(float(data[0].strip()))
                    val.append(float(data[1].strip()))

            self.df = pd.DataFrame(val, index=res_no)
            self.df['roll'] = self.df.rolling(5, center=True).mean()

    def md_plot(self, ax=None):
        '''
        Single RMSF plot
        '''
        if ax is None:
            ax = plt.gca()
        fig = ax.get_figure()

        ax.plot(self.df.iloc[:,0], label='Raw data', alpha=0.1)
        ax.plot(self.df.iloc[:,1], label='Rolling average', color='C0')
        ax.set(xlabel='Residue number',
               ylabel='RMSF (nm)',
               title=f'RMSF ({self.label})',
               )

        ax.legend(frameon=False)
        sns.despine()
        fig.tight_layout()

