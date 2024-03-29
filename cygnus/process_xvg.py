import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sns.set_style("ticks")
sns.set_palette("colorblind")


class XvgFile:
    '''For handling xvg files produced by GROMACS

    Attributes
    ----------
    xvg_file : str
        Name of the xvg file supplied
    x_data : np.array
        x_data from the xvg file, generated by process_data()
    y_data : np.array
        y_data from the xvg file, generated by process_data()
    x_label : str
        x_label from the xvg file, generated by process_data()
    y_label : str
        y_label from the xvg file, generated by process_data()

    Methods
    -------
    process_data()
        Reads the xvg file, extracts x_data, y_data and the respective labels
    plot(ax=None, format_plot=True, **kwargs)
        Provides a basic xy plot for the data (not suitable for all cases)
    '''

    def __init__(self, xvg_file):
        '''
        Parameters
        ----------
        xvg_file : str
            The name of the xvg file
        '''

        self.xvg_file = xvg_file

    def process_data(self):
        '''Extracts data

        Reads the xvg file and returns the x_data, y_data, x_label and y_label attributes

        Return self added to allow chaining
        '''

        with open(self.xvg_file, 'r') as file:
            x_data = []
            y_data = []
            for line in file:
                if line.startswith('#'):
                    pass
                elif line.startswith('@'):
                    if 'xaxis' in line:
                        data = line.strip().split('"')
                        self.x_label = data[1]
                    elif 'yaxis' in line:
                        data = line.strip().split('"')
                        self.y_label = data[1]
                else:
                    data = line.strip().split()
                    x_data.append(float(data[0]))
                    y_data.append(float(data[1]))

            self.x_data = np.array(x_data)
            self.y_data = np.array(y_data)

        return self

    def plot(self, ax=None, format_plot=True, **kwargs):
        '''Basic xy plot

        Plots the x and y data on the supplied axes
        Also does some basic formatting, which can be tuned further

        '''

        if ax is None:
            ax = plt.gca()

        ax.plot(self.x_data, self.y_data, **kwargs)

        try:
            ax.set(xlabel=self.x_label,
                   ylabel=self.y_label)
        except AttributeError:
            pass

        sns.despine()

        return ax


class EM(XvgFile):
    '''For plotting energy minimisations'''

    def plot(self, ax=None, format_plot=True, **kwargs):
        if ax is None:
            ax = plt.gca()

        ax.plot(self.x_data, self.y_data, color='black', **kwargs)

        if format_plot:
            ax.set(xlabel='Steps',
                   ylabel='Potential Energy (kJ/mol)',)
            sns.despine()

        return ax


class PullForce(XvgFile):
    '''For plotting pulling simulations

    Assumes kJ/mol/nm as units

    '''

    def plot(self, ax=None, format_plot=True, **kwargs):
        '''Plots the pull force (in black)'''
        if ax is None:
            ax = plt.gca()

        ax.plot(self.x_data, self.y_data, color='black', linewidth=0.8,
                **kwargs)

        if format_plot:
            ax.set(xlabel='Time (ps)',
                   ylabel='Force (kJ/mol/nm)',)
            sns.despine()

        return ax


class PMF(XvgFile):
    '''For potential of mean force plots

    Assumes kcal/mol as units

    '''

    def plot(self, ax=None, format_plot=True, **kwargs):
        '''Plots the PMF and fixes the x label'''
        if ax is None:
            ax = plt.gca()

        ax.plot(self.x_data, self.y_data, **kwargs)

        if format_plot:
            ax.set(xlabel=r'$\xi$ (nm)',
                   ylabel='PMF (kcal/mol)',)
            sns.despine()

        return ax

    def get_dG(self, units='kcal/mol'):
        '''Get the delta G for the PMF'''
        self.dG_value = min(self.y_data) - max(self.y_data)
        print(f'Calculated dG value is {self.dG_value:.2f} {units}')


class Histo(XvgFile):
    '''For gmx wham output

    Special case where the xvg file contains more than 2 columns of data
    Has different attributes from the parent class

    Attributes
    ----------
    bins : list
        Histogram bins
    counts : list
        Histogram counts

    '''

    def process_data(self):
        '''Processes xvg file with more than 2 columns'''

        with open(self.xvg_file, 'r') as file:
            self.bins = []
            self.counts = []
            for line in file:
                if line.startswith(('#', '@')):
                    pass
                else:
                    data = line.strip().split()
                    self.bins.append(float(data[0]))
                    self.counts.append(data[1:])

    def plot(self, ax=None, **kwargs):
        '''Histogram plot'''

        window = [f'window_{i}' for i in np.arange(1, len(self.counts[0])+1)]
        df = pd.DataFrame(self.counts, index=self.bins, columns=window)
        df = df.astype(float)

        if ax is None:
            ax = plt.gca()

        for i, col in enumerate(df.columns, start=1):
            sns.kdeplot(data=df, x=df.index,
                        weights=df[col], label=f'Window {i}', ax=ax,
                        **kwargs)

        ax.set(xlim=(0, 7),
               xlabel=r'$\xi$ (nm)')
        ax.legend(frameon=False)
        sns.despine()

        return ax
