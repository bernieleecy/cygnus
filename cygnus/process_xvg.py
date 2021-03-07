import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sns.set_style("ticks")
sns.set_palette("colorblind")


class XvgFile:
    '''
    Basic template
    Includes a basic plotting function
    '''
    def __init__(self, xvg_file):
        self.xvg_file = xvg_file

    def process_data(self):
        with open(self.xvg_file, 'r') as file:
            x_data = []
            y_data = []
            for line in file:
                if line.startswith(('#','@')):
                    pass
                else:
                    data = line.strip().split()
                    x_data.append(float(data[0]))
                    y_data.append(float(data[1]))

            self.x_data = np.array(x_data)
            self.y_data = np.array(y_data)

    def plot(self, ax=None, format_plot=True, **kwargs):
        if ax is None:
            ax = plt.gca()

        ax.plot(self.x_data, self.y_data)

        return ax

class EM(XvgFile):
    '''
    For handling energy minimisation plots
    '''
    def plot(self, ax=None, format_plot=True, **kwargs):
        if ax is None:
            ax = plt.gca()

        ax.plot(self.x_data, self.y_data, color='black')

        if format_plot:
            ax.set(xlabel='Steps',
                   ylabel='Potential Energy (kJ/mol)',)
            sns.despine()

        return ax

class PullForce(XvgFile):
    '''
    For plotting pullf.xvg files from pulling simulations
    Assumes kJ/mol/nm as units
    '''
    def plot(self, ax=None, format_plot=True, **kwargs):
        if ax is None:
            ax = plt.gca()

        ax.plot(self.x_data, self.y_data, color='black', linewidth=0.8)

        if format_plot:
            ax.set(xlabel='Time (ps)',
                   ylabel='Force (kJ/mol/nm)',)
            sns.despine()

        return ax

class PMF(XvgFile):
    '''
    For handling PMF plots
    Assumes kcal/mol as units
    '''
    def plot(self, ax=None, format_plot=True, **kwargs):
        if ax is None:
            ax = plt.gca()

        ax.plot(self.x_data, self.y_data)

        if format_plot:
            ax.set(xlabel=r'$\xi$ (nm)',
                   ylabel='PMF (kcal/mol)',)
            sns.despine()

        return ax

    def get_dG(self, units='kcal/mol'):
        '''
        Get the delta G (binding free energy)
        '''
        self.dG_value = min(self.y_data) - max(self.y_data)
        print(f'Calculated dG value is {self.dG_value:.2f} {units}')

class Histo(XvgFile):
    '''
    For plotting histograms from gmx wham
    '''
    def process_data(self):
        with open(self.xvg_file, 'r') as file:
            self.bins = []
            self.counts = []
            for line in file:
                if line.startswith(('#','@')):
                    pass
                else:
                    data = line.strip().split()
                    self.bins.append(float(data[0]))
                    self.counts.append(data[1:])

    def plot(self, ax=None, **kwargs):
        window = [f'window_{i}' for i in np.arange(1, len(self.counts[0])+1)]
        df = pd.DataFrame(self.counts, index=self.bins, columns=window)
        df = df.astype(float)

        if ax is None:
            ax = plt.gca()

        for i, col in enumerate(df.columns, start=1):
            sns.kdeplot(data=df, x=df.index,
                        weights=df[col], label=f'Window {i}', ax=ax)

        ax.set(xlim=(0,7),
               xlabel=r'$\xi$ (nm)')
        ax.legend(frameon=False)
        sns.despine()