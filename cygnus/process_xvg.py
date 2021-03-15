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
<<<<<<< HEAD
    Needs to be able to deal with xlabel and ylabel in basic cases
    '''

=======
    '''
>>>>>>> cfdf00559097178646e8a964c3cec77c3c1f7347
    def __init__(self, xvg_file):
        self.xvg_file = xvg_file

    def process_data(self):
        with open(self.xvg_file, 'r') as file:
            x_data = []
            y_data = []
            for line in file:
<<<<<<< HEAD
                if line.startswith('#'):
                    pass
                elif line.startswith('@'):
                    if 'xaxis' in line:
                        data = line.strip().split('"')
                        self.x_label = data[1]
                    elif 'yaxis' in line:
                        data = line.strip().split('"')
                        self.y_label = data[1]
=======
                if line.startswith(('#','@')):
                    pass
>>>>>>> cfdf00559097178646e8a964c3cec77c3c1f7347
                else:
                    data = line.strip().split()
                    x_data.append(float(data[0]))
                    y_data.append(float(data[1]))

            self.x_data = np.array(x_data)
            self.y_data = np.array(y_data)

    def plot(self, ax=None, format_plot=True, **kwargs):
        if ax is None:
            ax = plt.gca()

<<<<<<< HEAD
        ax.plot(self.x_data, self.y_data, **kwargs)
        ax.set(xlabel=self.x_label,
               ylabel=self.y_label)
        sns.despine()

        return ax


=======
        ax.plot(self.x_data, self.y_data)

        return ax

>>>>>>> cfdf00559097178646e8a964c3cec77c3c1f7347
class EM(XvgFile):
    '''
    For handling energy minimisation plots
    '''
<<<<<<< HEAD

=======
>>>>>>> cfdf00559097178646e8a964c3cec77c3c1f7347
    def plot(self, ax=None, format_plot=True, **kwargs):
        if ax is None:
            ax = plt.gca()

<<<<<<< HEAD
        ax.plot(self.x_data, self.y_data, color='black', **kwargs)
=======
        ax.plot(self.x_data, self.y_data, color='black')
>>>>>>> cfdf00559097178646e8a964c3cec77c3c1f7347

        if format_plot:
            ax.set(xlabel='Steps',
                   ylabel='Potential Energy (kJ/mol)',)
            sns.despine()

        return ax

<<<<<<< HEAD

=======
>>>>>>> cfdf00559097178646e8a964c3cec77c3c1f7347
class PullForce(XvgFile):
    '''
    For plotting pullf.xvg files from pulling simulations
    Assumes kJ/mol/nm as units
    '''
<<<<<<< HEAD

=======
>>>>>>> cfdf00559097178646e8a964c3cec77c3c1f7347
    def plot(self, ax=None, format_plot=True, **kwargs):
        if ax is None:
            ax = plt.gca()

<<<<<<< HEAD
        ax.plot(self.x_data, self.y_data, color='black', linewidth=0.8,
                **kwargs)
=======
        ax.plot(self.x_data, self.y_data, color='black', linewidth=0.8)
>>>>>>> cfdf00559097178646e8a964c3cec77c3c1f7347

        if format_plot:
            ax.set(xlabel='Time (ps)',
                   ylabel='Force (kJ/mol/nm)',)
            sns.despine()

        return ax

<<<<<<< HEAD

=======
>>>>>>> cfdf00559097178646e8a964c3cec77c3c1f7347
class PMF(XvgFile):
    '''
    For handling PMF plots
    Assumes kcal/mol as units
    '''
<<<<<<< HEAD

=======
>>>>>>> cfdf00559097178646e8a964c3cec77c3c1f7347
    def plot(self, ax=None, format_plot=True, **kwargs):
        if ax is None:
            ax = plt.gca()

<<<<<<< HEAD
        ax.plot(self.x_data, self.y_data, **kwargs)
=======
        ax.plot(self.x_data, self.y_data)
>>>>>>> cfdf00559097178646e8a964c3cec77c3c1f7347

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

<<<<<<< HEAD

=======
>>>>>>> cfdf00559097178646e8a964c3cec77c3c1f7347
class Histo(XvgFile):
    '''
    For plotting histograms from gmx wham
    '''
<<<<<<< HEAD

=======
>>>>>>> cfdf00559097178646e8a964c3cec77c3c1f7347
    def process_data(self):
        with open(self.xvg_file, 'r') as file:
            self.bins = []
            self.counts = []
            for line in file:
<<<<<<< HEAD
                if line.startswith(('#', '@')):
=======
                if line.startswith(('#','@')):
>>>>>>> cfdf00559097178646e8a964c3cec77c3c1f7347
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
<<<<<<< HEAD
                        weights=df[col], label=f'Window {i}', ax=ax,
                        **kwargs)

        ax.set(xlim=(0, 7),
=======
                        weights=df[col], label=f'Window {i}', ax=ax)

        ax.set(xlim=(0,7),
>>>>>>> cfdf00559097178646e8a964c3cec77c3c1f7347
               xlabel=r'$\xi$ (nm)')
        ax.legend(frameon=False)
        sns.despine()
