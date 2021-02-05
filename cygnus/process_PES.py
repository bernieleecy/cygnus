import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_style("ticks")
sns.set_palette("colorblind")


class OrcaData:
    '''
    Currently only has 1 child class (InterEnergy), but may need to expand in the
    future
    '''
    def __init__(self, energy_file):
        '''
        Energy file (already grepped)
        '''
        self.energy_file = energy_file

    def process_data(self):
        raise NotImplementedError

    def plot(self):
        raise NotImplementedError

class InterEnergy(OrcaData):
    '''
    Interaction energy
    Usually supplied as a text file of reaction energy data
    Currently assumes use of the counterpoise correction
    '''
    def __init__(self, energy_file, theory_level):
        self.energy_file = energy_file
        self.theory_level = theory_level

    def process_data(self):
        with open(self.energy_file, 'r') as file:
            # energies have already been converted to kcal/mol
            # float(data[-2]) is the dimer energy
            # float(data[-1]) is the BSSE correction (subtract it)

            self.dist = []
            self.energies = []

            for line in file:
                data = line.strip().split() # only need the first and last columns
                file_data = data[0].split('_')

                # remove the .out + convert to integer
                self.dist.append(float(file_data[1]))

                energy_kcal = float(data[-2]) - float(data[-1])
                self.energies.append(energy_kcal)

            min_index = self.energies.index(min(self.energies))
            self.min_dist = self.dist[min_index]
            self.min_energy = self.energies[min_index]
            print(f'The mimina is at {self.min_dist}, {self.min_energy}')

    def plot(self, ax=None, format_plot=True,
                  **kwargs):
        '''
        For plotting a 1D PES
        '''
        if ax is None:
            ax = plt.gca()
        fig = ax.get_figure()

        ax.plot(self.dist, self.energies, marker='x', label=self.theory_level)

        if format_plot:
            ax.set(xlabel='Distance (Å)',
                   ylabel='Interaction energy (kcal/mol)',
                   **kwargs)
            ax.legend(frameon=False)
            sns.despine()
            fig.tight_layout()
            print('Formatting specified in plot')

        else:
            print('Formatting not specified in plot')

        return fig, ax


class AmberData:
    '''
    Monomer energies hard coded in
    '''
    def __init__(self, energy_file, force_field='AMBER99sb-ILDN'):
        self.energy_file = energy_file
        self.force_field = force_field

    def process_data(self, monomer=None):
        convert_kcal = 4.184
        tma_energy = 385.581
        monomer_energies = {'indole': 58.634,
                            'benzene': 27.138}
        if monomer.lower() in monomer_energies.keys():
            complex_energy = (monomer_energies[monomer.lower()]
                              + tma_energy) / convert_kcal
        else:
            raise ValueError('Invalid monomer name supplied')

        self.dist = []
        self.energies = []

        with open(self.energy_file, 'r') as file:
            for line in file:
                if line.strip() and '.xvg' in line:
                    data = line.strip().split('_')
                    tma_dist = float(data[1][:4])
                    self.dist.append(tma_dist * 10)
                elif line.strip():
                    data = line.strip().split()
                    int_energy = (float(data[1]) / convert_kcal) \
                                 - complex_energy
                    self.energies.append(int_energy)

            min_index = self.energies.index(min(self.energies))
            self.min_dist = self.dist[min_index]
            self.min_energy = self.energies[min_index]
            print(f'The mimina is at {self.min_dist}, {self.min_energy}')

    def write_data(self, out_file='processed_data.txt'):
        with open(out_file, 'w') as save_file:
            for x, y in zip(self.dist, self.energies):
                print(f'{x:.2f}, {y:.2f}', file=save_file)

    def plot(self, ax=None, format_plot=True,
                **kwargs):
        '''
        For plotting a 1D PES
        '''
        if ax is None:
            ax = plt.gca()
        fig = ax.get_figure()

        ax.plot(self.dist, self.energies, marker='x', label=self.force_field)

        if format_plot:
            ax.set(xlabel='Distance (Å)',
                   ylabel='Interaction energy (kcal/mol)',
                   **kwargs)
            ax.legend(frameon=False)
            sns.despine()
            fig.tight_layout()
            print('Formatting specified in plot')

        else:
            print('Formatting not specified in plot')

        return fig, ax
