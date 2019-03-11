"""Module for plotting UV spectra from TD-DFT output of Gaussian09
Currently designed for use with Jupyter for inline plotting of data

Initial test version takes mangled txt input from existing extract scripts
Intention is to refactor for compatibility with new python extract API

simple mung available in jupyter_spec.sh to produce current input format
from TD-DFT .log file"""

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
DrawingOptions.bondLineWidth=2.5

from os import path
import numpy as np
import matplotlib.pyplot as plt

class Structure:
    """Takes filename of excited state information from g09 calculation
    Instantiates class with methods to read ES data and build spectrum
    Currently coded for output from getg09ES.awk script prepended with
    smiles info from babel
    Correctly formatted input can be generated as follows:
    babel file.log file.smi
    getg09ES.awk name.log >> file.smi && mv file.smi file.txt"""
    
    def __init__(self, filename, name=None):
        with open(filename) as file:
            self.content = file.readlines()
        self.smiles = self.content.pop(0).split()[0]
        self.logfile = self.content.pop(0)
        if name == None:
            self.name = path.splitext(self.logfile)[0]
        else:
            self.name = self.logfile
    
    def make_structure(self):
        '''Generate image of molecule'''
        self.mol = Chem.MolFromSmiles(self.smiles)
        return self.mol
    
    def read_states(self):
        '''Return list of tuples for singlet state energies and oscillator strength'''
        self.states = [ (float(line.split()[2]), float(line.split()[7])) for line in self.content if 'Singlet' in line]
        return self.states
    
    def gen_spectrum(self, X=np.linspace(1242/200, 1242/900, 1000), width=0.3):
        '''Generate spectrum from ES information in wavelength domain
        Sets and returns (X, Y) array tuple.  Default region is 200nm-900nm'''
        Y = np.zeros_like(X)
        for state in self.states:
            Y += self._transition(X, state[0], width) * state[1]
        self.spectrum = 1242 / X, Y
        return self.spectrum

    def gen_linespec(self):
        '''Generate linespectrum in wavelength domain from excited states
        Sets and returns (X, Y) array tuple'''
        X, Y = zip(*self.states)
        self.linespec = 1242 / np.array(X), np.array(Y)
        return self.linespec

    def gen_plot(self, content='all', plot=None):
        '''Generates pyplot from available data.  By default, instantiates new Figure
        and Axes objects.  Existing Plot, as (fig, ax), can be passed for adding overlays
        ---------
        Plotted content can be specified:
        content="all" - plots all available spectra
        content="line" - plots line spectrum
        content="spec" - plots homogeneously broadened spectrum
        ---------
        Returns (fig, ax) and list of added lines'''
        fig, ax = plot or plt.subplots()
        lines = []
        if content == 'all' or content == 'spec':
            lines.append(ax.plot(self.spectrum[0], self.spectrum[1], label=self.name))
        if content == 'all' or content == 'line':
            lines.append(ax.stem(self.linespec[0], self.linespec[1], markerfmt=' ', basefmt=' ', linefmt='k'))
        return fig, ax, lines
    
    @staticmethod
    def _transition(arr, energy, width):
        return np.exp(-np.square(arr - energy)/(2 * width ** 2))

class Pair:            
    # might easily generalise to n Structure container with this as special instance...
    # Should remove plotting function from Pair and put in underlying structure class
    # Pair just requests plots from contained classes
    """Class design for optical switch
    Takes filenames for UV data from 'open' and 'closed' forms
    
    Method to plot both UV spectra on single axes"""
    def __init__(self, op, cl):
        self.structures = [Structure(op), Structure(cl)]
    
    def spectrum(self):
        '''Builds list of spectra for contained Structures'''
        for struc in self.structures:
            try: struc.states
            except AttributeError: struc.read_states()
        self.spectra = [ s.gen_spectrum() for s in self.structures]
        return self.spectra

    def lines(self):
        '''Builds list of line spectra for contained Structures'''
        for struc in self.structures:
            try: struc.states
            except AttributeError: struc.read_states()
        self.linespectra = [ s.gen_linespec() for s in self.structures]
        return self.linespectra
    
    def plot_spec(self, content='all'):        # refactor for general names and loop self.structures
        '''Plots UV spectra of structure pair'''
        fig, ax = plt.subplots()
        lines = []
        for struc in self.structures:
            fig, ax, l = struc.gen_plot(content=content, plot=(fig, ax))
            lines.append(l)
        ax.legend()
        self.plot = fig, ax, lines
        return self.plot
    
    @staticmethod
    def plot_together(pair1, pair2):
        '''Method to include 2 pairs in single plot'''
        fig, ax= plt.subplots()
        for pair in [pair1, pair2]:
            name = [s.name for s in pair.structures]
            ax.plot(pair.spectra[0][0], pair.spectra[0][1], label=name[0])
            ax.plot(pair.spectra[1][0], pair.spectra[1][1], label=name[1])
        ax.legend()
