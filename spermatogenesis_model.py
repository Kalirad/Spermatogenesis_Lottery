"""
A simple model of stochastic distribution of proteins among sperms during spermatogenesis.
"""

__author__ = 'Ata Kalirad'

__version__ = '1.0'

import os
import uuid
import numpy as np
import pandas as pd

class Cell(object):
    
    def __init__(self, proteome):
        """Initialize an instance of Cell object.
        
        Arguments:
            proteome {dict} -- A dictionary where each keys refer to transcript ID and values refer to counts.
        """

        self.proteome = proteome

class Spermatogenesis(object):
    
    def __init__(self, ref, cutoff, crucial_prot, cv=0.2): 
        """Initialize object.
        
        Arguments:
            ref {dict} -- A dictionary of protein (or transcript) IDs and their counts/
            cutoff {float} -- The percentage of each functional proteins that has 
                              to be inherited from spermatogonia to a sperm to make the sperm viable
            crucial_prot {float} -- The percentage of protein types required for the functionality of a sperm.
        
        Keyword Arguments:
            cv {float} -- The coefficient of variation for cell-to-cell variability (default: {0.2})
        """
        assert type(cutoff) is float
        assert type(crucial_prot) is float
        self.ref = ref
        self.cv = cv
        random_proteome = self.randomize_proteome()
        self.generate_spermatogonia(random_proteome)
        self.cutoff = cutoff 
        self.crucial_prot = crucial_prot

    def randomize_proteome(self):
        """Randomly draw the number of protein based on the given mean and variance.
        
        Returns:
            dict -- randomized dictionary
        """
        random_proteome = {}
        for i in self.ref:
            mean = self.ref[i]
            std = mean * self.cv
            random_proteome[i] = int(np.random.normal(mean, std))
        return random_proteome

    def generate_spermatogonia(self, proteome):
        """Create a spermatogonium with a given proteome (or transcriptome)
        
        """
        self.spermatogonia = Cell(proteome)
        
    def fission(self, mother_prot):
        """Divide a cell into two daughter cells.
        
        Arguments:
            mother_prot {dict} -- The proteome (or transcriptome) of the mother cell.
        
        Returns:
            Cell
        """
        daughter_1_prot = [np.random.binomial(1, 0.5, i).sum() for i in mother_prot.values()]
        daughter_1_transcript = {}
        for i,j in zip(mother_prot, daughter_1_prot):
            daughter_1_transcript[i] = j
        daughter_1 = Cell(daughter_1_transcript)
        daughter_2_prot = list(np.array(mother_prot.values()) - np.array(daughter_1_prot))
        daughter_2_transcript = {}
        for i,j in zip(mother_prot, daughter_2_prot):
            daughter_2_transcript[i] = j
        daughter_2 = Cell(daughter_2_transcript)
        return daughter_1, daughter_2
    
    def do_meiosis(self):
        """Do two rounds of cell division to simulate meiosis
        
        Returns:
            Cell 
        """
        c1, c2 = self.fission(self.spermatogonia.proteome)
        s1, s2 = self.fission(c1.proteome)
        s3, s4 = self.fission(c2.proteome)
        return [s1, s2, s3, s4]
        
    def is_viable(self, sperm):
        """Determine if a sperm is viable.
        
        Arguments:
            sperm {Cell} 
        
        Returns:
            Boolean 
        """
        viable = False
        prot_count = 0
        for i in sperm.proteome:
            if sperm.proteome[i] >= self.cutoff * self.ref[i]:
                prot_count += 1
        prot_prop = prot_count/float(len(sperm.proteome))
        if prot_prop >= self.crucial_prot:
            viable = True
        return viable
    
    def make_sperms(self, t):
        """Make a number of sperms.
        
        Arguments:
            t {int} -- For a given t, 4t sperms will be produced.
        """
        sperms = []
        for i in range(t):
            one_round = self.do_meiosis()
            sperms.extend(one_round)
            random_proteome = self.randomize_proteome()
            self.generate_spermatogonia(random_proteome)
        self.load = sperms
        
    def get_viable_sprems(self):
        """Count viable sperms.
        
        Returns:
            int -- number of viable sperms
        """
        viable = [self.is_viable(i) for i in self.load]
        return sum(viable)

    def run_and_save(self, directory, num=1e6):
        """Make sperms and save the result.
        
        """
        self.make_sperms(num)
        a = self.get_viable_sprems()
        if not os.path.exists(directory):
            os.makedirs(directory)
        with open(directory + "/run_T_" + str(self.cutoff) +"_N_" + str(self.crucial_prot) + "_" + str(uuid.uuid4()), 'w') as file:
            file.write(str(a))

if __name__ == "__main__":
    import doctest
    doctest.testmod()

