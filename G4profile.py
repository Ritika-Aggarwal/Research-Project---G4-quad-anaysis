import sys
print(sys.version)
print(sys.executable)
import os, sys, getopt
import time
import shutil
import numpy as np
import bottleneck as bn
from matplotlib import pyplot
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import pandas as pd
class G4hunter(object):

    scoring_dict = dict(
        G=(0, 1, 2, 3, 4),
        C=(0, -1, -2, -3, -4),
        A=(0, 0),
        T=(0, 0)
    )
    removable_bases_by_strand = {
        '+': [k for k, v in scoring_dict.items() if v[1] <= 0],
        '-': [k for k, v in scoring_dict.items() if v[1] >= 0],
    }
    # k = key, v = value (in dictionary)

    def __init__(self, custom_scoring_dict=None):
        if custom_scoring_dict is not None:
            self.scoring_dict = custom_scoring_dict

    def set_scoring_dict(self, custom_scoring_dict):
        self.scoring_dict = custom_scoring_dict


    def GFinder_from_list(self, seqs, window_size):
        """Calculate scores from list of sequences"""
        scores = []
        for i in range(len(seqs)):
            score_list = self.compute_score_array(seqs[i])
            scores.append(self.compute_score_windows(score_list, window_size))

        return scores

    def compute_score_array(self, seq):
        idx = 0 
        scores = []
        seq = seq.upper()
        #calculate the idx of each base and the stock in list
        while idx < len(seq):
            current_base = seq[idx]
            # go to next base if current is ambiguous
            if current_base not in self.scoring_dict:
                idx += 1
                continue
            # find stretch of bases
            next_idx = idx + 1
            while next_idx < len(seq) and seq[next_idx] == current_base:
                next_idx += 1
            stretch_len = next_idx - idx
            if stretch_len >= len(self.scoring_dict[current_base]):
                base_score = self.scoring_dict[current_base][-1]
            else:
                base_score = self.scoring_dict[current_base][stretch_len]
            scores += [base_score] * stretch_len
            idx += stretch_len

        return np.array(scores)
 
    def compute_score_windows(self, scores, window_size):
        return bn.move_mean(scores, window=window_size)[window_size-1:]

    def CalScore(self,liste, k):
        #k = window size, liste = score list at base level
        Score_Liste=[]
        #calcule de la moynne des scores pour toutes les sequences - les derniers k bases
        for i in range (len (liste)-(k-1)):
            #print (len(liste)), i, k
            j,Sum=0,0
            while (j<k):
                #print j, i
                Sum=Sum+liste[i]
                j=j+1
                i=i+1
            Mean=Sum/float(k)
            Score_Liste.append(Mean) 
        return Score_Liste
        #score_list = list of scores for window size 

    def plot2(self,liste, repert, i):
        # make a little extra space between the subplots
        plt.subplots_adjust(wspace=1.0)
        dt = 1
        t = np.arange(0, len(liste), dt)
        figure= plt.figure()
        plt.plot(t, liste, 'b-')
        plt.xlim(0,len(liste))
        #figure.suptitle('Score of window nucleotides', fontsize=16)
        plt.xlabel('Position (ntS)')
        plt.ylabel('Score')
        plt.grid(True)
        figure.savefig(repert+'Score_plot_'+i+'.pdf', dpi=figure.dpi)           

    def get_G4_window_starts(self, seq, scores, score_threshold, window_size):
        return [start for start, score in enumerate(scores)]

    def get_G4_window_scores(self, seq, scores, score_threshold, window_size):
        return [score for score in scores]
            