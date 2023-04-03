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

def main(argv):
   if not argv:
       sys.stdout.write("Sorry: you must specify at least an argument\n")
       sys.stdout.write("More help avalaible with -h or --help option\n")
       sys.exit(1)

   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:w:s:",["help","ifile=","ofile="])
   except getopt.GetoptError:
      print('\033[1m' +'python G4Hunter.py -i <inputfile> -o <outputrepository> -w <window> -s <score threshold>\n'+'\033[0;0m')
      sys.exit(1)
      
   for opt, arg in opts:
      if opt in ('-h',"--help"):
          print('\033[1m' +'\t ----------------------'+'\033[0;0m')
          print('\033[1m' +'\n\t  Welcome To G4Hunter :'+'\033[0;0m')
          print('\033[1m' +'\t ----------------------'+'\033[0;0m')

          print('\n G4Hunter takes into account G-richness and G-skewness of a given sequence and gives a quadruplex propensity score as output.')
          print('To run G4Hunter use the commande line: \n')
          print('\033[1m' +'python G4Hunter.py -i <inputfile> -o <outputrepository> -w <window> -s <score threshold>\n'+'\033[0;0m')
          sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile= arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
      elif opt in ("-w", "--wind"):
         window = arg
      elif opt in ("-s", "--score"):
         score = arg
             
   return  inputfile, outputfile, int(window), float(score)
   print(outputfile)

#calcule le score de chaque base dans un fichier qui peut contenir plusieur sequences fasta
#le fichier doit comporter la nom de la seq + la sequence en une seul ligne donc pas de \n ou \r 
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

    def BaseScore(self,line):
        #line = sequence from input file
        item, liste=0, []
        #calcule le item de chaque base et la stock dans liste
        while ( item < len(line)):
            #a la fin d une sequence il est possible d avoir des GGG dans se cas
            # on verifie si la secore+1<len(line) car il ya un deuxieme G 
            #et 
            if (item < len(line) and (line[item]=="G" or line[item]=="g")):
                liste.append(1)
                #print liste
                if(item+1< len(line) and (line[item+1]=="G" or line[item+1]=="g")):
                    liste[item]=2
                    liste.append(2)
                    if (item+2< len(line) and (line[item+2]=="G" or line[item+2]=="g")):
                        liste[item+1]=3
                        liste[item]=3
                        liste.append(3)
                        if (item+3< len(line) and (line[item+3]=="G" or line[item+3]=="g")):
                            liste[item]=4
                            liste[item+1]=4
                            liste[item+2]=4
                            liste.append(4)
                            item=item+1
                        item=item+1
                    item=item+1
                item=item+1
                while(item < len(line) and (line[item]=="G" or line[item]=="g")):
                        liste.append(4)
                        item=item+1
        
            #elif (item < len(line) and (line[item]=="T" or line[item]=="A"  or line[item]=="t" or line[item]=="a" or line[item]=="U"or line[item]=="u" or
            #line[item]=="-" or line[item]=="N" or line[item]=="_" or line[item]=="Y"
            #or line[item]=="W" or line[item]=="R" or
            #line[item]=="K" or line[item]=="M"or #line[item]=="S" or line[item]=="B"
            #or line[item]=="V"or line[item]=="D"or line[item]=="H"or line[item]=="N")):
            #liste.append(0)
            #item=item+1
            elif (item < len(line) and line[item]!="G" and line[item]!="g" and line[item]!= "C" and line[item]!="c" ):
                        liste.append(0)
                        item=item+1
                
            elif(item < len(line) and (line[item]=="C" or line[item]=="c")):
                liste.append(-1)
                if(item+1< len(line) and (line[item+1]=="C" or line[item+1]=="c" )):
                    liste[item]=-2
                    liste.append(-2)
                    if (item+2< len(line) and (line[item+2]=="C" or line[item+2]=="c" )):
                        liste[item+1]=-3
                        liste[item]=-3
                        liste.append(-3)
                        if (item+3< len(line) and (line[item+3]=="C" or line[item+3]=="c"  )):
                            liste[item]=-4
                            liste[item+1]=-4
                            liste[item+2]=-4
                            liste.append(-4)
                            item=item+1
                        item=item+1   
                    item=item+1
                item=item+1
                while(item < len(line) and (line[item]=="C" or line[item]=="c")):
                    liste.append(-4)
                    item=item+1
            
            else:
                    item=item+1 #la fin du la ligne ou y a des entrers
        return line, liste
 
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
        return [
            start for start, score in enumerate(scores)
            if (score >= score_threshold) or (score <= -score_threshold)
        ]

    def GetG4(self,line,liste, Window,k):
        LG4=[]
        for i in range(len(liste)) :
            if (liste[i]>= float(Window) or liste[i]<= - float(Window)):
                LG4.append(i)
        return LG4

    def consolidate_G4s(self, seq, scores, g4_window_starts, window_size):
        seq = seq.upper()
        G4_data = []
        idx = 0
        length = len(g4_window_starts)
        while idx < length:
            current_start = pos = g4_window_starts[idx]
            next_idx = idx + 1
            while next_idx < length:
                if g4_window_starts[next_idx] == pos + 1:
                    pos += 1
                    next_idx += 1
                else:
                    break
            # remove leading AT
            while seq[current_start] in ('A', 'T'):
                current_start += 1
            # remove trailing AT
            while seq[pos+window_size-1] in ('A', 'T'):
                pos -= 1
            subseq = seq[current_start:pos+window_size]
            if pos == current_start: # no concatenation
                score = scores[current_start]
            else:
                score = self.compute_score_array(subseq).mean()
            G4_data.append(
                dict(
                    seq_name = seq,
                    start=current_start, #start should be current_start+1
                    end=pos+window_size,
                    score=score,
                    G4=subseq
                
                ))
            idx = next_idx
        
        return G4_data if G4_data else None

    def consolidate_G4_overlaps(
            self,seq_id, seq, scores, g4_window_starts, window_size):

        consolidated_G4_data = dict(
            n_windows=len(g4_window_starts),
            sequence_name = seq_id,
            sequence = seq,
            starts=[],
            ends=[],
            sizes=[],
            strands=[],
            scores=[],
            G4s=[]
        )

        if not g4_window_starts:
            return None

        starts_by_strand = {'+': [], '-': []}
        for start in g4_window_starts:
            if scores[start] >= 0:
                starts_by_strand['+'].append(start)
            else:
                starts_by_strand['-'].append(start)

        def trim_and_score(consolidated_G4_data, start, stop, score, seq, strand, removables):
            # remove leading bases that decrease score
            while seq[start] in removables:
                start += 1
                current_score = None
            # remove trailing bases that decrease score
            while seq[stop-1] in removables:
                stop -= 1
                current_score = None
            subseq = seq[start:stop]
            if score is None:
                score = \
                    self.compute_score_array(subseq).mean()
            if strand == '-':
                subseq = str(Seq(subseq).reverse_complement())

            consolidated_G4_data['starts'].append(start+1)
            consolidated_G4_data['ends'].append(stop)
            consolidated_G4_data['sizes'].append(stop-start)
            consolidated_G4_data['strands'].append(strand)
            consolidated_G4_data['scores'].append(abs(score))
            consolidated_G4_data['G4s'].append(subseq)
            
        for strand in ('+', '-'):

            if not starts_by_strand[strand]:
                # no data for this strand
                continue

            removables = self.removable_bases_by_strand[strand]            
            start = starts_by_strand[strand][0]
            stop = start + window_size
            score = scores[start]

            if len(starts_by_strand[strand]) > 1:
                # more than 1 window to merge
                for next_start in starts_by_strand[strand][1:]:
                    if next_start < stop:
                        stop = next_start + window_size
                        score = None
                    else:
                        # First treat the current data
                        trim_and_score(
                            consolidated_G4_data,
                            start, stop, score,
                            seq, strand, removables)
                        # Second reassign new data to current
                        start = next_start
                        stop = start + window_size
                        score = scores[start] #doubt = score of next start ?
    
            # process very last element
            trim_and_score(
                consolidated_G4_data, start, stop, score, seq, strand, removables)
            
        consolidated_G4_data['max_score'] = max(consolidated_G4_data['scores'])
        consolidated_G4_data['n_g4'] = len(consolidated_G4_data['scores'])
        consolidated_G4_data['total_length'] = sum(consolidated_G4_data['sizes'])

        return consolidated_G4_data


'''
    def WriteSeq(self,line,liste, LISTE, F ):
        if not LISTE:
            return None

        i,k,I=0,0,0
        a=b=LISTE[i]
        MSCORE=[]
    
        if (len(LISTE)>1):
            c=LISTE[i+1]
            while (i< len(LISTE)-2):
                if(c==b+1):
                    k=k+1
                    i=i+1
                else:
                    I=I+1
                    seq=line[a:a+F+k]
                    sequence,liste2=self.BaseScore(seq)
                    MSCORE.append(abs(round(np.mean(liste2),2)))
                    k=0
                    i=i+1
                    a=LISTE[i]
                b=LISTE[i] 
                c=LISTE[i+1] 
            I=I+1
            seq=line[a:a+F+k+1]
            sequence,liste2=self.BaseScore(seq)
            MSCORE.append(abs(round(np.mean(liste2),2)))
            
            
        #dans le cas ou il ya une seul sequence donc pas de chevauchement
        else:
            I=I+1
            seq=line[a:a+F]
            MSCORE.append(abs(liste[a]))
        return MSCORE




    


#calcule le score de chaque base dans un fichier qui peut contenir plusieur sequences fasta
#le fichier doit comporter la nom de la seq + la sequence en une seul ligne donc pas de \n ou \r  
if __name__ == "__main__":
    try:
        inputfile, outputfile , window, score = main(sys.argv[1:])
        fname=inputfile.split("/")[-1]
        name=fname.split(".")

    except ValueError:
        print('\033[1m' +"\n \t Oops! invalide parameters  \n" +'\033[0;0m')
        print("--------------------------------------------------------------------\n")
        sys.exit()
    except UnboundLocalError:
        print('\033[1m' +"\n \t Oops! invalide parameters  \n" +'\033[0;0m')
        print("--------------------------------------------------------------------\n")
        sys.exit()

    OPF= os.listdir(outputfile)
    flag=False
    for dir in OPF:
        DIR="Results_"+str(name[0])
        if dir== DIR:
            print ("true",DIR)
            flag=True
    if flag==True:
        shutil.rmtree(outputfile+"/"+DIR+"/")
        os.makedirs(outputfile+"/"+DIR+"/", mode=0o777)
        print('\033[1m' +"\n \t Re-evaluation of G-quadruplex propensity with G4Hunter " +'\033[0;0m')
        print("\n#####################################")
        print("#    New Results directory Created  #")
        print("#####################################\n")
    else:
        os.makedirs(outputfile+"/"+DIR+"/", mode=0o777)
        print("\n########################################################################")
        print("#                            Results directory Created                 #")
        print("########################################################################\n")
    
    #================================================================
    plot=[]
    fname=inputfile.split("/")[-1]
    filefasta=fname.split(".")
    filein=open(inputfile,"r")
    print("\n Input file:", '\033[1m' + filefasta[0]+'\033[0;0m')
    #repertoire des fichiers de sortie

    Res1file= open (outputfile +"/"+DIR+"/"+filefasta[0]+"-W"+ str(window)+"-S"+str(score)+".txt", "w")
    Res2file= open (outputfile +"/"+DIR+"/"+filefasta[0]+"-Merged.txt", "w")
    #=========================================
    
    startTime = time.time()
    
    
    soft1=G4hunter()
    ScoreListe, DNASeq, NumListe, HeaderListe=G4hunter.GFinder_from_list(filein, window)
    for i in range(len(DNASeq)):
        G4Seq=G4hunter.get_G4s_for_sequence(DNASeq[i],Res1file, ScoreListe[i], float(score), int(window), HeaderListe[i], len(NumListe[i]))
        if (len(G4Seq)>0):
            MSCORE=soft1.WriteSeq(DNASeq[i],Res2file,ScoreListe[i], G4Seq, HeaderListe[i], int(window), len(NumListe[i]))
    #plot.append(MSCORE)
    """
    malist, alllist=[], []
    #print ScoreListe[0]
    for jj in range (len(ScoreListe[0])):
        cc, mean=0, 0
        for kk in range(len(ScoreListe)-1):
            #print kk, jj, ScoreListe[kk][jj]
            cc+=ScoreListe[kk][jj]
        mean=cc/len(ScoreListe)
        alllist.append(mean)
        if abs(mean) >=score :
            malist.append(mean)
        else:
            malist.append(0)
    #soft1.plot2(ScoreListe[0], outputfile +"/Results/")
    soft1.plot2(malist, outputfile +"/"+DIR+"/", "sc")
    soft1.plot2(alllist, outputfile +"/"+DIR+"/", "all")
    """
    filein.close()
    fin=time.time()

    print("\n Results files and Score Figure are created in:   ")#,fin-startTime, "secondes""
    print('\033[1m' + outputfile,"/",DIR,"/","\n "+'\033[0;0m]')


    Res1file.close()
    Res2file.close()
'''