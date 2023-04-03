from tqdm import tqdm
from itertools import product
from Bio import SeqIO
import pandas as pd
from g4h import G4hunter
import pandas as pd
import matplotlib.pyplot as plt
import pickle
from collections import defaultdict

df = pd.read_excel('New_Data_for_G4_oligo.xlsx')
#Define the matrices for indivisual transcription factor binding site
#first element is BREu = SSRCGCC, length = 7 bases, size of matrix = 4*7, max score = 3.33
# first row = A, second row = T, third row = G and forth row = C
BREu = [
    [0.0, 0.0, 0.5, 0, 0, 0, 0],
    [0.0, 0.0, 0.0, 0, 0, 0, 0],
    [0.5, 0.5, 0.5, 0, 1, 0, 0],
    [0.5, 0.5, 0.0, 1, 0, 1, 1]]

#second element is ETS = MGGAART, length = 7 bases, size of matrix = 4*7, max score = 6
ETS = [
    [0.5, 0, 0, 1, 1, 0.5, 0],
    [0.0, 0, 0, 0, 0, 0.0, 1],
    [0.0, 1, 1, 0, 0, 0.5, 0],
    [0.5, 0, 0, 0, 0, 0.0, 0]]

#Third element is INR = YYANWYY, length = 7 bases, size of matrix = 4*7, max score = 3.75
INR = [
    [0.0, 0.0, 1, 0.25, 0.5, 0.0, 0.0],
    [0.5, 0.5, 0, 0.25, 0.5, 0.5, 0.5],
    [0.0, 0.0, 0, 0.25, 0.0, 0.0, 0.0],
    [0.5, 0.5, 0, 0.25, 0.0, 0.5, 0.5]]

#forth element is NFY = CCAAT, length = 5 bases, size of matrix = 4*5, max score = 5
NFY = [
    [0, 0, 1, 1, 0],
    [0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0],
    [1, 1, 0, 0, 0]]

#seventh element is SP1_NC = GGGNGGG, length = 7 bases, size of matrix = 4*7, max score = 7
#Gives some possibilities of having any bases at position 4 while still favoring G
#threshold score > 6 : keep all degenerated motif with degenerated = 6.5 and canonical = 7
SP1 = [
    [0, 0, 0, .5, 0, 0, 0],
    [0, 0, 0, .5, 0, 0, 0],
    [1, 1, 1, .5, 1, 1, 1],
    [0, 0, 0, 1., 0, 0, 0]]

#sixth element is TATA = TATAWAAG, length = 8 bases, size of matrix = 4*8, max score = 7.5
#threshold score > 20 : keep all degenerated motif with degenerated = 20 and canonical = 23
TATA = [
    [0, 4, 0, 4, 4, 1, 1, 0],
    [4, 0, 4, 0, 4, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 0]]

def find_motif_in_sequence(seq, motif, threshold):
    """find motif occurences in given sequence
    
    return list of (motif sequence, start, end, score)

    """
    seq = seq.upper()
    motif_size = len(motif[0])
    window_starts = len(seq) - motif_size
    bases_encodings = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    coded_seq = [
        bases_encodings[b] if b in 'ATGC' else None
        for b in seq]
    scores = [
        sum(
            0 if coded_b is None else motif[coded_b][i]
            for i, coded_b in enumerate(coded_seq[start:start+motif_size])
        ) for start in range(window_starts)]
    
    return [
        (seq[i:i+motif_size], i, i+motif_size, score)
        for i, score in enumerate(scores)
        if score >= threshold]

motifs = dict(
    BREu = BREu,
    ETS = ETS,
    INR = INR,
    NFY = NFY,
    SP1 = SP1,
    TATA = TATA)

thresholds = dict(
    BREu = 5.5,
    ETS = 6,
    INR = 3.75,
    NFY = 5,
    SP1 = 6,
    TATA = 20)

hits_per_sequence = [
    [seq,sure,{motif_name: find_motif_in_sequence(
        seq, motif_matrix, thresholds[motif_name])
        for motif_name, motif_matrix in motifs.items()}]
        for seq,sure in zip(df['corrected_seq'],df['SURE sense'])]

with open('/Users/ritika/Work/esnault et al/files-july/mapped_promoter_elements_with_seq_n_sure.pckl', 'bw') as h:
      pickle.dump(hits_per_sequence, h)

with open('/Users/ritika/Work/G4hunter scripts/motif2.pckl', 'rb') as h:
    new_hits_per_sequence = pickle.load(h)

wanted_motifs = ('BREu', 'ETS', 'INR', 'NFY', 'SP1', 'TATA')
all_motif_combinations = list(product(*[(True, False)] * len(wanted_motifs)))
#search for how product in itertools work
combination_results = {
    '-'.join(motif for motif, b in zip(wanted_motifs, combi) if b): [
        (idx, entry) for idx, entry in enumerate(new_hits_per_sequence)
        if (
            all(entry[motif] for motif, b in zip(wanted_motifs, combi) if b)
            and not any(entry[motif] for motif, b in zip(wanted_motifs, combi) if not b)
        )
    ] for combi in all_motif_combinations
}

assert(
    sum(len(item) for item in combination_results.values())
    == len(new_hits_per_sequence))

combination_results2 = {}
for combi in all_motif_combinations:
    # define motifs that we wants
    motifs_that_I_want = []
    motifs_that_I_do_not_want = []
    for motif, b in zip(wanted_motifs, combi):
        if b:
            motifs_that_I_want.append(motif)
        else:
            motifs_that_I_do_not_want.append(motif)
    # define key string
    present_motifs_string = '-'.join(motifs_that_I_want)
    combination_results2[present_motifs_string] = []
    # identify sequences that have all these motifs present
    for idx, entry in enumerate(new_hits_per_sequence):
        all_motifs_that_I_want_are_present = True
        for motif in motifs_that_I_want:
            if len(entry[motif]) == 0:
                all_motifs_that_I_want_are_present = False
                break
        all_motifs_that_I_do_not_want_are_absent = True
        for motif in motifs_that_I_do_not_want:
            if len(entry[motif]) != 0:
                all_motifs_that_I_do_not_want_are_absent = False
                break
        if (
            all_motifs_that_I_want_are_present
            and all_motifs_that_I_do_not_want_are_absent
        ):
            combination_results2[present_motifs_string].append((idx, entry))

assert(
    sum(len(item) for item in combination_results2.values())
    == len(new_hits_per_sequence))

assert(list(combination_results.keys()) == list(combination_results2.keys()))
assert(combination_results == combination_results2)


combination_results = {}
for combi in all_motif_combinations:
    motifs_that_I_want_present = []
    for motif, b in zip(wanted_motifs, combi):
        if b:
            motifs_that_I_want_present.append(motif)
    present_motifs_string = '-'.join(motifs_that_I_want_present)
    combination_results[present_motifs_string] = []
    for idx,entry in enumerate(new_hits_per_sequence):
        listl = []
        for k, v in entry.items():
            if len(v) != 0:
                listl.append(k)
        
        if set(listl) == set(motifs_that_I_want_present):
            if listl != motifs_that_I_want_present:
                print(listl)
                print(motifs_that_I_want_present)
            combination_results[present_motifs_string].append((idx, entry))

assert(
    sum(len(item) for item in combination_results.values())
    == len(new_hits_per_sequence))


final_data = [(k, len(v)) for k, v in combination_results.items()]
final_data.sort(key=lambda x: x[1], reverse=True)

    
wanted_motifs = ('BREu', 'ETS', 'INR', 'NFY', 'SP1', 'TATA')
desired_locations = {
    'BREu': (172, 192), 
    'ETS': (176, 196), 
    'INR': (187, 207),
    'NFY': (135, 155),
    'SP1': (143, 163),
    'TATA':(159, 179)
}

hits = [
    entry for entry in new_hits_per_sequence
    if all(
        any(
            desired_locations[motif][0] >= hit[1] >= desired_locations[motif][1]
            for hit in entry[motif][1]
        ) for motif in wanted_motifs
    )]


'''
list_elements = []
List_BREu = []
List_ETS = []
List_INR = []
List_NFY =[]
List_SP1all = []
List_TATAall = []
for element in all_hits_per_sequence:
    print(element)
    print(type(element))
    for element_key, element_value in element.items():
        if element_key in wanted_motifs():
            print(element_value)
            List_BREu.append(element_)

for i in len(element['INR']):
    List_INR.append(element['INR'][i][1])
    print(List_INR)
    pass

Startsbymotifs = defaultdict(list)
for dictionary in listofresults:
    For motif, data in dict.items():
         Startsbymotifs[motif].extend(data[1])

#dataframe[column][row][number of motifs found][1]
wanted_motifs = ('BREu', 'ETS', 'INR', 'NFY', 'SP1_all', 'TATA_all')
for col in wanted_motifs:
    l = []
    for row in range(len(df)):
        if len(df[col][row])>0:
            for e in range(len(df[col][row])):
                l.append(df[col][row][e][1])
    print(l)
'''
