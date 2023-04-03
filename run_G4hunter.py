from G4hunter import G4hunter
import pandas as pd
import matplotlib.pyplot as plt

inputfile = "../Shira's Data/15753_oligos_measurements_and_sequence.tab"
window_size = 25
# outputfile
score_threshold = 1.5



df = pd.read_csv(inputfile, sep='\t')
seq_ids = df['Oligo_name'].to_list()
seqs = df['Oligo_sequence'].to_list()

hunter = G4hunter()
#create a list to store the dictionaries
all_g4s = []

#get the results
for seq_id, seq in zip(seq_ids, seqs):
    score_array = hunter.compute_score_array(seq)
    scores = hunter.compute_score_windows(score_array, window_size)
    G4_windows_starts = hunter.get_G4_window_starts(seq, scores, score_threshold, window_size)
    all_g4s.append(hunter.consolidate_G4_overlaps(seq_id, seq, scores, G4_windows_starts, window_size))

df_g4 = pd.DataFrame([g4 for g4 in all_g4s if g4 is not None])
df = pd.merge(df, df_g4, left_on='Oligo_name', right_on='sequence_name', how='outer')
df.drop(['sequence_name', 'sequence'], axis=1, inplace=True)
df.to_csv(
    "../Shira's Data/15753_oligos_measurements_and_sequence_with_g4s.csv",
    header=True)
df1 = pd.read_csv("../Shira's Data/15753_oligos_measurements_and_sequence_with_g4s.csv")
