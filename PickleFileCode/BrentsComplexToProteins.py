import pickle
import pandas as pd
from collections import defaultdict

df = pd.read_csv('/Users/brent/Desktop/proteins/proteinGroups.txt', sep="\t")

lfq_df = df.set_index("Protein IDs", drop = False)

lfq_df = lfq_df.loc[:,"LFQ intensity 01OV007":"LFQ intensity 17OV026"]

row_list = []
for index, row in lfq_df.iterrows():
        row_list.append(index)

pickle_complex = open("theRealProteinToComplexDict.pickle","rb")
dict_complex = pickle.load(pickle_complex)

#Map the complex as the key with the protein as the value for each complex
complex_to_protein = defaultdict(list)
for k, v in dict_complex.items():
    for i in v:
    	if ((k not in complex_to_protein[i]) and (k in row_list)):
        	complex_to_protein[i].append(k)

complex_to_protein = {k: v for k, v in complex_to_protein.items() if v}

with open('BrentsComplexToProteins.pickle', 'wb') as handle:
    pickle.dump(complex_to_protein, handle, protocol=pickle.HIGHEST_PROTOCOL)
