import pickle, os, numpy as np
from collections import defaultdict

#Retrieve all dictionaries
pickle_complex = open("BrentsComplexToProteins.pickle","rb")
pickle_p_values = open("ProteinNamesToPValsDictionaryVariance.pickle","rb")

dict_complex = pickle.load(pickle_complex)
dict_p_values = pickle.load(pickle_p_values)

newDictionary = {}

total_count = 0
for key, value in dict_complex.items():
    total_count = 0
    for p_value in value:
        if p_value in dict_p_values:
            if dict_p_values[p_value] < 0.05:
                total_count += 1
    print(total_count/len(value))
    newDictionary[key] = (total_count/len(value))

total_count = 0
for key, value in newDictionary.items():
    if value > 0:
        print(key)
        total_count += 1
        
print(total_count)
print(dict_complex)

with open('ComplexDifferentialVariance.pickle', 'wb') as handle:
   pickle.dump(newDictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)