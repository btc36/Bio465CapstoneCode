import pickle
from collections import defaultdict

pickle_complex_pathway = open("complexToPathways.pickle","rb")
pickle_complex_pathway_root = open("complexToRootPathways.pickle","rb")

dict_complex_pathway = pickle.load(pickle_complex_pathway)
dict_complex_pathway_root = pickle.load(pickle_complex_pathway_root)

#Map the Pathway as the key with the complex as the value for each pathway
pathway_to_complex = defaultdict(list)
for k, v in dict_complex_pathway.items():
    for i in v:
        pathway_to_complex[i].append(k)
        
#Map the Pathway Root as the key with the complex as the value for each pathway root
pathwayRoot_to_complex = defaultdict(list)
for k, v in dict_complex_pathway_root.items():
    for i in v:
        pathwayRoot_to_complex[i].append(k)

#print(len(pathwayRoot_to_complex))

with open('BrentsPathwayToComplex.pickle', 'wb') as handle:
   pickle.dump(pathway_to_complex, handle, protocol=pickle.HIGHEST_PROTOCOL)
   
with open('BrentsPathwayRootToComplex.pickle', 'wb') as handle:
   pickle.dump(pathwayRoot_to_complex, handle, protocol=pickle.HIGHEST_PROTOCOL)