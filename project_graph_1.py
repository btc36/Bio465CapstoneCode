import seaborn as sns
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
import os
import networkx as nx
import pylab
import webbrowser
import mplcursors
from collections import defaultdict
from sklearn.decomposition import PCA

def checkStar(protein_name):

	protein_name_new = id_to_name.get(protein_name)

	if dict_p_values.get(protein_name_new) == None:
		return protein_name_new
	elif dict_p_values.get(protein_name_new) == 1:
		return protein_name_new
	elif dict_p_values.get(protein_name_new) < 0.001:
		return protein_name_new + "**"
	else:
		return protein_name_new + "*"

def refreshGraph():
    plt.clf()
    nx.draw_networkx_nodes(T, pos, nodelist=test_array, node_color='k', node_size=400, alpha=0.8)
    nx.draw_networkx_edges(T,pos,width=1.0, alpha=0.5)
    plt.axis('off')
    plt.axis((-8,8,-8,8))
    fig.patch.set_facecolor('white')
    plt.show()

def onClick(event):
        (x,y) = (event.xdata, event.ydata)
        for i in test_array:
                node = pos[i]
                distance = pow(x-node[0],2)+pow(y-node[1],2)
                if distance < 0.1:
                        webbrowser.open('http://google.com')
                        refreshGraph()

def getLFQData():
        df = pd.read_csv('/Users/brent/Desktop/proteins/proteinGroups.txt', sep="\t")
        lfq_df = df.set_index("Protein IDs", drop = False)
        return(lfq_df.loc[:,"LFQ intensity 01OV007":"LFQ intensity 17OV026"])

def createCancerNonCancerArrays(df):
        arrayCancer = []
        arrayNonCancer = []

        for col in df:
                if '_NM' not in col:
                        arrayCancer.append(col)
                else:
                        arrayNonCancer.append(col)
        return(arrayCancer, arrayNonCancer)

lfq_df = getLFQData()
arrayCancer, arrayNonCancer = createCancerNonCancerArrays(lfq_df)

#Retrieve all dictionaries
pickle_complex = open("BrentsComplexToProteins.pickle","rb")
pickle_p_values = open("ProteinNamesToPValsDictionary.pickle","rb")
pickle_happy_complex = open("TheHappyComplexFamilies.pickle", "rb")
pickle_complex_pathway = open("complexToPathways.pickle", "rb")
pickle_id_to_name = open("uniprotToGene.pickle", "rb")
pickle_pathwayRoot_complex = open("BrentsPathwayRootToComplex.pickle", "rb")
pickle_pathway_complex = open("BrentsPathwayToComplex.pickle", "rb")
pickle_allPathways = open("allPathways.pickle", "rb")

dict_complex = pickle.load(pickle_complex)
dict_p_values = pickle.load(pickle_p_values)
happy_complex = pickle.load(pickle_happy_complex)
complex_pathway = pickle.load(pickle_complex_pathway)
id_to_name = pickle.load(pickle_id_to_name)
pathwayRoot_complex = pickle.load(pickle_pathwayRoot_complex)
pathway_complex = pickle.load(pickle_pathway_complex)
allPathways = pickle.load(pickle_allPathways)

'''for key, value in pathway_complex.items():
        directory = "/Users/brent/Desktop/proteins/pathway/"+key
        if not os.path.exists(directory):
                for value_in_array in value:
                        directory = "/Users/brent/Desktop/proteins/pathway/"+key+"/"+value_in_array
                        if not os.path.exists(directory):
                                os.makedirs(directory)'''

#Each keyP will be a node on the network. Click on the node and it should bring
#up the corresponding complexes within the path. Need to draw an edge to all connecting pathways
#pathway_list = []
#pathway_index = []

#for keyP, valueP in pathway_complex.items():
#        pathway_list.append(keyP)

#for i in range(0, len(pathway_list)):
#        pathway_index.append(i)
print(len(pathway_complex))
run_size = 1
for keyP, valueP in pathway_complex.items():
        print(keyP, run_size)
        run_size += 1
        for key, value in dict_complex.items():
                directory = "/Users/brent/Desktop/proteins/pathway/"+keyP+"/"+key
                if not os.path.exists(directory):
                        continue
                if os.path.exists(directory+'/'+key+'.png'):
                        continue

                #test_array.append(int_i)
                value_array_cancer = []
                value_array_non_cancer = []
                protein_array = []
                #print(key, value)
                for value_in_array in value:
                        #print(value_in_array)
                        for j in lfq_df.loc[value_in_array, arrayCancer]:
                                value_array_cancer.append(j)
                        for d in lfq_df.loc[value_in_array, arrayNonCancer]:
                                value_array_non_cancer.append(d)
                        protein_array.append(value_in_array)
                        
                #Created these three arrays for the transform..
                ultimate_value_array = []
                ultimate_name_array = []
                ultimate_canc_nonCanc_array = []

                #Add values to the arrays
                for i in protein_array:
                        for j in lfq_df.loc[i, arrayCancer]:
                            if j > 0:
                                ultimate_value_array.append(j)
                                ultimate_name_array.append(i)
                                ultimate_canc_nonCanc_array.append('Yes')
                        for l in lfq_df.loc[i, arrayNonCancer]:
                            if l > 0:
                                ultimate_value_array.append(l)
                                ultimate_name_array.append(i)
                                ultimate_canc_nonCanc_array.append('No')

                #Set theme
                sns.set_style('darkgrid')

                ultimate = {'Protein IDs': ultimate_name_array,
                         'Protein Abundance': ultimate_value_array,
                         'Cancer': ultimate_canc_nonCanc_array}

                if (len(ultimate.get('Protein IDs')) == 0) or len(ultimate.get('Protein Abundance')) == 0 or len(ultimate.get('Cancer')) == 0:
                        continue

                df_ultimate = pd.DataFrame(data = ultimate)

                df_ultimate = df_ultimate.loc[(df_ultimate["Cancer"] == "No") | (df_ultimate["Cancer"] == "Yes")]

                df_ultimate = df_ultimate.dropna(axis=0)

                fig, ax = plt.subplots(figsize=(10, 6))

                my_pal = {"Yes": "#003CFF", "No": "#ffffff"}

                ax = sns.boxplot(x = "Protein IDs", y = "Protein Abundance",
                                     data = df_ultimate, hue = 'Cancer',
                                     palette=my_pal, linewidth = 0.5,showfliers=False)

                ax = sns.stripplot(data = df_ultimate, x= 'Protein IDs',
                                       y='Protein Abundance', hue='Cancer', size = 1.8,
                                       dodge=True, jitter=True, color='.3')

                protein_array_name = []
                for i in protein_array:
                        protein_array_name.append(checkStar(i))
                #print(protein_array_name)
                ax.set_xticklabels(protein_array_name)
                ax.set_ylabel('')
                ax.set_xlabel('')
                ax.tick_params(labelsize='5')

                #Adjust legend
                handles, labels = ax.get_legend_handles_labels()
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

                # Put a legend to the right of the current axis
                ax.legend()
                plt.legend(handles[0:3], ['Yes', 'No'], title='Cancer', fontsize='13', frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))
                print("saved for: " + keyP + " " + key)
                fig.savefig(directory+'/'+key+'.png', dpi=300)

'''fig1, ax1 = plt.subplots()
fig1.canvas.mpl_connect('button_press_event', onClick)

T = nx.Graph()
for node in test_array:
	T.add_node(node)

pos={}
pos[1]=np.array([0,0])
pos[2]=np.array([-2,1])
pos[3]=np.array([2,1])
pos[4]=np.array([-3,2])
pos[5]=np.array([-1,2])
pos[6]=np.array([1,2])
pos[7]=np.array([6,2])
pos[8]=np.array([-3,2])
pos[9]=np.array([5,2])
pos[10]=np.array([4,2])

### Draw nodes and edges
refreshGraph()'''
