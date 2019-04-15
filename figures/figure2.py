import seaborn as sns
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import mplcursors
import webbrowser

def saveI(i):
	print(i)
	#return i
        	
def refreshGraph():
    nx.draw(G, pos=pos, show_labels=False, node_size=8)
    plt.axis((-1,1,-1,1))
    
def onClick(event):
    (x,y) = (event.xdata, event.ydata)

    for i in pathways:            
        node = pos[i]
        distance = pow(x-node[0],2)+pow(y-node[1],2)
        if distance < 0.0001:
        	print(i)
        	#webbrowser.open('http://google.com')
			
def onHover(event):
	(x,y) = (event.xdata, event.ydata)
	
	mplcursors.cursor(hover = True).connect("add", lambda sel: sel.annotation.set_text(pathways[sel.target.index]))
			
df = pd.read_csv('ReactomePathwaysRelation.txt', sep="\t")

headers = ["Parent","Child"]
df.columns = headers

# keep only HSA rows --------------------------------
Parents = df["Parent"]
Children = df["Child"]

# grab the names of R-HSA- pathways
pathways = []
for pathway in Parents:
  if pathway.startswith("R-HSA-"):
    pathways.append(pathway)
for pathway in Children:
  if pathway.startswith("R-HSA-"):
    pathways.append(pathway)
    
# make unique
pathways = list(set(pathways))

#print(pathways)

# remove non-R-HSA pathways
df = df[df["Parent"].isin(pathways)]
df = df[df["Child"].isin(pathways)]
  
fig, ax = plt.subplots()
fig.canvas.mpl_connect('button_press_event', onClick)
fig.canvas.mpl_connect('motion_notify_event', onHover)

# make the graph
G = nx.Graph()

# add nodes -------------------------------------------
G.add_nodes_from(pathways)

# add edges -------------------------------------------
for index, row in df.iterrows():
  G.add_edge(row[0],row[1])

pos = nx.spring_layout(G)

#print(pos)

refreshGraph()
    
plt.show()