import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle

pickle_complex = open("BrentsComplexToProteins.pickle","rb")
dict_complex = pickle.load(pickle_complex)
df = pd.DataFrame.from_dict(dict_complex, orient='index')
print(df)
sns.set()
#ax = sns.heatmap(df)
#plt.show()
