import pandas as pd
import numpy as np
from pathlib import Path
from scipy import integrate
import matplotlib.pyplot as plt


thisCase=Path('.')
allForcesFiles=[x for x in thisCase.glob('**/*Ux_1.38.*.csv')]
allForcesAsDataframes=[]
blades = []
#Manually impose the header because the file header will be omitted.
header=["Ux","Uy","Uz","p","vtkValidPointMask","arc_length","Points:0","Points:1","Points:2"]

for file in allForcesFiles:
    toDrop=[]
    thisDatabase=pd.read_csv(file)
    #As the parsing introduces NaN when reading '(' or ')' scan the columns and remove them.
    for column in thisDatabase:
          if thisDatabase[column].isnull().values.any():
              toDrop.append(column)
    thisDatabase=thisDatabase.drop(toDrop,axis=1)
    thisDatabase.columns=header
    allForcesAsDataframes.append(thisDatabase)
Umean=np.array([0]*len(allForcesAsDataframes[0]['Ux'].values))
for i in range(0,len(allForcesAsDataframes)):
    Umean=Umean+allForcesAsDataframes[i]['Ux'].values
Umean=Umean/(0.3*len(allForcesAsDataframes))
x=allForcesAsDataframes[0]['arc_length'].values-0.27/2
plt.plot(x,Umean)
