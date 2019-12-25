import pandas as pd
import numpy as np
from pathlib import Path
#This script loads every forces file into pandas dataframes. Databases are stored at allForcesAsDataframes and their original path at allForcesFiles.

#Get all the forces files.
thisCase=Path('.')
allForcesFiles=[x for x in thisCase.glob('**/*forces.dat')]
allForcesAsDataframes=[]
#Manually impose the header because the file header will be omitted.
header=['Time', 'Force pressure x','Force pressure y','Force pressure z', 'Force viscous x','Force viscous y','Force viscous z','Force porous x','Force porous y','Force porous z', 'Torque pressure x','Torque pressure y','Torque pressure z', 'Torque viscous x','Torque viscous y','Torque viscous z','Torque porous x','Torque porous y','Torque porous z']

for file in allForcesFiles:
    toDrop=[]
    thisDatabase=pd.read_csv(file, header=None, skiprows=range(3), delimiter="\s+|\(+|\)+",engine='python')
    #As the parsing introduces NaN when reading '(' or ')' scan the columns and remove them.
    for column in thisDatabase: 
          if thisDatabase[column].isnull().values.any(): 
              toDrop.append(column)
    thisDatabase=thisDatabase.drop(toDrop,axis=1)
    thisDatabase.columns=header
    allForcesAsDataframes.append(thisDatabase)

