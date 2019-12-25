import pandas as pd
import numpy as np
from pathlib import Path
from scipy import integrate
import matplotlib.pyplot as plt
#This script loads every forces file into pandas dataframes. Databases are stored at allForcesAsDataframes and their original path at allForcesFiles.

def bindAngle(theta):
    #Bind between 0 and 360
    while theta>=360.:
        theta=theta-360.
    while theta<0:
        theta=theta+360.
    return theta
#Get all the forces files.
thisCase=Path('.')
allForcesFiles=[x for x in thisCase.glob('**/*forces.dat')]
allForcesAsDataframes=[]
blades = []
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
    blades.append(int(str(file.absolute())[str(file.absolute()).find('forcesBlade')+len('forcesBlade')]))

# m  s
Rairfoil=0.03414
Lchord=0.0254
length_z=0.001
rho=1e3
Uinf=0.3
omega=7.4 #rad/s
offset=[0,120,-120]
A=2*Rairfoil*length_z
Tavg=np.array([0]*len(thisDatabase))
for i in range(0,len(blades)):
    vars()['blade'+str(blades[i])+'forces']=allForcesAsDataframes[i]
    vars()['blade'+str(blades[i])+'_t']=allForcesAsDataframes[i]['Time'].values
    vars()['blade'+str(blades[i])+'_Fpx']=allForcesAsDataframes[i]['Force pressure x'].values
    vars()['blade'+str(blades[i])+'_Fpy']=allForcesAsDataframes[i]['Force pressure y'].values
    vars()['blade'+str(blades[i])+'_Fvx']=allForcesAsDataframes[i]['Force viscous x'].values
    vars()['blade'+str(blades[i])+'_Fvy']=allForcesAsDataframes[i]['Force viscous y'].values
    vars()['blade'+str(blades[i])+'_Fx']=vars()['blade'+str(blades[i])+'_Fpx']+vars()['blade'+str(blades[i])+'_Fvx']
    vars()['blade'+str(blades[i])+'_Fy']=vars()['blade'+str(blades[i])+'_Fpy']+vars()['blade'+str(blades[i])+'_Fvy']
    vars()['blade'+str(blades[i])+'_Tpz']=allForcesAsDataframes[i]['Torque pressure z'].values
    vars()['blade'+str(blades[i])+'_Tvz']=allForcesAsDataframes[i]['Torque viscous z'].values
    vars()['blade'+str(blades[i])+'_Tz']=vars()['blade'+str(blades[i])+'_Tpz']+vars()['blade'+str(blades[i])+'_Tvz']
    vars()['blade'+str(blades[i])+'_theta']=vars()['blade'+str(blades[i])+'_t']*omega*180/np.pi+offset[i]

for i in range(len(blade1_t)):
    if blade1_t[i] < 2 * np.pi / omega and blade1_t[i + 1] > 2 * np.pi / omega:
        tindex = i

for i in range(0,len(blades)):
    vars()['blade' + str(blades[i]) + '_Tavg'] = [0] * len(vars()['blade' + str(blades[i]) + '_t'])
    for j in range(tindex, len(vars()['blade' + str(blades[i]) + '_t']) + 1):
        t = vars()['blade' + str(blades[i]) + '_t'][j - tindex:j]
        Tz = vars()['blade' + str(blades[i]) + '_Tz'][j - tindex:j]
        vars()['blade' + str(blades[i]) + '_Tavg'][j - 1] = integrate.simps(Tz, t) / (t[-1] - t[0])
    Tavg = Tavg + np.array(vars()['blade' + str(blades[i]) + '_Tavg'])

Cp=Tavg*omega/(0.5*rho*Uinf**3*A)
plt.plot(blade1_t[tindex:-1],Cp[tindex:-1])
#plt.plot(blade1_t,blade1_Tz+blade2_Tz+blade3_Tz)
