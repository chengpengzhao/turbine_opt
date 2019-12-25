
        #Note: theta must be inserted in degrees. Theta=0 corresponds with this
#        |n           |y
#        |            |
#        |            |
#        ---------t   ------x
        F_xy=np.matrix([Fx,Fy]).T
        F_tn=rotMatx(theta)*F_xy
        return F_tn[0,0],F_tn[1,0]
k
class forcesFile:
    #Class containing the forces as a dataframe, the path of the file and the blade to which it belongs
    filePath=Path('Not initialized')
    def __init__(self,path):
        self.forces=pd.DataFrame()
        self.blade=np.NaN
        self.startTime=np.NaN
        self.endTime=np.NaN
        self.filePath=path
        #Get the character after forcesBlade in the path of the file.
        strPath=str(self.filePath.absolute())
        self.blade=int(strPath[strPath.find('forcesBlade')+len('forcesBlade')])
        self.forces=pd.read_csv(self.filePath, header=None, skiprows=range(3), delimiter="\s+|\(+|\)+",engine='python')
        #As the parsing introduces NaN when reading '(' or ')' scan the columns and remove them.
        toDrop=[]
        for column in self.forces: 
              if self.forces[column].isnull().values.any(): 
                  toDrop.append(column)
        self.forces=self.forces.drop(toDrop,axis=1)
        self.forces.columns=header
        self.startTime=self.forces.Time[0]
        self.endTime=self.forces.Time[len(self.forces)-1]

class bladeForces:
    #Class that extracts and orders data for a specified blade from the list of all data.
    def __init__(self,number,allForcesFiles):
        self.number=int(number) #Number of the blade.
        self.omega=7.4
        self.forces=pd.DataFrame(columns=header)
        self.avgFx=[]
        self.avgFy=[]
        self.avgTz=[]
        self.P=[]
        self.T=2*np.pi/self.omega

        #Create a dataframe containing things from all runs.
        forcesFiles=pd.DataFrame(allForcesFiles,columns=['Objects'])
        #Get start, end and time span, and blade.
        forcesFiles['startTime']=np.nan
        forcesFiles['endTime']=np.nan
        forcesFiles['timeSpan']=np.nan
        forcesFiles['blade']=np.nan
        forcesFiles.startTime=forcesFiles.Objects.apply(lambda x: x.startTime)
        forcesFiles.endTime=forcesFiles.Objects.apply(lambda x: x.endTime)
        forcesFiles.blade=forcesFiles.Objects.apply(lambda x: x.blade)
        forcesFiles.timeSpan=forcesFiles.endTime-forcesFiles.startTime
        #Remove the files not corresponding to this blade.
        forcesFiles=forcesFiles[forcesFiles.blade==self.number]
        #Get different timesteps.
        uniqueStartingTimes=forcesFiles.startTime.unique()
        uniqueStartingTimes.sort()
        rightIdx=[]
        for i in uniqueStartingTimes:
            rightIdx.append(forcesFiles[forcesFiles.startTime==i].timeSpan.idxmax()) #Get the index of the entriy with the same starting time and maximum timespan
            #Now, append all data cutting the previous table to the starting of the next.
        for i in range(len(rightIdx)-1):
            #rows=forcesFiles.iloc[rightIdx[i]].Objects.forces.Time<forcesFiles.iloc[rightIdx[i+1]].startTime #List of bools                
            rows=forcesFiles.loc[rightIdx[i]].Objects.forces.Time<forcesFiles.loc[rightIdx[i+1]].startTime
            self.forces=self.forces.append(forcesFiles.loc[rightIdx[i]].Objects.forces[rows], ignore_index=True)
        #And attach the last file.
        self.forces=self.forces.append(forcesFiles.loc[rightIdx[len(rightIdx)-1]].Objects.forces,ignore_index=True)
        self.calculateAvgF()
        self.calculateAngle()
#        self.calculateAvgTz() Commented out since they may need to be flip first
#        self.calculatePower()
#        self.nonDimensionalize()
            
    def calculateAvgF(self):
        #self.avgFx=self.forces['Force pressure x'].mean()+self.forces['Force viscous x'].mean()
        Fpx=self.forces['Force pressure x']
        Fvx=self.forces['Force viscous x']
        Fpy=self.forces['Force pressure y']
        Fvy=self.forces['Force viscous y']
        F=pd.DataFrame()
        F['Fx']=Fpx+Fvx
        F['Fy']=Fpy+Fvy
        F['Time']=self.forces['Time']
        t=0
        while t<F.Time.iloc[-1]:
            rowsThisRevolution=(F.Time>=t) & (F.Time<(t+self.T))
            self.avgFx.append(F[rowsThisRevolution].Fx.mean())
            self.avgFy.append(F[rowsThisRevolution].Fy.mean())
            t=t+self.T
        #Delete the last, since the revolution won't be complete
        del self.avgFx[len(self.avgFx)-1], self.avgFy[len(self.avgFy)-1] 
        self.avgFx=np.array(self.avgFx)
        self.avgFy=np.array(self.avgFy)
        
    def calculateAvgTz(self):
        #self.avTzy=self.forces['Torque pressure z'].mean()+self.forces['Force viscous z'].mean()
        Tpz=self.forces['Torque pressure z']
        Tvz=self.forces['Torque viscous z']
        T=pd.DataFrame()
        T['Tz']=Tpz+Tvz
        T['Time']=self.forces['Time']
        t=0
        while t<T.Time.iloc[-1]:
            rowsThisRevolution=(T.Time>=t) & (T.Time<(t+self.T))
            self.avgTz.append(T[rowsThisRevolution].Tz.mean())
            t=t+self.T
        #Delete the last, since the revolution won't be complete
        del self.avgTz[len(self.avgTz)-1]
        self.avgTz=np.array(self.avgTz)
        
    def calculatePower(self):
        torque=self.forces[['Time','Torque pressure z', 'Torque viscous z']].copy()
        torque['Torque']=torque['Torque pressure z']+torque['Torque viscous z']
        torque.drop(['Torque pressure z', 'Torque viscous z'], axis=1)
        #Iterate over all revolutions.
        t=0
        while t<torque.Time.iloc[-1]:
            rowsThisRevolution=(torque.Time>=t) & (torque.Time<(t+self.T))
            self.P.append( integrate.simps( torque.Torque[rowsThisRevolution], torque.Time[rowsThisRevolution] )*self.omega )
            t=t+self.T
        #Delete the last, since the revolution won't be complete
        del self.P[len(self.P)-1]
        self.P=np.array(self.P)
        
    def calculateAngle(self):
        #NOTE: this is angle OF THE TURBINE. It is corrected in the next function.
        self.forces['Angle']=self.forces.Time.apply(lambda t: bindAngle( t/self.T*360. ) )
        
    def correctAngle(self,offset):   
        #After the call of this function the Angle inside bladeForces is Angle OF THE BLADE.
        self.forces['Angle']=self.forces.Angle.apply(lambda theta: bindAngle( theta+offset ) )
    
    def nonDimensionalize(self):
        self.r=0.03414
        self.chord=0.0254
        self.u_infty=0.3
        self.h=0.01
        self.gap=0.085
        self.rho=1e3
        self.Force_ref=0.5*self.rho*self.u_infty**2*(2*self.r*self.h)
        self.Torque_ref=self.Force_ref * self.r
        self.P_ref=0.5*self.rho*self.u_infty**3*(2*self.r*self.h) #Pref=0.5 rho u^3 (2rh)
        
        #IMPORTANT NOTE: AVERAGE VALUES OF REVOLUTIONS CORRESPOND TO REVOLUTIONS OF THE TURBINE.
        #This means that they are taken from angles 0-360 for Blade1, from 120-479 for Blade2....        
        self.P=self.P/self.P_ref
        self.avgTz=self.avgTz/self.Torque_ref
        self.avgFx=self.avgFx/self.Force_ref
        self.avgFy=self.avgFy/self.Force_ref
    def flipTorque(self):
        self.forces['Torque pressure z']=-self.forces['Torque pressure z']
        self.forces['Torque viscous z']=-self.forces['Torque viscous z']
    def createLastRev(self):
        t=0
        while t<self.forces.Time.iloc[-1]:
            t=t+self.T
        #The last revolution is incomplete, then we need to take times between t-2*T,t-T
        rowsLastRevolution=(self.forces.Time>=t-2*self.T) & (self.forces.Time<(t-self.T))
        self.forcesLastRev=self.forces[rowsLastRevolution].copy()

    def calcTangentNormalForces(self,kind,position):
        self.forcesLastRev['Fpt']=np.nan
        self.forcesLastRev['Fpn']=np.nan
        self.forcesLastRev['Fvt']=np.nan
        self.forcesLastRev['Fvn']=np.nan
        self.forcesLastRev['Tz']=self.forcesLastRev['Torque pressure z']+self.forcesLastRev['Torque viscous z']
        if kind=='Forward' and position=='Bottom':
            offset=180
            for index, row in self.forcesLastRev.iterrows():
                theta=row.Angle+offset
                self.forcesLastRev.loc[index,['Fpt','Fpn']]=rotateVector(-theta,row['Force pressure x'],row['Force pressure y'])
                self.forcesLastRev.loc[index,['Fvt','Fvn']]=rotateVector(-theta,row['Force viscous x'],row['Force viscous y'])
        elif  kind=='Single' or (kind=='Forward' and position=='Top'):
            offset=90
            for index, row in self.forcesLastRev.iterrows():
                theta=row.Angle+offset
                self.forcesLastRev.loc[index,['Fpn','Fpt']]=rotateVector(theta,row['Force pressure x'],row['Force pressure y'])
                self.forcesLastRev.loc[index,['Fvn','Fvt']]=rotateVector(theta,row['Force viscous x'],row['Force viscous y'])
        elif kind=='Backward' and position=='Top':
            offset=180
            for index, row in self.forcesLastRev.iterrows():
                theta=row.Angle+offset
                self.forcesLastRev.loc[index,['Fpt','Fpn']]=rotateVector(-theta,row['Force pressure x'],row['Force pressure y'])
                self.forcesLastRev.loc[index,['Fvt','Fvn']]=rotateVector(-theta,row['Force viscous x'],row['Force viscous y'])
        elif kind=='Backward' and position=='Bottom':
            offset=90
            for index, row in self.forcesLastRev.iterrows():
                theta=row.Angle+offset
                self.forcesLastRev.loc[index,['Fpn','Fpt']]=rotateVector(theta,row['Force pressure x'],row['Force pressure y'])
                self.forcesLastRev.loc[index,['Fvn','Fvt']]=rotateVector(theta,row['Force viscous x'],row['Force viscous y'])
        
        
class turbineData:
    
    def getRevValues(self):
        #Add the averages of each blade.
        self.avgP=np.zeros(len(self.blades[0].P))
        self.avgFx=self.avgP.copy()
        self.avgFy=self.avgP.copy()
        self.avgTz=self.avgP.copy()
        
        for blade in self.blades:
            self.avgP=self.avgP+blade.P
            self.avgFx=self.avgFx+blade.avgFx
            self.avgFy=self.avgFy+blade.avgFy
            self.avgTz=self.avgTz+blade.avgTz
                
    def __init__(self, threeBlades):
        #Get the path of the file.
        casePath=Path('.')
        self.blades=threeBlades
        #Get the angle correction.
        self.getType(casePath)
        self.getAngleCorrection()
        for i,angle in enumerate(self.correctAngle):
            self.blades[i].correctAngle(angle)
            if self.getTorqueFlip():
                self.blades[i].flipTorque()
            self.blades[i].calculateAvgTz()
            self.blades[i].calculatePower()
            self.blades[i].nonDimensionalize()
        #for i,angle in enumerate(self.blades):
            self.blades[i].createLastRev()
            self.blades[i].calcTangentNormalForces(self.kind,self.position)
        self.getRevValues()
        
        


