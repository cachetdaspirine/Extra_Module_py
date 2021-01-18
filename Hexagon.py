#!/home/hugo/anaconda3/bin/python3
#!/usr/bin/python3

import numpy as np
from System import * # Modul that create the system and return its elastic energy
from Shape import * # Modul that allow to make an Hexagonal shape or a fiber

SimNum=2

SizeMax=50

nurange=1
deltanu=2
nu0=0.
kmain=1.
Kappa=1.
Eps=0.1

os.system('rm -rf Sim'+str(SimNum))
os.system('mkdir Sim'+str(SimNum))

#---------------------------------------------------------------
# Open the output file, and write the name of the column
#---------------------------------------------------------------
with open('Sim'+str(SimNum)+'/Energy.out','w') as myfile:
    myfile.write('Size Np ')
    for i in range(nurange):
        myfile.write('nu'+str(i)+' ')
    myfile.write('\n')

#---------------------------------------------------------------
# Start the loop over the size of the system hexagon of size
# 2n, and 2n+1 has the same size, so size step is 2, started
# with initial size 5 : 6 particle hexagone
#---------------------------------------------------------------
for S in range(5,SizeMax,2):
    #---------------------------------------------------------------
    # write the size, and number of particle
    #---------------------------------------------------------------
    with open('Sim'+str(SimNum)+'/Energy.out','a') as myfile:
        myfile.write(str(S)+' '+str(Np(Parallel(S)))+' ')
    for k in range(nurange) :
        #---------------------------------------------------------------
        # set the value of nu : the volumique energy coefficient
        #---------------------------------------------------------------
        nu=nu0+k*deltanu
        print(S)
        #---------------------------------------------------------------
        # Create the system, the two following line, are about plot the
        # system spring by spring to check if it's what we want
        #---------------------------------------------------------------
        system=System(Parallel(S),eps=Eps,Kmain=kmain,Kcoupling=Kappa,Kvol=nu)
        #system.PrintPerSpring('Sim'+str(SimNum)+'/Size_'+str(S)+'.out')
        #system.PlotPerSpring()
        #---------------------------------------------------------------
        # Output the value of the Energy
        #---------------------------------------------------------------
        with open('Sim'+str(SimNum)+'/Energy.out','a') as myfile:
            myfile.write(str(system.Energy)+' ')
    #---------------------------------------------------------------
    # Next line
    #---------------------------------------------------------------
    with open('Sim'+str(SimNum)+'/Energy.out','a') as myfile:
        myfile.write('\n')
