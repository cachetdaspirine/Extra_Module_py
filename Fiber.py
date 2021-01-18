#!/home/hugo/anaconda3/bin/python3

import numpy as np
import os
from System import * # Module that compute the elastic Energy
from Shape import * # Module that create shapes of 0 and 1

SimNum=2

Wmax=10
LengthRange=20

nurange=0
deltanu=0
nu0=8.
k=1
Kappa=0.1
eps=0.1

#---------------------------------------------------------------
# Delete and remake the Result folder named : Sim"i"
#---------------------------------------------------------------
os.system('rm -rf Sim'+str(SimNum))
os.system('mkdir Sim'+str(SimNum))

#---------------------------------------------------------------
# First loop over the width
#---------------------------------------------------------------
for w in range(2,Wmax):
    #---------------------------------------------------------------
    # Open a new file, notice the 'w' that means write, and thus
    # create a new file if the given file already exist
    #---------------------------------------------------------------
    with open('Sim'+str(SimNum)+'/Energy_width'+str(w)+'.out','w') as myfile:        
        myfile.write('Width Length Np ')
        for i in range(nurange):
            myfile.write('nu'+str(i)+' ')
        myfile.write('\n')
    #---------------------------------------------------------------
    # Loop over the length    
    #---------------------------------------------------------------
    for L in range(LengthRange):
        Length=w+L*w # Set the value of the Length
        #---------------------------------------------------------------
        #---------------------------------------------------------------
            #---------------------------------------------------------------
            # Loop over the value of nu
            #---------------------------------------------------------------
        for k in range(nurange+1) :
        # Output the width/length/number of particle notice the 'a' that
        # means append : thus write, at the end of the existing file
            with open('Sim'+str(SimNum)+'/Energy_width'+str(w)+'.out','a') as myfile:
                myfile.write(str(w)+' '+str(Length)+' '+str(Np(Fiber(w,Length)))+' ')
            nu=nu0+k*deltanu # set the value of nu
            #---------------------------------------------------------------
            # Create the system and output its elastic Energy
            #---------------------------------------------------------------
            print(str(w)+" "+str(Length))
            system=System(Fiber(w,Length),eps=eps,Kmain=k,Kcoupling=Kappa,Kvol=nu)
            with open('Sim'+str(SimNum)+'/Energy_width'+str(w)+'.out','a') as myfile:
                myfile.write(str(system.Energy)+' ')
        with open('Sim'+str(SimNum)+'/Energy_width'+str(w)+'.out','a') as myfile:
            myfile.write('\n')


