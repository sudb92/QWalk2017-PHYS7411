'''
Program to print pretty plots of the 1D QWalk simulations

B.Sudarsan
March 14,2017
'''
import numpy as np
import matplotlib.pyplot as plt

fsim = open("QWalkSimulation(N=100),symmetric.dat","r")
#fsim = open("QWalkSimulation(N=100),right-biased.dat","r")

#ftheo = open("QWalkTheoretical(N=100),left-biased.dat","r")

#Read the files, skip the first line with titles
linesim = fsim.readlines()
linesim = linesim[1:]

#Prepare container variables to store data in
xsim = np.empty([])
probsim = np.empty([])
errsim = np.empty([])

#Process the string data read from strings, box them into arrays
for line2 in linesim:
    line21 = line2.split()
   
    try:
        xsim = np.append(xsim,float(line21[0]))
        probsim = np.append(probsim,float(line21[1]))
    
    except ValueError,e:
        print 'error',e,'on line'

'''
#Plot the results
plt.plot(xtheo,probtheo,label='Theoretical result')
plt.plot(xsim,probsim,label='Simulation result')
plt.title("Comparison between theory and simulation, 1D random walk")
plt.legend(loc='lower right')
plt.show()

#Plot the deviation
chi2_mag = np.dot(chi2.T,chi2)
plt.plot(xsim,chi2,label='chi2 between results')
plt.title("Chi2 deviation between results:"+str(chi2_mag))
plt.show()
'''
print xsim[1::2], xsim
plt.plot(xsim[1::2],probsim[1::2])#,label='Simulation result')
#plt.title("$n=100, Initial$ $conditions=(Coin=(|0>+i|1>)/\sqrt{2}, Lattice=|0>)$")
plt.title("$n=100, Initial$ $conditions=(Coin=|1>, Lattice=|0>)$")
plt.legend(loc='lower right')
plt.xlabel('Position, x')
plt.ylabel('Probability P(x,t=100)')
plt.show()
