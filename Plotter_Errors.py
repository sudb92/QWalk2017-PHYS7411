'''
Another program to print pretty plots of the 1D QWalk simulations for Errors

B.Sudarsan
March 28,2017
'''
import numpy as np
import matplotlib.pyplot as plt

#fsim = open("Error-checker-Biased-Sim-analytic-N=50.dat","rb")

out = np.loadtxt("Error-checker-Biased-Sim-analytic-N=50.dat")

#print out[:,0]

timestep = out[:,0]
probr = out[:,1]
probi = out[:,2]
probsimr = out[:,3]
probsimi = out[:,4]
errr = out[:,5]
erri = out[:,6]
errsimr = out[:,7]
errsimi = out[:,8]

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
plt.plot(timestep,probr,'x-',label='Analytic,Re(P)')
plt.plot(timestep,probi,'x-',label='Analytic,Im(P)')
plt.plot(timestep,probsimr,label='Simulation,Re(P)')
plt.plot(timestep,probsimi,label='Simulation,Im(P)')
plt.title("Tot. probability P versus time, initial state: Coin=|0>, Lattice=|0>")
plt.legend(loc='lower right')
plt.xlabel('Timestep t')
plt.ylabel('Total probability')
ax = plt.gca()
ax.set_ylim([-1.5,1.5])
plt.show()

plt.plot(timestep,errr,'x-',label='Re(Integ.err) from quad()')
plt.plot(timestep,erri,'x-',label='Im(Integ.err) from quad()')
plt.plot(timestep,errsimr,label='Re|Psi(sim)-Psi(analyt)|')
plt.plot(timestep,errsimi,label='Im|Psi(sim)-Psi(analyt)|')
plt.title("Error vs time, initial state: Coin=|0>, Lattice=|0>")
plt.legend(loc='lower right')
plt.xlabel('Timestep t')
plt.ylabel('Error')
ax = plt.gca()
ax.set_ylim([-1.5e-14,1.5e-14])
plt.show()
