'''
Program to numerically plot the theoretical model of the 1D quantum walk
with an initial bias to the left-side, and the actual simulation on top.

Uses quad() to perform integration over real and imaginary part of both the halves of the wavefunction.

B.Sudarsan
March 28, 2017
'''
import scipy
from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt
import time
from lib import *

#The real and imaginary parts of the inverse Fourier transform, written out in detail
def Psi_Left_Re(k,n,t):
    omegak = np.arcsin(np.sin(k)/np.sqrt(2.0))
    out = 0.5*(1 + (-1)**(t+n))*(1.0/(2.0*np.pi))*(1.0+(np.cos(k)/np.sqrt(1+np.cos(k)**2)))*np.cos(omegak*t+k*n)
    return out

def Psi_Left_Im(k,n,t):
    omegak = np.arcsin(np.sin(k)/np.sqrt(2.0))
    out = -0.5*(1 + (-1)**(t+n))*(1.0/(2.0*np.pi))*(1.0+(np.cos(k)/np.sqrt(1+np.cos(k)**2)))*np.sin(omegak*t+k*n)
    return out

def Psi_Right_Re(k,n,t):
    omegak = np.arcsin(np.sin(k)/np.sqrt(2.0))
    out = 0.5*(1 + (-1)**(t+n))*(1.0/(2.0*np.pi))*(1.0/np.sqrt(1+np.cos(k)**2))*np.cos(k-omegak*t-k*n)
    return out

def Psi_Right_Im(k,n,t):
    omegak = np.arcsin(np.sin(k)/np.sqrt(2.0))
    out = 0.5*(1 + (-1)**(t+n))*(1.0/(2.0*np.pi))*(1.0/np.sqrt(1+np.cos(k)**2))*np.sin(k-omegak*t-k*n)
    return out

#Test routines
#print quad(Psi_Left_Re, -np.pi,np.pi,args=(2,52))
#print quad(Psi_Left_Im, -np.pi,np.pi,args=(2,52))
#print quad(Psi_Right_Re, -np.pi,np.pi,args=(2,52))
#print quad(Psi_Right_Im, -np.pi,np.pi,args=(2,52))

'''
Begin actual simulation of QWalk in 1D.
'''

#LSIZE represents the total number of points in the lattice.
#Lattice runs from -N to N including zero
N = 50
LSIZE=2*N+1
COINS = 1

#State of the system, initialize with zeros on every element.
state0 = np.zeros(shape=(LSIZE*COINS*2,1),dtype=complex)

def return_array_indices(lattice_site):
    '''
    Function to return the array indices corresponding to a given lattice site |x> given as input.
    This would allow us to compute probabilities easily.
    '''
    x = lattice_site
    indices_list = np.array([])
    i = 1
    for coin_id in range(0,2*COINS):    #Leads to 0,1 for 1 coin, 0,1,2,3 for 2 coins etc.
        indices_list = np.append(indices_list, x + N + coin_id*(2*N+1))
    return (indices_list)

'''
Things needed for simulation.
'''


#State of the system
state0_sim = np.zeros(shape=(LSIZE*COINS*2,1),dtype=complex)
locs = return_array_indices(0)
state0_sim[locs[1]] = 1


#TOSS is the cointoss operator, STEP is the conditional step operator
TOSS = np.zeros(shape=(LSIZE*COINS*2,LSIZE*COINS*2),dtype=complex)
STEP = np.zeros(shape=(LSIZE*COINS*2,LSIZE*COINS*2),dtype=complex)

for row in range(len(state0)):
    for col in range(len(state0)):
        if (row == col) and (row<LSIZE): #Submatrix A's diagonal in figure above
            TOSS[row,col] = 1/np.sqrt(2)
        elif (row == col) and (row>=LSIZE): #Submatrix D's diagonal
            TOSS[row,col] = -1/np.sqrt(2)
        elif(row == col + LSIZE):              #Submatrix B's diagonal
            TOSS[row,col] = 1/np.sqrt(2)
        elif(col == row + LSIZE):            #Submatrix D's diagonal
            TOSS[row,col] = 1/np.sqrt(2)
        else:
            TOSS[row,col] = 0

        if (row == col +1) and (row < LSIZE): #Submatrix A's subdiagonal elements below
            STEP[row,col] = 1
        elif (col == row +1) and (row >= LSIZE): #Submatrix D's subdiagonal elements above
            STEP[row,col] = 1    
        else:
            STEP[row,col]=0

#print TOSS
#quit()

#Position varies from -N to N via zero, i.e. 0 to LSIZE-1 = 0 to 2N including 2N
#Define a probability vector to print
#Define an error
lattice = range(LSIZE)
x_probability = np.ndarray(shape=(LSIZE,1),dtype=float)
error_vec = np.ndarray(shape=(LSIZE,1),dtype=float)

x_probability_sim = np.ndarray(shape=(LSIZE,1),dtype=float)
error_vec_sim = np.ndarray(shape=(LSIZE,1),dtype=float)


ferr = open("Error-checker-Biased-Sim-analytic-N="+str(N)+".dat","ab")

#Setup output figure properties
plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111)
meanlabel = ax.annotate(' ',xy=(N/2,1.0),xytext=(N/2,1.0))
line1, = ax.plot(np.arange(-N,N+1), x_probability,'r-',label='P(x,t)-analyt')
line2, = ax.plot(np.arange(-N,N+1), x_probability_sim,'g-',label='P(x,t)-simul')
line3, = ax.plot(np.arange(-N,N+1), error_vec_sim,'b-',label='difference')
plt.legend(loc='upper right')

#Loop in time!
timestep = 0
while(timestep<=N):
    
    for x in lattice:
        #Find array indices corresponding to the lattice position
        locs = return_array_indices(x-N)
        n_left = locs[0]
        n_right = locs[1]
        
        #Calculate the integrals corresponding to each position
        LR,LRerr = quad(Psi_Left_Re, -np.pi,np.pi,args=(x-N,timestep))
        LI,LIerr = quad(Psi_Left_Im, -np.pi,np.pi,args=(x-N,timestep))
        RR,RRerr = quad(Psi_Right_Re, -np.pi,np.pi,args=(x-N,timestep))
        RI,RIerr = quad(Psi_Right_Im, -np.pi,np.pi,args=(x-N,timestep))

        #Store the values in the state vector
        state0[n_left] = LR + 1j*LI
        state0[n_right] = RR + 1j*RI

        #Dump all the computational errors due to the quad() calls into the error_vec at the current location.
        error_vec[x] = np.sqrt(np.absolute(LRerr)**2 + np.absolute(LIerr)**2 + np.absolute(RRerr)**2 + np.absolute(RIerr)**2)

        #Store the probability of each position separately in the requisite array.
        x_probability[x] = np.absolute(state0[locs[0]])**2 + np.absolute(state0[locs[1]])**2

        #Do it for the simulation as well
        x_probability_sim[x] = np.absolute(state0_sim[locs[0]])**2 + np.absolute(state0_sim[locs[1]])**2
        

    #Calculate both the probability norm as well as the error norm to ensure that there is no buildup of computational error
    x_probability_norm = np.dot(np.conj(state0.T),state0)
    error_vec_norm = np.dot(np.conj(error_vec.T),error_vec)

    x_probability_sim_norm = np.dot(np.conj(state0_sim.T),state0_sim)

    error_vec_sim = x_probability - x_probability_sim
    error_vec_sim_norm = np.sqrt(np.dot(np.conj(error_vec_sim.T),error_vec_sim))

    #Write probabilities and error parameters to output file
    outline=np.array([timestep,np.real(x_probability_norm),np.imag(x_probability_norm),np.real(x_probability_sim_norm),np.imag(x_probability_sim_norm),np.real(error_vec_norm),np.imag(error_vec_norm),np.real(error_vec_sim_norm),np.imag(error_vec_sim_norm)])
    np.savetxt(ferr,[outline],fmt='%1.10e')
    
    
    #time.sleep(0.1)    #As it is, the running speed is low!
    line1.set_ydata(x_probability)
    line2.set_ydata(x_probability_sim)
    line3.set_ydata(error_vec_sim)
    ax.set_title("t:"+str(timestep)+" Prob.sum(theo,sim):"+str(x_probability_norm)+","+str(x_probability_sim_norm))
    meanlabel.xy = (0,max(x_probability)/2)
    meanlabel.xytext = (0,max(x_probability)/2)
    meanlabel.set_text("Intgn err:"+str(error_vec_norm)+"\nd(simul-analytic):"+str(error_vec_sim_norm))
    ax.set_xlabel("Position in lattice")
    ax.set_ylabel("Probability P(x,t)")
    ax.relim()
    ax.autoscale_view()
    fig.canvas.draw()

    #Update simulated state's state
    state0_sim = np.dot(TOSS,state0_sim)
    state0_sim = np.dot(STEP,state0_sim)
    
    timestep = timestep + 1

f = open("QWalkTheoretical(N="+str(N)+"),left-biased.dat","w")
f.write("#x \t prob[x] \t error_vec[x]")
for i in np.arange(len(x_probability)):
    f.write("\n"+str(i-N)+"\t"+str(x_probability[i][0])+"\t"+str(error_vec[i][0]))
f.close()
ferr.close()                           
plt.show(fig)
quit()






