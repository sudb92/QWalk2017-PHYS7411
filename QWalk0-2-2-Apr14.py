'''
Phase-2 implements a Quantum Random walk on a 2d-lattice under the Grover operator.

B.Sudarsan
18 April 2017
'''
import numpy as np
import matplotlib.pyplot as plt
import math as m
import time

#LSIZE represents the total number of points in the lattice.
#Lattice runs from -N to N including zero
N = 40
LSIZE=2*N+1
COINS = 2

#State of the system, needs two coin and two space indices. Dimensionality is 2x2 times LSIZExLSIZE
state0 = np.zeros(shape=(2,2,LSIZE,LSIZE),dtype=np.complex_)
state1 = np.copy(state0)

#Define coin for the Grover walk.
GCOIN = 0.5*np.ones(shape=(2,2,2,2),dtype=np.complex_)
GCOIN[0,0,0,0] = GCOIN[1,1,1,1]= GCOIN[1,0,1,0]= GCOIN[0,1,0,1] = -1./2.
 
#Check the definition of GCOIN
for i in range(2):
    for j in range(2):
        for k in range(2):
            for l in range(2):
                print GCOIN[i,j,k,l],
        print '\n'

np.set_printoptions(linewidth=1000,suppress=True)

#Start the actual simulation

#Initialize state0 to an appropriate initial state

#BIASED INITIAL CONDITIONS
#Case-1 : Biased initial condition, |coin> is |0,0> to start-off with - left/down bias
#state0[0,0,N,N] = 1 #state0[i,j,0,0] corresponds to lattice-point |-N,-N>, state0[i,j,m,n] to |m-N,n-N>
#print state0

#Case-2 : Biased initial condition, |coin> is |1,1> to start-off with - right/up bias
#state0[1,1,N,N] = 1
#print state0

#Case-3 : Hadamard-Symmetric initial condition, |coin> is 0.5*(|00>+i|01>-i|10>+|11>) to start-off with
#state0[0,0,N,N] = 0.5
#state0[0,1,N,N] = 0.5j
#state0[1,0,N,N] = -0.5j
#state0[1,1,N,N] = 0.5
#print state0

#Case-4 : Grover-Symmetric initial condition, |coin> is 0.5*(|00>+i|01>-i|10>+|11>) to start-off with
state0[0,0,N,N] = 0.5
state0[0,1,N,N] = -0.5
state0[1,0,N,N] = -0.5
state0[1,1,N,N] = 0.5


#Position varies from -N to N via zero, i.e. 0 to LSIZE-1 = 0 to 2N including 2N
#Account for two axes, use a multidimensional array to store the probabilities.
lattice = np.array(range(LSIZE))
xy_probability = np.ndarray(shape=(LSIZE,LSIZE),dtype=float)

def matrix_norm(mat):
    '''
    Function to calculate the norm of an input matrix - the sum of mod of all elements squared
    '''
    norm = 0.0
    for i in np.arange(len(mat)):
        for j in np.arange(len(mat)):
            norm += mat[i,j]
    return norm


def toss(state0,j,k,x,y):
    #Tosses the state-element by matrix multiplication with toss matrix
    out = 0.0 + 0j
    for jp in np.arange(2):
        for kp in np.arange(2):
            out += GCOIN[jp,kp,j,k]*state0[jp,kp,x,y]
    return out

#Calculate the initial probability distribution
for x in lattice:
    for y in lattice:
        xy_probability[x,y] = np.absolute(state0[0,0,x,y])**2 + np.absolute(state0[0,1,x,y])**2 + np.absolute(state0[1,0,x,y])**2 + np.absolute(state0[1,1,x,y])**2

#Set up plot
plt.ion()
fig = plt.figure()
ax = plt.gca()
line1 = ax.imshow(xy_probability,cmap='coolwarm',animated=True,origin='lower',extent=[-N,N,-N,N])
cb = fig.colorbar(line1,ax=ax)
ax.set_xlabel("Lattice position along x")
ax.set_ylabel("Lattice position along y")


#Check that probabilities add to 1
xy_probability_norm = matrix_norm(xy_probability)


#Start the time-loop
timestep = 0
while(timestep<=N):

    norm = 0.0
    for x in lattice:
        for y in lattice:
            for j in np.arange(2):
                for k in np.arange(2):
                    norm += np.conj(state0[j,k,x,y])*state0[j,k,x,y]

    #Update probability array    
    for x in lattice:
        for y in lattice:
            xy_probability[x,y]=0.0
            xy_probability[x,y] = np.absolute(state0[0,0,x,y])**2 + np.absolute(state0[0,1,x,y])**2 + np.absolute(state0[1,0,x,y])**2 + np.absolute(state0[1,1,x,y])**2
    
    #Check that probabilities add to 1
    xy_probability_norm = matrix_norm(xy_probability)
    print xy_probability_norm, norm

    #Find the minimum and maximum to rescale plot with
    arrmin = np.min(xy_probability)
    arrmax = np.max(xy_probability)

    #Update plot data with probability, rescale colorbar
    line1.set_data(xy_probability)
    line1.set_clim(arrmin,arrmax)

    #Set titles, rescale axes
    ax.set_title("Time-step:"+str(timestep)+"   Total Prob:"+str(xy_probability_norm))
    
    #Draw the plot again
    fig.canvas.draw()

    #Apply the toss operators to state0, update the state of the system
    for x in lattice:
        for y in lattice:
            for j in np.arange(2):
                for k in np.arange(2):
                    state1[j,k,x,y] = toss(state0,j,k,x,y)

    state0 = np.copy(state1)

    '''
    for j in np.arange(2):
        for k in np.arange(2):
            if (j,k)==(0,0):
                for x in lattice:
                    for y in lattice:
                        if (x-1<LSIZE)and(y-1<LSIZE)and(x>=1)and(y>=1):
                            state1[j,k,x,y] = state0[j,k,x-1,y-1]
            
            if (j,k)==(0,1):
                for x in lattice:
                    for y in lattice:
                        #print j,k,x,y,state0[j,k,x,y]
                        if (x-1<LSIZE)and(y+1>=1.)and(x-1>=0)and(y+1<LSIZE):
                            state1[j,k,x,y] = state0[j,k,x-1,y+1]
                        #print j,k,x,y,state0[j,k,x,y]
    
            if (j,k)==(1,0):
                for x in lattice:
                    for y in lattice:
                        if (x+1>=0.)and(y-1<LSIZE)and(x+1<LSIZE)and(y-1>=0.):
                            state1[j,k,x,y] = state0[j,k,x+1,y-1]

            if (j,k)==(1,1):
                for x in lattice:
                    for y in lattice:
                        if (x+1>=0.)and(y+1>=0.)and(x+1<LSIZE)and(y+1<LSIZE):
                            state1[j,k,x,y] = state0[j,k,x+1,y+1]
    '''
    #Apply the STEP operator to the state
    for x in lattice:
        for y in lattice:
            for j in np.arange(2):
                for k in np.arange(2):
                    if(x-(-1)**j in lattice) and (y-(-1)**k in lattice):
                            state1[j,k,x,y] = state0[j,k,x-(-1)**j,y-(-1)**k]

    #Increment timestep
    state0= np.copy(state1)            
    timestep = timestep + 1

plt.show(fig)
