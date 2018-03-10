'''
Program simulating a 1D quantum walk on a line, for different initial conditions

B.Sudarsan
17 Feb 2017
'''

import numpy as np
import matplotlib.pyplot as plt
import math as m
import time

#LSIZE represents the total number of points in the lattice.
#Lattice runs from -N to N including zero

N = 100
LSIZE=2*N+1
COINS = 1

#State of the system
state0 = np.zeros(shape=(LSIZE*COINS*2,1),dtype=complex)

#TOSS is the cointoss operator, STEP is the conditional step operator
TOSS = np.zeros(shape=(LSIZE*COINS*2,LSIZE*COINS*2),dtype=complex)
STEP = np.zeros(shape=(LSIZE*COINS*2,LSIZE*COINS*2),dtype=complex)
'''
0th element maps to lattice point -10 means 10 maps to 0th lattice point, 20 maps to
10th latpoint, with cointoss 0 basis. 21 maps to latticepoint -10 again, `32 maps to
0th latpoint, 42 maps to +10.

'''

#print len(state0)
#print len(TOSS)
#print len(STEP)


def return_array_indices(lattice_site):
    '''
    Function to return the array indices corresponding to a given lattice site |x> given as input.
    This would allow us to compute probabilities easily.
    '''
    x = lattice_site
    indices_list = np.array([])
    i = 1
    for coin_id in range(0,2*COINS):    #Leads to 0,1 for 1 coin, 0,1,2,3 for 2 coins etc.
        indices_list = np.append(indices_list, int(x + N + coin_id*(2*N+1)))
    return (indices_list)

#Define operators TOSS and STEP
row=0
col=0
TOSS[row,col]=1

#print TOSS

'''
[A | B]
-------
[C | D}

Fig.1) The submatrices involved in operators used, STEP and TOSS.
'''

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
#'''
#Testing the initialized matrices
np.set_printoptions(linewidth=100,suppress=True)

#Start the actual simulation


#print state0 #Will be all zeros

#Initialize state0 to an appropriate initial state

#Case-1 : Biased initial condition, |coin> is |0> to start-off with
#locs = return_array_indices(0)
#state0[np.int(locs[0])] = 1
#print state0

#Case-2 : Biased initial condition, |coin> is |1> to start-off with
#locs = return_array_indices(0)
#state0[np.int(locs[1])] = 1
#print state0

#Case-3 : Symmetric initial condition, |coin> is (|0>+i|1>)/sqrt(2) to start-off with
locs = return_array_indices(0)
state0[np.int(locs[0])] = 1/np.sqrt(2)
state0[np.int(locs[1])] = (0+1j)/np.sqrt(2)
print state0

#Position varies from -N to N via zero, i.e. 0 to LSIZE-1 = 0 to 2N including 2N
lattice = range(LSIZE)
x_probability = np.ndarray(shape=(LSIZE,1),dtype=float)

for x in lattice:
    # x runs from 0 to 2N, x-N runs from -N to N via 0
    locs = return_array_indices(x-N)
    x_probability[x] = np.absolute(state0[np.int(locs[0])])**2 + np.absolute(state0[np.int(locs[1])])**2
    #print x_probability[x]

#Set up plot
plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(np.arange(-N,N+1), x_probability,'r-')

#Check that probabilities add to 1
x_probability_norm = np.dot(np.conj(state0.T),state0)
print x_probability_norm


#Start the time-loop
timestep = 0
while(timestep<=N):
    
    for x in lattice:
        locs = return_array_indices(x-N)
        x_probability[x] = np.absolute(state0[np.int(locs[0])])**2 + np.absolute(state0[np.int(locs[1])])**2

    #Check that probabilities add to 1
    x_probability_norm = np.dot(np.conj(state0.T),state0)

    #Plot slowly, so as to improve visibility
    time.sleep(0.1)
    #Update plot-line

    line1.set_ydata(x_probability)
    #Set titles, rescale axes
    ax.set_title("Time-step:"+str(timestep)+"   Total Prob:"+str(x_probability_norm[0,0]))
    ax.relim()
    ax.autoscale_view()

    #Draw the plot again
    #fig.canvas.draw()
    plt.draw()
    
    #if timestep == 0:
    #    break
    plt.show(fig)

    #Apply the timestep operators to state0, update the state of the system
    state0 = np.dot(TOSS,state0)
    state0 = np.dot(STEP,state0)
    timestep = timestep + 1

plt.show(fig)

f = open("QWalkSimulation(N="+str(N)+"),symmetric.dat","w")
f.write("#x \t prob[x] \t prob_norm")
for i in np.arange(len(x_probability)):
    f.write("\n"+str(i-N)+"\t"+str(x_probability[i][0])+"\t"+str(x_probability_norm[0]))
f.close()

quit()







        
