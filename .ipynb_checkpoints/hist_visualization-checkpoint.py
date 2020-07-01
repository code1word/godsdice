import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import pageRankSim

# Some example run code for pageRankSim
"""
#initialize run conditions
n_runs = 10 # number of runs in final iteration
shots = 50 # number of shots on IBM qasm simulator
n_frames = 10 # number of frames in gif n_runs >= n_frames or else it wont work

#Adjacency Matrix

A = np.array([  [0,1,0,0,1,1,0,0],
                [0,0,0,0,0,0,0,0],
                [1,1,0,0,0,0,1,0],
                [0,0,1,0,1,1,0,0],
                [0,0,0,0,0,0,1,0],
                [0,0,0,1,0,0,0,0],
                [0,0,0,0,0,1,0,0],
                [0,0,0,0,0,0,0,0]])

# for graph node 8 disconnected : index = 7
disc_indices = [7]

nqubits=int(np.log2(A.shape[0]))
disc_indices2 = np.zeros(len(disc_indices))

# Using indices from disconnected vertices to decide complementary indices from
# left-concatenating 1 to the binary number
for i in range(len(disc_indices)):
    disc_indices2[i]=disc_indices[i] + 2**nqubits

n_states = 2**nqubits - len(disc_indices)
initsize = A.shape[0]*2
init_vector = np.ones(initsize)/np.sqrt(initsize-2*len(disc_indices))

    # making parts of init vector =0
for i in range(len(disc_indices)):
    init_vector[int(disc_indices[i])] = 0
    init_vector[int(disc_indices2[i])] = 0

# initialize to superposition without |up>|111> or |down>|111>

rankmatrix=np.zeros((n_frames,n_states))

for i in range(1,n_frames+1,1):
    rankrow = pageRankSim.simulate(A,init_vector,int((n_runs/n_frames)*i),shots)
    rankmatrix[i-1,:]=rankrow
"""

# fname should be .mp4!
def hist_viz(rankmatrix, fps, fname):
    fig, ax = plt.subplots()

    xpos = np.arange(1,rankmatrix[0].shape[0]+1)
    ax.set_xticks( xpos, [ str(x) for x in xpos ]) # name them 1 - 7

    hist_ylim = np.min([ np.around(np.max(rankmatrix), decimals=1) + 0.05, 1.0 ])

    # animation
    def update_hist(row, data):
        ax.clear()
        ax.set_ylim(0, hist_ylim )
        ax.bar(xpos, data[row,:])
        for x in xpos:
            ax.text( x - 0.325, data[row,x-1] + 0.005, "%.3f" % data[row,x-1] )

    update_hist(0, rankmatrix)

    anim = animation.FuncAnimation(fig, update_hist, rankmatrix.shape[0], fargs=(rankmatrix, ) )

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=fps)
    anim.save(fname, writer=writer)
