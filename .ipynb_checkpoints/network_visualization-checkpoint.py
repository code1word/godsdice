import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import networkx as nx

import pageRankSim

# Some example run code for pageRankSim below
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

rankmatrix=np.zeros((n_frames,n_states))

for i in range(1,n_frames+1,1):
    rankrow = pageRankSim.simulate(A,init_vector,int((n_runs/n_frames)*i),shots)
    rankmatrix[i-1,:]=rankrow
"""

def connected_node_mask(A):
    node_mask = [] # true for connected nodes, false for disconnected nodes
    for i in range(A.shape[0]):
        if False not in np.equal( A[i,:], np.zeros(A.shape[0]) ) and False not in np.equal( A[:,i], np.zeros(A.shape[0]) ):
            node_mask.append(False)
        else:
            node_mask.append(True)
    return node_mask

# fname should be .mp4!
def network_viz(A, rankmatrix, fps, fname):
    node_mask = connected_node_mask(A)

    G = nx.from_numpy_matrix( A[node_mask][:,node_mask] )

    fig, ax = plt.subplots()

    pos=nx.circular_layout(G)

    def update_nx(row, data):
        ax.cla()
        #print(pos)
        for k in pos.keys():
            #print(k)
            x = pos[k][0]
            y = pos[k][1]
            circ = plt.Circle( ( x, y ), 0.095, fill=False, edgecolor="black" )
            ax.add_artist(circ)

            ax.text( x + 0.05, y + 0.05, int(k) + 1 )
        ax.axis("equal")
        nx.draw_networkx(G, pos=pos, ax=ax, cmap="gray", node_color=data[row,:], with_labels=False)
        ax.axis('off')

    update_nx(0, rankmatrix) # initialization

    sm = plt.cm.ScalarMappable(cmap="gray", norm=plt.Normalize(vmin=np.min(rankmatrix), vmax=np.max(rankmatrix)))
    plt.colorbar(sm, ax=ax)

    anim = animation.FuncAnimation(fig, update_nx, rankmatrix.shape[0], fargs=(rankmatrix, ))

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=fps)
    anim.save(fname, writer=writer)
