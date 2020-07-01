# %% codecell
# initialize
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
from qiskit import Aer, execute
from qiskit.extensions import UnitaryGate
from qiskit.visualization import plot_histogram
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
import Testing_Unitaries
######################
# MAKE SURE Testing_Unitaries IS UP TO DATE!
######################

# %% codecell

def connected_node_mask(A):
    node_mask = [] # true for connected nodes, false for disconnected nodes
    for i in range(A.shape[0]):
        if False not in np.equal( A[i,:], np.zeros(A.shape[0]) ) and False not in np.equal( A[:,i], np.zeros(A.shape[0]) ):
            node_mask.append(False)
        else:
            node_mask.append(True)
    return node_mask

def beginsim (A, n_frames):
    n_runs = n_frames # number of runs in final iteration
    shots = 50 # number of shots on IBM qasm simulator
    # for graph node 8 disconnected : index = 7
    disc_indices = []

    for i in range(len(A[1,:])):
        if(not connected_node_mask(A)[i]):
            disc_indices.append(i)

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
        rankrow = simulate(A,init_vector,int((n_runs/n_frames)*i),shots,n_states)
        rankmatrix[i-1,:]=rankrow

    return(rankmatrix)

    # %% codecell
    # run circuit

    # initialize to superposition without |up>|111> or |down>|111>
    
def simulate(A, init_vector, n_runs, shots, n_states):

    # number of qubits needed for graph
    n_qubits = int(np.log2(A.shape[0]))

    # make C and S operators
    C = Testing_Unitaries.constructC(A)
    c_gate = UnitaryGate(C)

    S = Testing_Unitaries.constructS(A)
    s_gate = UnitaryGate(S)



    # run the circuits for many timesteps and average the results
    # the paper started at 50 and ran up to 200-500, depending on the size of the
    # graph
    # I'm doing it in steps of 10

    runs = [10*i + 50 for i in range(n_runs)]
    # all the probabilities are added up here to yield the quantum rank

    # creates as many dictionary slots as we need
    rank = {}
    for i in range(n_states):
        rank[bin(i)[2:].zfill(n_qubits)] = 0

    for run in runs:

        # initialize qc
        qc = QuantumCircuit(n_qubits+1,n_qubits)

        qc.initialize(init_vector, [i for i in range(n_qubits+1)])

        # add on appropriate #'s of C and S gates
        for i in range(run):
            qc.append(c_gate, [i for i in range(n_qubits+1)])
            qc.append(s_gate, [i for i in range(n_qubits+1)])

        # measure output probabilities
        # we're 90% sure that the last qubit is the coin
        qc.measure([i for i in range(n_qubits)], [i for i in range(n_qubits)])

        # run circuit and store output probabilities
        simulator = Aer.get_backend("qasm_simulator")
        job = execute(qc, simulator, shots=shots)
        result = job.result()
        counts_run = result.get_counts(qc)
        # add up probabilities for all runs
        for key in counts_run.keys():
            rank[key] += counts_run[key]

        rankdata = list(rank.items())
        rankvect = np.array(rankdata)[:,1]
        rankvectint = [int(i) for i in rankvect]
        rankvectintnorm = [i/np.sum(rankvectint) for i in rankvectint]

    return rankvectintnorm
    # normalize & plot quantum rank
    #plot_histogram(rank)

# %% codecell
############################
# DECOMPOSING UNITARY NOT WORKING
############################
#
# decomp = QuantumCircuit(4,3)
# decomp.append(c_gate.definition, [0,1,2,3])
# decomp.append(s_gate.definition, [0,1,2,3])
# decomp.decompose().draw(output='mpl')

# # %% codecell
########################
# I don't think this block works, but it may be useful later
# Makes a subcircuit out of a UnitaryGate
########################

# # make a circuit with this gate
#
# testckt=unit_gate.definition
# testckt.name = 'test'
# testckt.draw(output='text')
########################
