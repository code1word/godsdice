import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from scipy.stats import gaussian_kde

# n is the number of points to evaluate the kde at
def eval_kde(data, n):
    # we use a bandwidth of 1 because we have discrete positions!
    kde = gaussian_kde( data, bw_method=1 )
    pos = np.linspace( np.min(data), np.max(data), n )
    return kde(pos)

# cwalk_dist is a matrix of results from a classical walk, where
# each row corresponds to some number of repetions R of the classical walk experiment with some number of steps
# (each row corresponds to trials w/ a diff. number of steps)

# for right now i'm assuming the qwalk will produce a similar matrix as cwalk
# will adjust as nessecary later
def dist_viz(cwalk_dist, qwalk_dist, fps, fname):
    fig, ax = plt.subplots()

    N = cwalk_dist.shape[0]
    R = cwalk_dist.shape[1]

    cwalk_kde = np.zeros( cwalk_dist.shape )
    # convert cwalk_dist data into kde data
    for i,rowdata in enumerate(cwalk_dist):
        cwalk_kde[i] = eval_kde( rowdata, len(rowdata) )

    qwalk_kde = np.zeros( qwalk_dist.shape )
    # convert cwalk_dist data into kde data
    for i,rowdata in enumerate(qwalk_dist):
        qwalk_kde[i] = eval_kde( rowdata, len(rowdata) )

    min_xval = np.min( np.min(cwalk_dist), np.min(qwalk_dist) )
    max_xval = np.max( np.max(cwalk_dist), np.max(qwalk_dist) )

    dist_xlim_min = min_xval - 1
    dist_xlim_max = max_xval + 1

    # determine whether qwalk or cwalk has the largest maximum val
    max_yval = np.max( np.max(cwalk_kde), np.max(qwalk_kde) )
    dist_ylim = np.min([ np.around(max_yval, decimals=1) + 0.05, 1.0 ])

    def update_dist(row):
        ax.clear()

        # assumes that we go up by 1 step for each row
        ax.text( max_xval - 2, max_yval - 2, row + 1 )

        ax.set_xlim( dist_xlim_min, dist_xlim_max )
        ax.set_ylim( 0, dist_ylim )

        crow = cdata[row]
        qrow = qdata[row]

        cwalk_xpos = np.linspace( np.min( cwalk_dist[row] ), np.max( cwalk_dist[row] ), R )
        qwalk_xpos = np.linspace( np.min( qwalk_dist[row] ), np.max( qwalk_dist[row] ), R )

        ax.plot( cwalk_pos, cwalk_kde[row], label="Classical Walk" )
        ax.plot( qwalk_pos, qwalk_kde[row], label="Quantum Walk" )

        plt.legend()

    update_dist(0)

    anim = animation.FuncAnimation(fig, update_dist, N )

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=fps)
    anim.save(fname, writer=writer)
