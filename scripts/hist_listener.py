import numpy as np
import pylab as plt
import zmq
import time

port = "5555"
# Socket to talk to server
context = zmq.Context()
socket = context.socket(zmq.SUB)
socket.connect ("tcp://686tupac:%s" % port)
socket.setsockopt(zmq.SUBSCRIBE, "energy_from_filt_value")

lasthistdict = {}
histdict = {}

plt.close("all")
plt.ion()
fig=plt.figure()
ax=plt.gca()
tplotlast = time.time()
total_value = 0
while True:
    try:
        countedname, channum, edges, counts, totalseen, new_timestamps = socket.recv_multipart()
        channum = int(channum)
        start, step, stop = map(float,edges.split(":"))
        bin_edges = np.arange(start, stop+step, step)
        bin_centers = (bin_edges[1:]+bin_edges[:-1])*0.5
        counts = np.fromstring(counts[1:-1],sep=",") # counts has form "[1,2,3]" fromstring wansts "1,2,3"
        totalseen = int(totalseen)
        new_timestamps = np.fromstring(new_timestamps[1:-1], sep=",") # new_timestamps has form "[1,2,3]" fromstring wansts "1,2,3"
        lasthistdict[channum] = histdict.get(channum, (bin_edges, counts, totalseen, new_timestamps))
        histdict[channum] = (bin_edges, counts, totalseen, new_timestamps)

        if (time.time()-5)>tplotlast:
            print("plotting")
            tplotlast = time.time()
            totaltotal = 0
            totalcounts = np.zeros_like(counts)
            numchannels = 0
            timestamps = []
            ax.clear()
            for k,v in histdict.iteritems():
                (bin_edges, counts, totalseen, new_timestamps) = v
                totaltotal+=totalseen
                timestamps.append(new_timestamps[1])
                totalcounts += counts
                numchannels+=1

            ax.plot(bin_centers, totalcounts, drawstyle="steps-mid")
            ax.set_title("mean timestamp %g, included triggers %g\n all triggers %g, num channels %g"%(np.mean(timestamps), np.sum(totalcounts),totaltotal,numchannels))
            ax.set_xlabel("energy (eV)")
            ax.set_ylabel("counts per %0.2f eV bin"%(bin_edges[1]-bin_edges[0]))
            plt.draw()
    except KeyboardInterrupt:
        print("KeyboardInterrupt")
        break
