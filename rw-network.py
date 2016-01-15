"""
/**
 *  rw-network.py
 *
 *
 *  This file is test suite for retinal waves model. It uses retinal 
 *  module to model 2D network of Starburst Amacrine cells with
 *  Cholinergic connections modeled by convolution of voltage neighbor cells
 *
 *  Copyright (C) Ruben Tikidji-Hamburyan (rath@gwu.edu)
 *
 */
"""

import sys
import numpy as np
import numpy.random as rnd
import nest
import nest.topology as nt
import matplotlib.pyplot as plt
from matplotlib.animation import *


size = 10,34 # x,y
hstp = 1.	 # grid step (minimal node distance) in microns
gnse = 0.1	 # deviation of node location in hstp
### Create Retinal cell positions.
pos  = [ ]
for y in xrange(size[1]):
	offset = 0.5 if y%2 else 0.
	for x in xrange(size[0]):
		pos.append( [hstp*np.sqrt(3)*(float(x)+offset) + hstp/2, hstp*float(y)/2. + hstp/2]  )
pos = np.array(pos)
lmax = np.max(pos[:,0])+hstp/2, np.max(pos[:,1])+hstp/2
lmax = lmax[0] if lmax[0] > lmax[1] else lmax[1]

### Import  module and set up simulation step.
nest.Install("retinamodule")
nest.SetKernelStatus({"resolution": 1.0})

### Create layer of retinal Star-Burst Amacrin cells and set up parameters
ret_stbrt=nt.CreateLayer({'positions':list(pos),'elements':"stbrst_gc_conv",'extent':[lmax,lmax],'center':[lmax/2., lmax/2]})
nest.SetStatus(nest.GetLeaves(ret_stbrt)[0],
	[ {	'totAHPinit':i,	'seed':s,} for i,s in zip(rnd.random(size[0]*size[1]),rnd.randint(32535,size=size[0]*size[1]) ) ] )

### Connect layer to it self
nt.ConnectLayers(ret_stbrt, ret_stbrt, {
	'connection_type':'divergent',
	'mask': { 'circular': { 'radius': 3.2 } },
	'weights': { 'gaussian': {'p_center':0.24, 'sigma': 1.6} },
})

### Record voltage and synapses
res=["Vm","synconv"]#,"spont","conCa","fAHP","sAHP"]#DB>>, "ICa", "IAHP","ISyn"]
multimeter = nest.Create("multimeter")
nest.SetStatus(multimeter, {"withtime":True, "record_from":res})
nest.Connect(multimeter,nest.GetLeaves(ret_stbrt)[0])


### Simulate 10 sec
nest.Simulate(30000.0)
#nest.Simulate(10000.0)
print "Simulation is done"
### Get Data
Rec = nest.GetStatus(multimeter)[0]["events"]
### Plot traces
for nid in nest.GetLeaves(ret_stbrt)[0]:
	plt.subplot(211)
	V = Rec["Vm"][np.where( Rec["senders"] == nid )]
	t = Rec["times"][np.where( Rec["senders"] == nid )]
	plt.plot(t,V)

	plt.subplot(212)
	V = Rec["synconv"][np.where( Rec["senders"] == nid )]
	t = Rec["times"][np.where( Rec["senders"] == nid )]
	plt.plot(t,V)
plt.show()

ims = []
f3 = plt.figure(3)

for t in Rec["times"][np.where( Rec["senders"] == nest.GetLeaves(ret_stbrt)[0][0] ) ][1000::5]:
	sys.stderr.write("=== {} ===\r".format(t))
	Vnorm = [ (nid,V) for nid,V in zip(Rec["senders"][np.where(Rec["times"] == t)],Rec["Vm"][np.where(Rec["times"] == t)] ) ]
	Vnorm.sort()
	Vnorm = np.array(Vnorm)
	#Vnorm[:,1] = (Vnorm[:,1]+100.)/130.
	#print Vnorm
	#exit(0)
	ims.append((plt.scatter(pos[:,0],pos[:,1],c=Vnorm[:,1],vmin=-90,vmax=20,s=100, cmap=plt.get_cmap('jet')),))

print
print "Making Video"
ims = ArtistAnimation(f3, ims, interval=30, repeat_delay=100, blit=True)
print "Saving Video"
writer = writers['ffmpeg'](fps=15, metadata=dict(artist='Me'), bitrate=1800)
ims.save("RetWaves.mp4", metadata={'artist':'Retinal Waves on NEST', 'copyright':'Ruben Tikidji-Hamburyan'})
print "Show Video"
plt.show()

	
