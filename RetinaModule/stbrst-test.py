import numpy as np
import numpy.random as rnd
import nest
import matplotlib.pyplot as plt

nest.Install("retinamodule")
nest.SetKernelStatus({"resolution": 1.0})

neurons = nest.Create("stbrst_gc_conv",2)
params = [ 
	{
	'totAHPinit':i,
	'seed':s,
	'gnoise':gnoise,
	} for i,s,gnoise in zip(0.2*rnd.random(2),rnd.randint(32535,size=2),[2E-7,0.]) 
]

nest.SetStatus(neurons,params )

res=["Vm","synconv","spont","conCa","fAHP","sAHP"]#DB>>, "ICa", "IAHP","ISyn"]
multimeter = nest.Create("multimeter")
nest.SetStatus(multimeter, {"withtime":True, "record_from":res})
nest.Connect(multimeter,neurons)

#sd = nest.Create("spike_detector")
#nest.Connect(neurons,sd)

#nest.Connect((neurons[0],),(neurons[1],),syn_spec={"weight": 0.24, "delay":2., 'model':"Convolution"})
nest.Connect((neurons[0],),(neurons[1],),syn_spec={"weight": 1., "delay":2., 'model':"Convolution"})


nest.Simulate(2000.0)

dmm = nest.GetStatus(multimeter)[0]

#plt.figure("Voltage and Synaptic Current")

for idx,name in enumerate(res):
	if idx == 0:
		tax=plt.subplot(3,3,idx+1)
	else:
		plt.subplot(3,3,idx+1,sharex=tax)
	Vms = dmm["events"][name][np.where( dmm["events"]["senders"] == 1 )]
	ts = dmm["events"]["times"][np.where( dmm["events"]["senders"] == 1 )]
	plt.ylabel(name)
	plt.plot(ts,Vms,'k-',lw=1.5)


	Vms = dmm["events"][name][np.where( dmm["events"]["senders"] == 2 )]
	ts = dmm["events"]["times"][np.where( dmm["events"]["senders"] == 2 )]
	plt.plot(ts,Vms,'r--',lw=1.5)


plt.show()


