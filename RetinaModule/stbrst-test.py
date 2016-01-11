import numpy as np
import numpy.random as rnd
import nest
import matplotlib.pyplot as plt

nest.Install("retinamodule")

neuron1 = nest.Create("stbrst_gc_conv")
neuron2 = nest.Create("stbrst_gc_conv")
nest.SetStatus(neuron1,{'totAHPinit':rnd.random(),'seed':rnd.randint(32535) } )
nest.SetStatus(neuron2,{'totAHPinit':rnd.random(),'seed':rnd.randint(32535) } )

multimeter = nest.Create("multimeter",2)
nest.SetStatus(multimeter, {"withtime":True, "record_from":["Vm","synconv"]})
#nest.Connect(multimeter[0], neuron1)
#nest.Connect(multimeter[1], neuron2)

nest.Simulate(1000.0)

#dmm = nest.GetStatus(multimeter)[0]
#Vms = dmm["events"]["Vm"][np.where( dmm["events"]["senders"] == 1 )]
#ts = dmm["events"]["times"][np.where( dmm["events"]["senders"] == 1 )]
#plt.plot(ts,Vms,'b--')

#Vms = dmm["events"]["Vm"][np.where( dmm["events"]["senders"] == 2 )]
#ts = dmm["events"]["times"][np.where( dmm["events"]["senders"] == 2 )]
#plt.plot(ts,Vms,'r--')

#plt.show()


