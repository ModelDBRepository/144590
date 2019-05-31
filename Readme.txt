This is the readme for the model associated with the paper:

Ghanim Ullah*, John R Cressman Jr, Ernest Barreto, and Steven J
Schiff.  (2009) "The Influence of Sodium and Potassium Dynamics on
Excitability, Seizures, and the Stability of Persistent States:
II. Network and Glia Dynamics". J. Computational Neuroscience,
26:171-183.

Abstract:

In these companion papers, we study how the interrelated dynamics of
sodium and potassium affect the excitability of neurons, the
occurrence of seizures, and the stability of persistent states of
activity. We seek to study these dynamics with respect to the
following compartments: neurons, glia, and extracellular space. We are
particularly interested in the slower time-scale dynamics that
determine overall excitability, and set the stage for transient
episodes of persistent oscillations, working memory, or seizures. In
this second of two companion papers, we present an ionic current
network model composed of populations of Hodgkin-Huxley type
excitatory and inhibitory neurons embedded within extracellular space
and glia, in order to investigate the role of micro-environmental
ionic dynamics on the stability of persistent activity. We show that
these networks reproduce seizure-like activity if glial cells fail to
maintain the proper micro-environmental conditions surrounding
neurons, and produce several experimentally testable predictions to
better understand such dynamics. Our work suggests that the stability
of persistent states to perturbation is set by glial activity, and
that how the response to such perturbations decays or grows may be a
critical factor in a variety of disparate transient phenomena such as
working memory, burst firing in neonatal brain or spinal cord, up
states, seizures, and perhaps spreading depression.

----------------------------------------------------------------------

Model files provided by the authors.

Usage:

Programming language: Fortran 90.

Download the program files from ModelDB and compile the main program
file Network.f90 using a Fortran 90 compiler.

Network.f90 couples 100 inhibitory neurons and 100 excitatory neurons
where the membrane potential dynamics of these neurons is taken from
Gutkin et al. model, 2001, J. Computational Neuroscience, 11, 121-134.

The synaptic currents here are modified from that given in Gutkin et
al., 2001 model. The model also includes dynamic potassium and sodium
concentrations that build on the model from companion paper

John R Cressman Jr, Ghanim Ullah, Jokubas Ziburkus, Steven J Schiff,
and Ernest Barreto. (2009) "The Influence of Sodium and Potassium
Dynamics on Excitability, Seizures, and the Stability of Persistent
States: I. Single Neuron Dynamics". J. Computational Neuroscience,
26:159-170.

The results from Network.f90 are stored into data files (see comments
in Network.f90) that include activity and raster plots for two network
types, the membrane potentials, extracellular potassium and
intracellular sodium of excitatory and inhibitory neurons.

The data files containing the activity of the network are read into
another program file called "Activity.f90". Activity.f90 simply
averages the activity of the network over the desired time window (50
or 100ms).

After the simulations for Network.f90 are complete, compile
Activity.f90 with a Fortran 90 compiler and plot the result using data
visualization package of your interest. This will produce graphs
similar to Fig.5 in the paper.
    

*Contact: 
212 Earth Engineering Science Building, 
The Pennsylvania State University, 
University Park, PA, 16802, USA 
Email: ghanim.phy@gmail.com
Voice: (814) 865 6951
Fax: (814) 865 6161
