# NetworkResonance
This repository creates the sub-figures that compose Figures 1-8 in  Stark, Levi, Rotstein, 2022, PLOS Computational Biology

## Overview
How one part of the brain responds to periodic input from another part depends on resonant circuit
properties. Resonance is a basic property of physical systems, and has been experimentally observed
at various levels of neuronal organization both in vitro and in vivo. Yet how resonance is generated in
neuronal networks is largely unknown. In particular, whether resonance can be generated directly at 
the level of a network of spiking neurons remains to be determined. Using detailed biophysical
modeling, we develop a conceptual framework according to which resonance at a given level of
organization is generated by the interplay of low- and high-pass filters, implemented at either the
same or across levels of neuronal organization. We tease apart representative, biophysically-plausible
generative mechanisms of resonance at four different levels of organization: membrane potential
fluctuations, single neuron spiking, synaptic transmission, and neuronal networks. We identify
conditions under which resonance at one level can be inherited to another level of organization,
provide conclusive evidence that resonance at each level can be generated without resonance at any
other level, and describe a number of representative routes to network resonance. The proposed
framework facilitates the investigation of resonance in neuronal systems.

## Wrappers for each Figure:
- NR_make_Fig1: 
  - Figure 1. Cycle-averaged firing rate resonance and spike timing resonance.
- NR_make_Fig2: 
  - Figure 2. Resonance generated at the level of membrane potential fluctuations can be inherited to the network level.
- NR_make_Fig3: 
  - Figure 3. Resonance can be generated directly at the spiking level.
- NR_make_Fig4: 
  - Figure 4. Resonance generated at the spiking level can be inherited to the network level.
- NR_make_Fig5: 
  - Figure 5. Resonance generated at the level of postsynaptic potentials can be inherited to the network level.
- NR_make_Fig6: 
  - Figure 6. Intrinsic network resonance can be generated by combining frequency-dependent mechanisms at the level of postsynaptic potentials and at the spiking level.
- NR_make_Fig7: 
  - Figure 7. Inhibition-induced network resonance can be inherited from the level of membrane potential fluctuations.
- NR_make_Fig8: 
  - Figure 8. Inhibition-induced network resonance is sharpened by presynaptic high-pass filtering.

## Analysis
- st_coherence:                     
  - coherence between spike trains and analog signal.
- st_fingerprint:                   
  - firing rate map (frequency-phase) for spike trains and analog signal.
- NR_calc_z_spectral                          
  - impedance and phase in the frequency domain (any input)
- NR_st_coherence_demo:             
  - demonstrate the difference between coherence and firing rate metrics
- NR_sinusoids_to_cmodel:           
  - sinusoid input to single cell/network models
- NR_sinusoids_to_cmodel_run:      
  - wrapper for sinusoids_to_cmodel
- NR_sinusoids_to_cmodel_images:    
  - plot two matrices

## Simulations
- NR_einet                          
  - network of E and I cells: LIF/Inap, all-to-all chemical connectivity
- NR_eisim                          
  - an E-cell and an I-cell that receives input to one of the cells
- NR_syntransmit                    
  - synaptic transmission simulation (LIF/HH model)

## Utilities
- inranges                          
  - determine which elements of a vector are in which range
- myjet                             
  - modified jet with extreme values in pure R,B
- ParseArgPairs                    
  - flexible argument assigning
- plot_raster                       
  - raster display for spike trains
- plotOneSpectrogram               
  - plot a spectrogram/phasogram
- plotTraces                       
  - in a matrix
- sortranges                       
  - to be a set of non-overlapping [ small large ] pairs

 ### To run the code
- Download all routines.
- In MATLAB, write NR_make_Fig#, where # is the number of the figure to be simulated and plotted.
