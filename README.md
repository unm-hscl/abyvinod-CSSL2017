# Open-loop underapproximation of stochastic reach-avoid problem

This directory contains code to recreate the examples from

  Abraham P. Vinod and Meeko M. K. Oishi, "Scalable Underapproximation for the
  Stochastic Reach-Avoid Problem for High-Dimensional LTI Systems Using Fourier
  Transforms", IEEE Control Systems Letters, 2017

# Requirements
    - MATLAB with Optimization toolbox and Global optimization toolbox installed  

# Instructions

For computations that take time, an option to load the matfile or rerun the
computation is provided.

- Run Figure1.m to generate Figure 1 
- Run Figure2.m to generate Figure 2 
- Run scriptForTableI to generate Table 1
- Run scriptForTableII to generate Table 2

## Auxillary instructions

- scriptForChainDI.m
    - Computes FTBU for 20 random points for n in [1,40] and also performs DPBDA
      for n<=3
    - Two options:
        - Rerun the codes: Implement FTBU either as patternsearch and fmincon
          (further compute time reqd) 
        - Load the matfile curseOfDim_full_10_5_20points_n1to40_reqdOnly.mat 
- scriptForComparison.m 
    - DBPDA and FTBU done over a grid with 1681 points
    - Two options:
        - Rerun the codes: Implement FTBU either as patternsearch and fmincon
          (further compute time reqd) 
        - Load the matfile Figure2_0x05_onX_bothPSandFM.mat 
- The parameters for the matfiles --- experiment are S=[-10,10]^n, T=[-5,5]^n, 20 points per
  dimension, n\in[1,40], and has the elapsed times and number of testing points
  only

## Contact details

* Abraham P. Vinod ([aby.vinod@gmail.com](mailto:aby.vinod@gmail.com))
* Meeko M. K. Oishi ([oishi@unm.edu](mailto:oishi@unm.edu))
