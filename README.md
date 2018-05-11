# Open-loop underapproximation of stochastic reach-avoid problem

This directory contains code to recreate the examples from

  Abraham P. Vinod and Meeko M. K. Oishi, "Scalable Underapproximation for the
  Stochastic Reach-Avoid Problem for High-Dimensional LTI Systems Using Fourier
  Transforms", IEEE Control Systems Letters, 2017

# Requirements
    - MATLAB:
        - Statistics and Machine Learning (for qscmvnv.m)
        - Optimization toolbox (for fmincon)
        - Global optimization toolbox installed (for pattersearch) 

# Instructions

For computations that take time, an option to load the matfile or rerun the
computation is provided.

- Run Figure1.m to generate Figure 1 
- Run Figure2.m to generate Figure 2 
- Run scriptForTableI to generate Table 1
- Run scriptForTableII to generate Table 2

## Other main scripts

- scriptForChainDI.m
    - Computes FTBU for 20 random points for n in [1,40] and also performs DPBDA
      for n<=3
- scriptForComparison.m 
    - DBPDA and FTBU done over a grid with 1681 points

## Contact details

* Abraham P. Vinod ([aby.vinod@gmail.com](mailto:aby.vinod@gmail.com))
* Meeko M. K. Oishi ([oishi@unm.edu](mailto:oishi@unm.edu))
