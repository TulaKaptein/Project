# Project
Project for Scientific Computing and Programming

How to use the program:

In main.f90 you can change the potential used for the calculation. The options are:
+ "ParticleInBox", which uses the particle in a box model
+ "GaussianPotWell", which uses the Gaussian potential well model

When you run the program, it will ask for an amount of gridpoints and an interval. 
The input for the interval should be like: "lowerbound upperbound". 

The gridpoints are saved in Gridpoints.txt, the results of the three-point scheme are saved in Eigenvalues.txt and Eigenvectors.txt and the results of the shooting algorithm are saved in Results.txt. The values for the Gaussian potential well are saved in GaussianPotWell.txt, if this potential is used. 
