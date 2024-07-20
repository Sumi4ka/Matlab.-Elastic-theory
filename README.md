# Name: Elastic theory
## Description:
This code is written for the numerical calculation of problems in the theory of elasticity.
The main equations I use are the stationary equations of equilibrium in displacements. 
Static press conditions or free surface conditions are established at the boundary nodes.
The calculation grid on which calculations are performed is built by the method of thickening a uniform grid.
This project consist of such class as mesh (grid), zone, solver, plt.

    mesh - this class creates and edits (thickens nodes) grid;
    zone - contains information about boundary conditions;
    solver - the main class where all calculations are performed;
    plt - the class that generates illustrations.
    Example - an example where the following Problem
## Problem
an elastic (steel) cylinder is located between non-deformable planes, 
there is no friction between them. 
The entire side surface of the cylinder is compressed by a certain amount.

## Instruction:
Install Matlab 2017a (Matlab 9.2) or newer;
Launch Example.m.
