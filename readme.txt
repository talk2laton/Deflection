This file explains how to use the Deflection.m code.

This program computes the deflection and slope of a laterally loaded 
statically indeterminate beams. 
To use this program, you call the function placing the arguments in cells
with keywords at the beginning of each cell except for the first 4 arguments.

First Argument
The first argument is the name of the problem as a string e.g.: 'PROB 1'.

Second Argument
The second argument is the flexural rigidity (EI) of the beam. This is the 
product of the modulus of elasticity, and moment of intertial about normal
in/out of page(through the centroid of the cross section).

Third Argument 
This is the argument that describes the supports on the beam. This is a 
matrix holding the location of a support in the first column and the 
support type in the second column along the same row. For example, a fixed
support that restrict only vertical movement of the beam is classed as type 1
while supports that restrict both vertical and rotational movement is classed
 as type 2.

Fourth Argument
This is the vector of points on the beam to used for the computation. This 
nodes divide the beam into N-1 members(N is the number of nodes).

Fifth argument and on
From the fifth argument and onward, we use cells. The first element of the
cell contains a keyword describing what type of load is inside the argument.
The second element is the magnitude of the load while, the third element of
a cell argument is its location.
Keywords: Point Load = 'CF'
Moment = 'M'
Distributed Load = 'DF'
To add a downward point load of magnitude 5N at location 4m, the argument
would be {'CF',-5,4}. Note the negative sign. If the force is acting upward
the argument would be {'CF',5,4};

To add a clockwise moment of magnitude 10N-m at location 14m, the argument
would be {'M',-10,14}. Note the negative sign. If the moment is anticlockwise
the argument would be {'M',10,14};

To add distributed load we need to describe all of them with the minimum
number of point required to describe the profile with the highest
complexity. For example, a linear profile can be described as {'DF',[5,5],[2,10]}
meaning uniform force per unit length of 5N/m from point 2m to 10m. If the
values of the profile were given at 3 points, the code will automatically
assume it to be quadratic. If profile is uniform, the coefficient of the
second and first degrees would be zero.Hence describing the constant 5N/m
from 2m to 10m as {'DF',5,[2,10]}, {'DF',[5,5],[2,10]}, {'DF',[5,5,5],[2,8,10]}
will make no difference. But in case where the values in the force vector
are different, SFBM will generate a polynomial fit for the forces as a
function of position. For instance {'DF',[1,5,5],[2,8,10]} will generate a
quadratic function, while {'DF',[1,4,5],[2,8,10]} will generate a linear
expression and {'DF',[5,5,5],[2,8,10]}will generate a degree zero expression

There is no limit to the number degree of polynomial that can be used. 

For example:
Deflection('Prob 200',3e7,[0,2;4,1],[0,4],{'M',2000,2},{'CF',-2000,3},{'DF',-1000,[0,1]})
Name : Prob 200
Flexural Rigidity EI = 3e7
Supports: location = 0, type = 2
Supports: location = 4, type = 1
Nodes: 0 and 4
Torque(Moment): Constant -2000Nm at 2m
Concentrated Load: Constant -2000N at 3m
Distributed Load: Constant -1000N/m from 0m to 1m

Solution:
v = a vector of deflection at the given nodes
t = a vector of slope at the given nodes