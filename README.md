# FEM
Finite Element Method from scratch

This repository is dedicated to creating a simulation of different systems using the Finite Element Method.
It goes from the very start of the process, from the variational form of the differential equation, all the way to the actual programming required to output a result. The main language used is MATLAB but it can be transposed to others such as Python or C.

The final goal of this project is to create a simulation software similar to COMSOL or Ansys, where the user can simulate magnetic fields for materials with non linear magnetic behaviour, heat transfer with temperature dependent conductivity and more, both in 2D and 3D.

### Magnetic fields
Considering magnetostatics, the 2 equations that govern the system are:
$$\nabla \cdot \vec{B} = 0$$
And
$$\nabla \times \vec{H} = \vec{j}$$
While the two are related by
$$\vec{B} = \mu \vec{H}$$

Using the vector potential formulation $\nabla \times \vec{A} = \vec{B}$

the equations result in: $\nabla \times ( 1/\mu \nabla \times \vec{A} ) = \vec{j}$, which can be further simplified by using the Coulomb Gauge $\nabla \cdot \vec{A} = 0$:
$$-\nabla \cdot (1/\mu \nabla \vec{A}) = \vec{j}$$ (set of Poisson equations)

In 2D, $\vec{A} = A(x,y) \vec{e_z}$ makes the gauging used immediatly valid. Each component of A has its own Poisson equation.

### The Poisson Equation and you
The deformation of a bar can be described by the Poisson equation in 1D, so it was used to first try an implementation of FEM for a simple case. There are analytical solutions depending on the external force used. After Achieving good results, the 2D magnetostatics problem mentioned above was up next. Additionally in 2D, the steady-state heat equation is also a poisson equation! 

## Solvign 2D magnetostatic problems
Magnetic materials can be described by their magnetic permeability - how easy is it for the magnetic field to magnetize the object, or you can think of it as how easy it is for the external magnetic field to want to be in the region occupied by the material. When the permeability is linear, we have $$- 1/\mu \nabla^2  A = j$$. When its not, we have the previous form of the equation with the divergent of the gradient. In this repository you have the 2D and 3D magnetostatic non-linear problem solved for finding the magnetic field H when a magnetizable material of random shape is present.
