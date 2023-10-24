# FEM
Finite Element Method from scratch

This repository is dedicated to creating a simulation of different systems using the Finite Element Method.
It goes from the very start of the process, from the variational form of the differential equation, all the way to the actual programming required to output a result. The main language used is MATLAB but it can be transposed to others such as Python or C.

The final goal of this project is to simulate the magnetic field of a solid metal with the most general boundary conditions possible in 2D. There are already free software available for this purpose such as FreeFEM or the FEMM, but my objective is to create a program from start to finish including mesh creation.

### Goal
Simulate the magnetic field of solids.

The 2 equations that govern such system are:
$$\nabla \cdot \vec{B} = 0$$
And
$$\nabla \times \vec{H} = \vec{j}$$
While the two are related by
$$\vec{B} = \mu \vec{H}$$

Using the vector potential formulation $\nabla \times \vec{A} = \vec{B}$

the equations result in: $\nabla \times ( 1/\mu \nabla \times \vec{A} ) = \vec{j}$, which can be further simplified by using the Coulomb Gauge $\nabla \cdot \vec{A} = 0$:
$$-\nabla \cdot (1/\mu \nabla \vec{A}) = \vec{j}$$ (set of Poisson equations)

In 2D, $\vec{A} = A(x,y) \vec{e_z}$ makes the gauging used immediatly valid.

Each component of A has its own Poisson equation. So a 1D system is simpler to start constructing the code for, and only after testing for that simpler case, move on to the two dimensional problem, and eventually solve 3D systems.

The deformation of a bar can be described by the Poisson equation in 1D, so it was used to first try an implementation of FEM for a simple case. There are analytical solutions depending on the external force used. After Achieving good results, the 2D case was tackled: Magnetostatics. When the permeability is linear, we have $$- 1/\mu \nabla^2  A = j$$. When its not, we have the previous form of the equation with the divergent of the gradient.
