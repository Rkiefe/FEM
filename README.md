# FEM
Finite Element Method from scratch

This repository is dedicated to creating a simulation of different systems using the Finite Element Method.
It goes from the very start of the process, from the variational form of the differential equation, all the way to the actual programming required to output a result. The main language used is MATLAB but can easily be transposed to others such as Python or C.

The final goal of this project is to simulate the magnetic field of a solid metal with the most general boundary conditions possible in 2D. There are already free software available for this purpose such as FreeFEM or the FEMM, but my objective is to create a program from start to finish including mesh creation.

### Goal
Simulate the magnetic field of solids.

As such, starting from a simpler situation than a dynamic one, results in tackling magnetostatics. The 2 equations that govern such system are:
$$\nabla \cdot \vec{B} = 0$$
And
$$\nabla \times \vec{H} = \vec{j}$$

Which using the Coulomb Gauge $\nabla \cdot \vec{A} = 0$ and the usual vector potential $\vec{A}$, the two equations result in one: $\nabla \times ( \nabla \times \vec{A} ) = \mu \vec{j}$, which can be further simplified into:
$$-\nabla ^2 \vec{A} = \mu \vec{j}$$

Each component of A has its own Poisson equation. So a 1D system is simpler to start constructing the code for, and only after testing for that simpler case, move on to the two dimensional problem.

#### Current Status
Able to solve the *Poisson* equation for the 1D case. This is equivelent to obtaining the displacement of an elastic rod over its length (the bend). The program uses Robin boundary conditions, meaning we take as input the displacement at the ends of the rod and its spring constant in each end, allowing a more general system to be simulated.

> This models a situation where the bar is connected to a spring with spring constant k
Mats G.Larson - Theory Implementation and Applications


