# 3D tetrahedral mesh
This section is dedicated to creating a 3d mesh from an stl file for open boundary problems.
## What it does
It takes a "sample.stl" file and generates the 3d tetrahedral mesh of that object placed inside a container

## Requirements
It only needs pyvista and tetgen. Both can be installed using pip,
```shell
pip install pyvista && pip install tetgen
```

### Notes
This is not currently in use, it's a work in progress. To become fully functional I need a way to retrieve the surface triangles and label each tetrahedral element according to its original mesh.
For example, an element that belongs to the container should have the label "0" and an element that belongs to the "sample" should be labeled "1".
This way you can import multiple objects and give each different properties such as thermal conductivity and specific heat.
