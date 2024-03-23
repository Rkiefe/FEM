import pyvista as pv
import tetgen

from scan import *


sphere_a = pv.Sphere().triangulate()
sphere_b = pv.Sphere(center=(0.5, 0, 0)).triangulate()

result = sphere_a.boolean_union(sphere_b) #.triangulate()

tet = tetgen.TetGen(result)
tet.tetrahedralize(order=1, mindihedral=10, minratio=1.5)
tet_mesh = tet.grid

tet_mesh.plot(show_edges=True, color='w')
