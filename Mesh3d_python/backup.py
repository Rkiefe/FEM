from scan import *

import pyvista as pv
import tetgen

def main(showPlot=False):
    filename = "suporte.stl"
    mesh = pv.read(filename)
    
    # Tetrahedralize the mesh
    tet = tetgen.TetGen(mesh)
    nodes,elements = tet.tetrahedralize(mindihedral=10, minratio=1.5)
    
    tet_mesh = tet.grid

    if showPlot:
        tet_mesh.plot(show_edges=True, color='w') # tet_mesh


if __name__ == "__main__":
    main(True)
    # main(False)
