from scan import *

import pyvista as pv
import tetgen

def main(showPlot=False):
    
    # Import object
    filename = "suporte.stl"
    mesh = pv.read(filename) # read the mesh from stl file
    # mesh.triangulate()

    # Determine bounding box of the mesh
    bounds = mesh.bounds

    # Extra space between object and box
    padding = 200  
    box_bounds = [
        bounds[0] - padding, bounds[1] + padding,
        bounds[2] - padding/2, bounds[3] + padding/2,
        bounds[4] - padding/2, bounds[5] + padding/2
    ]

    # Create a cube mesh around the bounding box
    box = pv.Box(box_bounds).triangulate()
    
    combined = mesh + box

    tet = tetgen.TetGen(combined)
    nodes,elements = tet.tetrahedralize() # elem_growth_ratio=0.1
    # nodes, elements, tet = solidMesh(combined)

    print("Number of elements: ",len(elements))

    if showPlot:
        tet.plot(show_edges=True,opacity=0.5) # color='w'

if __name__ == "__main__":
    main(True)
    # main(False)
