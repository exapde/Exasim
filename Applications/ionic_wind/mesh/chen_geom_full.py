import gmsh
import sys
import numpy as np
from gmsh import model as gm

# THIS GENERATES THE COARSE, 435 ELEMENT MESH USED FOR THE INITIAL ROUND OF SIMULATIONS

#############################################
volume_gen=True
# volume_gen=False

case_name = 'chen_geom'
#############################################

gmsh.initialize()
gm.add(case_name)

############## CREATE PLANE GEOMETRY ##############

gm.occ.importShapes('chen_geom.STEP')
# gm.occ.synchronize()
# gmsh.fltk.run()
# exit()

# print(gm.occ.getEntities())
# exit()
# Exraction of the bbox parameters
xmin, ymin, zmin, xmax, ymax, zmax = gm.occ.getBoundingBox(2, 1)
size_scale_factor = .1/(xmax-xmin)

# Scale and rotate mesh -> This is to align the imported CAD geometry with the desired CSYS and scaling for the problem.
dimtags = gm.occ.getEntities()
gmsh.model.occ.dilate(dimtags, 0,0,0, size_scale_factor, size_scale_factor, size_scale_factor)

gmsh.model.occ.rotate(dimtags, 0,0,0, 1, 0, 0, np.pi/2)
xmin, ymin, zmin, xmax, ymax, zmax = gm.occ.getBoundingBox(2, 1)    # Need to over-write the BBox values so that we can generate an accurate length scale later

gm.occ.synchronize()
# print(gm.getEntities())
# gmsh.fltk.run()
# exit()

xlength = xmax-xmin
ylength = ymax-ymin
zlength = zmax-zmin

ref_len = max(xlength, ylength, zlength)

xoffset = 2*ref_len
yoffset = 2*ref_len
zoffset = 2*ref_len

x0 = xmin - xoffset
y0 = ymin - yoffset
z0 = zmin - zoffset

dx = xlength + 2*xoffset
dy = ylength + 2*yoffset
dz = zlength + 2*zoffset

# Outer boundary of box
hmax=dx/50

# Needle tip
gm.mesh.field.add("Distance", 4)
gm.mesh.field.setNumbers(4, "CurvesList", [6])
gm.mesh.field.setNumbers(4, "PointsList", [7])
gmsh.model.mesh.field.add("MathEval", 13)
gmsh.model.mesh.field.setString(13, "F", "0.047619*F4 + 9.523e-7")

# Needle body
gm.mesh.field.add("Distance", 3)
gm.mesh.field.setNumbers(3, "CurvesList", [2, 3, 4])
gmsh.model.mesh.field.add("MathEval", 12)
gmsh.model.mesh.field.setString(12, "F", "5*(F3/{})^2 + .0005".format(ref_len))

# Outflow and grounded cylinder
gm.mesh.field.add("Distance", 2)
gm.mesh.field.setNumbers(2, "CurvesList", [8, 9, 10, 11])
gmsh.model.mesh.field.add("MathEval", 11)
gmsh.model.mesh.field.setString(11, "F", "5*(F2/{})^2 + .002".format(ref_len))

# Axisymmetry boundary
gm.mesh.field.add("Distance", 1)
gm.mesh.field.setNumbers(1, "CurvesList", [7])
gmsh.model.mesh.field.add("MathEval", 10)
gmsh.model.mesh.field.setString(10, "F", "5*(F1/{})^2 + .002".format(ref_len))

# Let's use the minimum of all the fields as the background mesh field:
gm.mesh.field.add("Min", 20)
gm.mesh.field.setNumbers(20, "FieldsList", [12, 13])
# gm.mesh.field.setNumbers(20, "FieldsList", [12])

gmsh.model.addPhysicalGroup(1, [2,3,4,5,6], name='needle')
gmsh.model.addPhysicalGroup(1, [7], name='symmetry_axis')
gmsh.model.addPhysicalGroup(1, [8], name='outflow')
gmsh.model.addPhysicalGroup(1, [9,10,11,12], name='ground_and_cylinder')
gmsh.model.addPhysicalGroup(1, [1, 13], name='farfield')
gmsh.model.addPhysicalGroup(2, [1], name='omega')

gm.mesh.field.setAsBackgroundMesh(20)
gmsh.option.setNumber("Mesh.MeshSizeMax", .001)
gm.mesh.setSmoothing(3, 3, 10)

gm.mesh.generate(2)
gmsh.write('chen_geom_coarse3.msh')
gmsh.fltk.run()
exit()

"mahogany row, that's mit central"
if volume_gen:
    gm.mesh.generate(3)
else:
    gm.mesh.generate(2)

if volume_gen:
    if assign_pg:
        gmsh.write(case_name+'PG.msh3')
        gmsh.write(case_name+'PG.su2')

        # Copy physical groups from one SU2 file to the other - workaround due to gmsh bug - HAS TO BE DONE AFTER THE MESH W/O PHYSICAL GROUPS HAS BEEN GENERATED
        with open(case_name+'PG.su2', 'r') as file:
            lines = file.readlines()

        marker_start = 'NMARK= 7\n'

        idx = lines.index(marker_start)

        pg = lines[idx:]

        with open(case_name+'.su2', 'a') as file:
            file.writelines(pg)
    else:
        gmsh.write(case_name+'.msh3')
        gmsh.write(case_name+'.su2')

else:
    gmsh.write(case_name+'_surface_processed.msh3')

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
