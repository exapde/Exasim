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

# ref_pt1 = 85
# ref_pt2 = 212
# ref_len = np.linalg.norm(gm.getValue(0, ref_pt1, [0]) - gm.getValue(0, ref_pt2, [0]))

# # Global scale factor - base linear point density
# # base_factor = 100/1849
# base_factor = 100/1849

# # Mesh resolution along curves
# # print('Assigning transfinite curves')
# for curve in gm.occ.getEntities(1):
#     curve_idx = curve[1]

#     length = gm.occ.getMass(1, curve[1])
#     numpts = max(int(length*base_factor), 3)

#     # Custom setting the element densities for important surfaces
#     if curve_idx in [306, 308, 103, 105]:     # engine pylons
#         numpts *= 5

#     if curve_idx in [294, 251, 301, 293, 300, 285, 287, 299]:     # Tail box base
#         numpts *= 2

#     if curve_idx in [1, 256, 253, 255, 43, 45, 44, 47, 12, 17, 252, 6]:     # Fan forward nubs
#         numpts *= 2

#     if curve_idx in [95, 34, 80, 87, 77, 78, 96, 97]:     # Engine sides
#         numpts *= 2

#     if curve_idx in [92, 94, 298]:     # Between engines
#         numpts *= 3

#     if curve_idx in [130, 324, 121, 318, 116, 117, 312, 110, 110, 312]:   # Leading edges
#         numpts *= 4

#     if curve_idx in [302, 99]:  # Engine ducts
#         numpts *= 2

#     if curve_idx in [322, 323, 316, 317, 127, 128, 113, 114]:  # Wingtip chords
#         numpts *= 10

#     if curve_idx in [186, 191, 134, 198, 190]:  # Engine ducts
#         numpts *= 2

#     if curve_idx in [328, 133, 111, 115, 118]:     # Airfoil chords
#         gm.mesh.setTransfiniteCurve(curve_idx, int(numpts*3), coef=1.05)
#         continue

#     if curve_idx in [132, 123, 313, 314, 108, 124, 326]:     # Airfoil chords
#         gm.mesh.setTransfiniteCurve(curve_idx, int(numpts*3), coef=1/1.05)
#         continue

#     gm.mesh.setTransfiniteCurve(curve_idx, numpts)

# aircraft_surfaces = [a[1] for a in gm.getEntities(2)]

# # Assigning farfield boundary box
# ################## Make box ##################

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

gm.mesh.field.add("Distance", 4)
gm.mesh.field.setNumbers(4, "PointsList", [7])

gm.mesh.field.add("Distance", 3)
gm.mesh.field.setNumbers(3, "CurvesList", [3, 4])
gm.mesh.field.setNumbers(3, "PointsList", [5, 6])

gm.mesh.field.add("Distance", 2)
gm.mesh.field.setNumbers(2, "CurvesList", [8, 9, 10, 11])

gm.mesh.field.add("Distance", 1)
gm.mesh.field.setNumbers(1, "CurvesList", [7])

gmsh.model.mesh.field.add("MathEval", 13)
gmsh.model.mesh.field.setString(13, "F", "10*(F4/{})^1.5 + .00015".format(ref_len))

gmsh.model.mesh.field.add("MathEval", 12)
gmsh.model.mesh.field.setString(12, "F", "5*(F3/{})^2 + .001".format(ref_len))

gmsh.model.mesh.field.add("MathEval", 11)
gmsh.model.mesh.field.setString(11, "F", "5*(F2/{})^2 + .002".format(ref_len))

gmsh.model.mesh.field.add("MathEval", 10)
gmsh.model.mesh.field.setString(10, "F", "5*(F1/{})^2 + .002".format(ref_len))

# Let's use the minimum of all the fields as the background mesh field:
gm.mesh.field.add("Min", 20)
gm.mesh.field.setNumbers(20, "FieldsList", [10, 11, 12, 13])

gm.mesh.field.setAsBackgroundMesh(20)
# 
gmsh.option.setNumber("Mesh.MeshSizeMax", hmax)
gm.mesh.setSmoothing(3, 3, 10)

# TODO: convert units in mm by nondimensionalizing and then changing

gm.mesh.generate(2)
gmsh.write('chen_geom_coarse2.msh3')
gmsh.fltk.run()
exit()


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
