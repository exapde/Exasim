import gmsh
import numpy as np

def read_gmsh(mshfile, dim, pgbool):
    gmsh.initialize()
    gmsh.merge(mshfile)

    nodes=gmsh.model.mesh.getNodes()
    p = np.array(nodes[1]).reshape(-1,3)        # nodes[0] is just a list of the nodes, 1,...,len(nodes)

    elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements()

    # Iterate over element types to find the one that matches the mesh dimension
    elemdims = np.zeros_like(elementTypes)
    for itag, tag in enumerate(elementTypes):
        elemdims[itag] = gmsh.model.mesh.getElementProperties(tag)[1]

    tags_matching_dim = elementTypes[elementTypes==dim]
    if tags_matching_dim.size > 1:
        raise ValueError(f'There is more than one type of element of dimension {dim} in the mshfile, exiting!')
    else:
        elemtag = tags_matching_dim[0]      # This is 1,2,4,15, etc
        elemidx = np.where(elementTypes==elemtag)[0][0]     # In this case, the triangles are located at index 1: elementTypes=[1,2,15]

    numelem = elementTags[elemidx].size      # size -- 1D array
    t = nodeTags[elemidx].reshape(numelem,-1) - 1   # Subtracting 1 to make it 0-indexed

    # Get physical groups
    if pgbool:
        num_physgrp1d = [i[1] for i in gmsh.model.getPhysicalGroups(dim-1)]   # For a 2D mesh, faces are 1D (last argument)
        pg_nodes = {}
        pg_names = {}
        for pg in num_physgrp1d:
            # Note that the actual PG value isn't important for the user, it only matters that it matches up between the PG names (human readable) and the nodes they correspond nto
            pg_nodes[pg] = gmsh.model.mesh.getNodesForPhysicalGroup(dim-1,pg)[0] - 1   # Subtracting 1 to make the nodes 0 indexed
            pg_names[gmsh.model.getPhysicalName(dim-1,pg)] = pg     # Lookup table mapping PG names to their index
    else:
        pg_nodes = None
        pg_names = None

    return p, t, pg_nodes, pg_names

if __name__ == '__main__':
    p,t,pg_nodes,pg_names = read_gmsh('chen_geom_coarse3.msh3', 2, True)
    print(p.shape)
    print(t.shape)
    print(pg_nodes.keys())
    print(pg_names)