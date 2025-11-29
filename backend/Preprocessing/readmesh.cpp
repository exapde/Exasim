
/*
 * Mesh File Readers and Utilities
 * ==============================
 * 
 * This file provides functions and data structures for reading mesh data from various file formats
 * into a unified Mesh structure. Supported formats include binary, text, Gmsh (.msh v2/v4), VTK (.vtk),
 * and VTU (.vtu). The Mesh structure stores mesh topology, geometry, and auxiliary data in flattened
 * column-major arrays for efficient processing.
 * 
 * Main Components:
 * ----------------
 * 
 * struct Mesh
 *   - Stores mesh points, connectivity, dimensions, element types, and auxiliary arrays for advanced features.
 * 
 * Mesh Readers:
 *   - void readMeshFromBinaryFile(const std::string& filename, Mesh& mesh)
 *       Reads mesh data from a custom binary format.
 *   - void readMeshFromTextFile(const std::string& filename, Mesh& mesh)
 *       Reads mesh data from a custom text format.
 *   - void readMeshFromMshV2File(const std::string& filename, Mesh& mesh)
 *       Reads mesh data from Gmsh v2.x ASCII files.
 *   - void readMeshFromMshV4File(const std::string& filename, Mesh& mesh)
 *       Reads mesh data from Gmsh v4.x ASCII files.
 *   - void readMeshFromGmshFile(const std::string& filename, Mesh& mesh)
 *       Dispatches to the appropriate Gmsh reader based on file version.
 *   - void readMeshFromVTKFile(const std::string& filename, Mesh& mesh)
 *       Reads mesh data from legacy VTK ASCII files.
 *   - void readMeshFromVTUFile(const std::string& filename, Mesh& mesh)
 *       Reads mesh data from VTU XML files (ASCII only).
 *   - void readMeshFromFile(const std::string& filename, Mesh& mesh)
 *       General mesh reader that dispatches based on file extension.
 * 
 * Utilities:
 * ----------
 *   - std::string getFileExtension(const std::string& filename)
 *       Extracts and normalizes the file extension.
 *   - void readFieldFromBinaryFile(const std::string& filename, vector<double>& xdg, vector<int>& ndims)
 *       Reads field data from a binary file.
 *   - void readPartitionFromFile(const std::string& filename, vector<int>& elem2cpu, int ne)
 *       Reads mesh partition data from binary or text files.
 * 
 * Notes:
 * ------
 * - All mesh readers assume homogeneous element types within a mesh file.
 * - Mesh points and connectivity are stored in column-major order for compatibility with downstream processing.
 * - Error handling is performed via exceptions or program termination for unsupported formats or malformed files.
 * - Only ASCII mesh files are supported for Gmsh, VTK, and VTU formats.
 * - The Mesh structure is extensible for advanced mesh features (boundary conditions, curved boundaries, etc.).
 */

struct Mesh {
    std::vector<double> p; // flattened nd × np array, column-major
    std::vector<int> t;    // flattened nve × ne array, column-major
    int nd, dim;           // number of spatial dimensions
    int np;                // number of points
    int nve;               // number of vertices per element
    int nvf;               // number of vertices per face
    int nfe;               // number of faces per element
    int ne;                // number of elements    
    int nf;                // number of faces    
    int elemtype;          // elemtype        
    int npe, npf;    
    int nbndexpr, nbcm, nprdexpr, nprdcom;        
    
    // For each local element e, elemGlobalID[e] is its global element index
    std::vector<int> elemGlobalID;
    std::vector<int> nodeGlobalID;  // [np]  global node IDs

    vector<int> f, f2t, t2t, t2f, t2lf, inte, intl, localfaces;    
    std::vector<double> xdg, udg, vdg, wdg, uhat;
    std::vector<int> xdgdims, udgdims, vdgdims, wdgdims, uhatdims, elem2cpu;    
    
    std::vector<int> interfaceConditions;
    std::vector<int> boundaryConditions;
    std::vector<int> curvedBoundaries;
    std::vector<int> periodicBoundaries1;
    std::vector<int> periodicBoundaries2;
    std::vector<int> cartGridPart;
    
    char** boundaryExprs = nullptr;
    char** curvedBoundaryExprs = nullptr;
    char** periodicExprs1 = nullptr;
    char** periodicExprs2 = nullptr;        
};

void readMeshFromBinaryFile(const std::string& filename, Mesh& mesh) 
{
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (!in) error("Unable to open file " + filename);

    int* ndims = readiarrayfromdouble(in, 4);
    mesh.nd = ndims[0];
    mesh.np = ndims[1];
    mesh.nve = ndims[2];
    mesh.ne = ndims[3];
    free(ndims);
    
    readarray(in, mesh.p, mesh.np * mesh.nd);    
    readiarrayfromdouble(in, mesh.t, mesh.ne * mesh.nve);
            
    in.close();   
}

void readMeshFromTextFile(const std::string& filename, Mesh& mesh)
{    
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    // Read mesh dimensions
    infile >> mesh.nd >> mesh.np >> mesh.nve >> mesh.ne;
    if (mesh.nd <= 0 || mesh.np <= 0 || mesh.nve <= 0 || mesh.ne <= 0) {
        throw std::runtime_error("Invalid mesh dimensions.");
    }

    mesh.p.resize(mesh.nd * mesh.np);
    mesh.t.resize(mesh.nve * mesh.ne);

    // Read points (column-major)
    for (int j = 0; j < mesh.np; ++j) {
        for (int i = 0; i < mesh.nd; ++i) {
            if (!(infile >> mesh.p[i + mesh.nd * j])) {
                throw std::runtime_error("Error reading point at (" + std::to_string(i) + ", " + std::to_string(j) + ")");
            }
        }
    }

    // Read connectivity (column-major)
    for (int j = 0; j < mesh.ne; ++j) {
        for (int i = 0; i < mesh.nve; ++i) {
            if (!(infile >> mesh.t[i + mesh.nve * j])) {
                throw std::runtime_error("Error reading element at (" + std::to_string(i) + ", " + std::to_string(j) + ")");
            }
        }
    }    
    
    infile.close();   
}

void readMeshFromMshV2File(const std::string& filename, Mesh& mesh) 
{
    std::ifstream infile(filename);
    if (!infile.is_open()) error("Error opening file: " + filename);

    // Check for MeshFormat and version
    std::string line;
    bool foundMeshFormat = false;
    while (std::getline(infile, line)) {
        if (line == "$MeshFormat") {
            foundMeshFormat = true;
            break;
        }
    }

    if (!foundMeshFormat) error("Error: $MeshFormat section not found.");    

    double version;
    int fileType, dataSize;
    infile >> version >> fileType >> dataSize;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip to end of line

    if (version < 2.0 || version >= 3.0) 
      error("Error: Unsupported Gmsh version. Only version 2.x is supported.");        

    // Skip $EndMeshFormat
    while (std::getline(infile, line)) {
        if (line == "$Nodes") break;
    }

    int numNodes;
    infile >> numNodes;
    std::vector<std::array<double, 3>> points(numNodes);

    int dim = 1;

    for (int i = 0; i < numNodes; ++i) {
        int index;
        double x, y, z;
        infile >> index >> x >> y >> z;

        // Infer dimension conservatively
        if (std::abs(z) > 1e-14) dim = std::max(dim, 3);
        else if (std::abs(y) > 1e-14) dim = std::max(dim, 2);

        points[i] = {x, y, z};
    }

    // Skip to $Elements
    while (std::getline(infile, line)) {
        if (line == "$Elements") break;
    }

    int numElements;
    infile >> numElements;

    std::vector<std::vector<int>> elements;
    int expectedElementType = -1;
    int numVerticesPerElement = -1;

    const std::unordered_map<int, int> elementTypeToNVE = {
        {1, 2}, // line
        {2, 3}, // triangle
        {3, 4}, // quad
        {4, 4}, // tetra
        {5, 8}  // hex
    };

    for (int i = 0; i < numElements; ++i) {
        int index, etype, numTags;
        infile >> index >> etype >> numTags;

        for (int j = 0; j < numTags; ++j) {
            int dummy;
            infile >> dummy;
        }

        if (elementTypeToNVE.find(etype) == elementTypeToNVE.end()) 
        {
            std::cerr << "Error in readMeshFromMshV2File: Unsupported element type " << etype << "." << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (expectedElementType == -1) {
            expectedElementType = etype;
            numVerticesPerElement = elementTypeToNVE.at(etype);
        } else if (etype != expectedElementType) {
            std::cerr << "Error in readMeshFromMshV2File: Mesh contains multiple element types (" 
                      << expectedElementType << " and " << etype << ")." << std::endl;
            std::exit(EXIT_FAILURE);
        }

        std::vector<int> conn(numVerticesPerElement);
        for (int j = 0; j < numVerticesPerElement; ++j) {
            infile >> conn[j];
            conn[j] -= 1; // 0-based index
        }
        elements.push_back(conn);
    }

    mesh.nd = dim;
    mesh.np = numNodes;
    mesh.ne = static_cast<int>(elements.size());
    mesh.nve = numVerticesPerElement;

    // Flatten p (column-major)
    mesh.p.resize(mesh.nd * mesh.np);
    for (int j = 0; j < mesh.nd; ++j)
        for (int i = 0; i < mesh.np; ++i)
            mesh.p[j + i * mesh.nd] = points[i][j];

    // Flatten t (column-major)
    mesh.t.resize(mesh.nve * mesh.ne);
    for (int j = 0; j < mesh.nve; ++j)
        for (int i = 0; i < mesh.ne; ++i)
            mesh.t[j + i * mesh.nve] = elements[i][j];

    infile.close();   
}

void readMeshFromMshV4File(const std::string& filename, Mesh& mesh) 
{
    std::ifstream infile(filename);
    if (!infile.is_open()) error("Error opening file: " + filename);

    std::string line;
    double version;
    bool foundMeshFormat = false;

    while (std::getline(infile, line)) {
        if (line == "$MeshFormat") {
            infile >> version;
            int fileType, dataSize;
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            foundMeshFormat = true;
            break;
        }
    }

    if (!foundMeshFormat || version < 4.0 || version >= 5.0) {
        std::cerr << "Error in readMeshFromMshV4File: Unsupported Gmsh version " << version << ". Only version 4.x is supported." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    while (std::getline(infile, line)) {
        if (line == "$Nodes") break;
    }

    size_t numEntityBlocks, totalNumNodes, minNodeTag, maxNodeTag;
    infile >> numEntityBlocks >> totalNumNodes >> minNodeTag >> maxNodeTag;

    std::vector<std::array<double, 3>> points(totalNumNodes);
    int dim = 1;
    size_t counter = 0;

    for (size_t b = 0; b < numEntityBlocks; ++b) {
        int entityDim, entityTag, parametric;
        size_t numNodesInBlock;
        infile >> entityDim >> entityTag >> parametric >> numNodesInBlock;

        std::vector<size_t> tags(numNodesInBlock);
        for (size_t i = 0; i < numNodesInBlock; ++i) {
            infile >> tags[i];
        }

        for (size_t i = 0; i < numNodesInBlock; ++i) {
            double x, y, z;
            infile >> x >> y >> z;
            points[counter] = {x, y, z};
            if (std::abs(z) > 1e-14) dim = std::max(dim, 3);
            else if (std::abs(y) > 1e-14) dim = std::max(dim, 2);
            ++counter;
        }
    }

    while (std::getline(infile, line)) {
        if (line == "$Elements") break;
    }

    size_t numElemBlocks, totalNumElements, minElemTag, maxElemTag;
    infile >> numElemBlocks >> totalNumElements >> minElemTag >> maxElemTag;

    std::vector<std::vector<int>> elements;
    int expectedElementType = -1;
    int numVerticesPerElement = -1;

    const std::unordered_map<int, int> elementTypeToNVE = {
        {1, 2}, {2, 3}, {3, 4}, {4, 4}, {5, 8}
    };

    for (size_t b = 0; b < numElemBlocks; ++b) {
        int entityDim, entityTag, etype;
        size_t numElementsInBlock;
        infile >> entityDim >> entityTag >> etype >> numElementsInBlock;

        if (elementTypeToNVE.find(etype) == elementTypeToNVE.end()) {
            std::cerr << "Error in readMeshFromMshV4File: Unsupported element type: " << etype << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (expectedElementType == -1) {
            expectedElementType = etype;
            numVerticesPerElement = elementTypeToNVE.at(etype);
        } else if (etype != expectedElementType) {
            std::cerr << "Error in readMeshFromMshV4File: Mixed element types not supported." << std::endl;
            std::exit(EXIT_FAILURE);
        }

        for (size_t i = 0; i < numElementsInBlock; ++i) {
            int elemTag;
            infile >> elemTag;
            std::vector<int> conn(numVerticesPerElement);
            for (int j = 0; j < numVerticesPerElement; ++j) {
                infile >> conn[j];
                conn[j] -= 1;
            }
            elements.push_back(conn);
        }
    }
    
    mesh.nd = dim;
    mesh.np = totalNumNodes;
    mesh.ne = static_cast<int>(elements.size());
    mesh.nve = numVerticesPerElement;

    mesh.p.resize(mesh.nd * mesh.np);
    for (int j = 0; j < mesh.nd; ++j)
        for (int i = 0; i < mesh.np; ++i)
            mesh.p[j + i * mesh.nd] = points[i][j];

    mesh.t.resize(mesh.nve * mesh.ne);
    for (int j = 0; j < mesh.nve; ++j)
        for (int i = 0; i < mesh.ne; ++i)
            mesh.t[j + i * mesh.nve] = elements[i][j];    
    
    infile.close();   
}

void readMeshFromGmshFile(const std::string& filename, Mesh& mesh) 
{
    std::ifstream infile(filename);
    if (!infile.is_open()) error("Error opening file: " + filename);

    std::string line;
    while (std::getline(infile, line)) {
        if (line == "$MeshFormat") {
            break;
        }
    }

    if (!infile) {
        std::cerr << "Error in readMeshFromGmshFile: $MeshFormat not found in file: " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::getline(infile, line); // read version line
    std::istringstream versionStream(line);
    double version;
    int fileType;
    versionStream >> version >> fileType;

    if (fileType != 0) {
        std::cerr << "Error in readMeshFromGmshFile: Only ASCII .msh files are supported (binary detected)" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (version >= 4.0) {
        std::cout << "Detected Gmsh format v" << version << " → using V4 reader\n";
        readMeshFromMshV4File(filename, mesh);
    } else if (version >= 2.0 && version < 4.0) {
        std::cout << "Detected Gmsh format v" << version << " → using V2/V3 reader\n";
        readMeshFromMshV2File(filename, mesh);
    } else {
        std::cerr << "Error in readMeshFromGmshFile: Unsupported Gmsh mesh format version " << version << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void readMeshFromVTKFile(const std::string& filename, Mesh& mesh) 
{
    std::ifstream infile(filename);
    if (!infile.is_open()) error("Error opening file: " + filename);

    std::string line;
    while (std::getline(infile, line)) {
        if (line.find("POINTS") == 0) break;
    }

    std::istringstream header(line);
    std::string token;
    int numPoints;
    header >> token >> numPoints >> token; // skip "POINTS <np> float"

    std::vector<double> coordinates;
    double x, y, z;
    while (coordinates.size() < 3 * numPoints && infile >> x >> y >> z) {
        coordinates.push_back(x);
        coordinates.push_back(y);
        coordinates.push_back(z);
    }

    int dim = 3;
    for (int i = 0; i < numPoints; ++i) {
        if (std::abs(coordinates[3*i+2]) > 1e-14) {
            dim = 3;
            break;
        } else if (std::abs(coordinates[3*i+1]) > 1e-14) {
            dim = std::max(dim, 2);
        } else {
            dim = std::max(dim, 1);
        }
    }

    while (std::getline(infile, line)) {
        if (line.find("CELLS") == 0) break;
    }

    std::istringstream cellHeader(line);
    int numCells, totalIndices;
    cellHeader >> token >> numCells >> totalIndices;

    std::vector<std::vector<int>> elements;
    for (int i = 0; i < numCells; ++i) {
        int nverts;
        infile >> nverts;
        std::vector<int> conn(nverts);
        for (int j = 0; j < nverts; ++j) {
            infile >> conn[j];
        }
        elements.push_back(conn);
    }

    int nve = elements[0].size();
    for (const auto& elem : elements) {
        if ((int)elem.size() != nve) {
            std::cerr << "Error in readMeshFromVTKFile: Mixed element types not supported in VTK reader.\n";
            std::exit(EXIT_FAILURE);
        }
    }

    mesh.nd = dim;
    mesh.np = numPoints;
    mesh.ne = elements.size();
    mesh.nve = nve;

    // Store points in column-major nd × np format
    mesh.p.resize(mesh.nd * mesh.np);
    for (int j = 0; j < mesh.nd; ++j) {
        for (int i = 0; i < mesh.np; ++i) {
            mesh.p[j + i * mesh.nd] = coordinates[3 * i + j];
        }
    }

    // Store connectivity in column-major nve × ne format
    mesh.t.resize(mesh.nve * mesh.ne);
    for (int j = 0; j < mesh.nve; ++j) {
        for (int i = 0; i < mesh.ne; ++i) {
            mesh.t[j + i * mesh.nve] = elements[i][j];
        }
    }    
    
    infile.close();   
}

void readMeshFromVTUFile(const std::string& filename, Mesh& mesh) 
{
    std::ifstream infile(filename);
    if (!infile.is_open()) error("Error opening file: " + filename);
    
    std::string line;
    bool readingPoints = false, readingConnectivity = false, readingOffsets = false, readingTypes = false;
    std::vector<std::array<double, 3>> points;
    std::vector<int> connectivity;
    std::vector<int> offsets;
    std::vector<int> types;

    while (std::getline(infile, line)) {
        std::istringstream iss(line);

        if (line.find("<DataArray") != std::string::npos && line.find("Points") != std::string::npos) {
            readingPoints = true;
            continue;
        }
        if (line.find("<DataArray") != std::string::npos && line.find("connectivity") != std::string::npos) {
            readingConnectivity = true;
            continue;
        }
        if (line.find("<DataArray") != std::string::npos && line.find("offsets") != std::string::npos) {
            readingOffsets = true;
            continue;
        }
        if (line.find("<DataArray") != std::string::npos && line.find("types") != std::string::npos) {
            readingTypes = true;
            continue;
        }

        if (line.find("</DataArray>") != std::string::npos) {
            readingPoints = readingConnectivity = readingOffsets = readingTypes = false;
            continue;
        }

        if (readingPoints) {
            double x, y, z;
            while (iss >> x >> y >> z) {
                points.push_back({x, y, z});
            }
        }

        if (readingConnectivity) {
            int val;
            while (iss >> val) {
                connectivity.push_back(val);
            }
        }

        if (readingOffsets) {
            int val;
            while (iss >> val) {
                offsets.push_back(val);
            }
        }

        if (readingTypes) {
            int val;
            while (iss >> val) {
                types.push_back(val);
            }
        }
    }

    // Determine dimension
    mesh.np = points.size();
    mesh.nd = 3;
    for (int d = 2; d >= 1; --d) {
        bool nonzero = std::any_of(points.begin(), points.end(), [d](const auto& pt) { return pt[d] != 0.0; });
        if (nonzero) {
            mesh.nd = d + 1;
            break;
        }
    }

    // Validate homogeneous cell types
    if (!types.empty() && !std::all_of(types.begin(), types.end(), [&](int t) { return t == types[0]; })) {
        std::cerr << "Error in readMeshFromVTUFile: VTU file contains mixed cell types, which is not supported.\n";
        std::exit(EXIT_FAILURE);
    }

    // Determine number of vertices per element
    mesh.nve = offsets[0];
    for (size_t i = 1; i < offsets.size(); ++i)
        mesh.nve = std::min(mesh.nve, offsets[i] - offsets[i - 1]);

    // Assemble connectivity matrix
    mesh.ne = offsets.size();
    mesh.t.resize(mesh.nve * mesh.ne);
    for (int e = 0; e < mesh.ne; ++e) {
        int start = (e == 0) ? 0 : offsets[e - 1];
        for (int j = 0; j < mesh.nve; ++j) {
            mesh.t[j + e * mesh.nve] = connectivity[start + j];
        }
    }

    // Assemble point coordinates in column-major order
    mesh.p.resize(mesh.nd * mesh.np);
    for (int j = 0; j < mesh.nd; ++j)
        for (int i = 0; i < mesh.np; ++i)
            mesh.p[j + i * mesh.nd] = points[i][j];    
    
    infile.close();   
}

std::string getFileExtension(const std::string& filename) 
{
    size_t dot = filename.find_last_of('.');
    if (dot == std::string::npos)
        return "";
    std::string ext = filename.substr(dot + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext;
}

void readMeshFromFile(const std::string& filename, Mesh& mesh) 
{
    std::string ext = getFileExtension(filename);

    if (ext == "bin") {
        readMeshFromBinaryFile(filename, mesh);
    } else if (ext == "txt") {
        readMeshFromTextFile(filename, mesh);
    } else if (ext == "msh") {
        readMeshFromGmshFile(filename, mesh);
    } else if (ext == "vtk") {
        readMeshFromVTKFile(filename, mesh);
    } else if (ext == "vtu") {
        readMeshFromVTUFile(filename, mesh);
    } else {
        error("Unsupported mesh file format: " + ext);        
    }
        
    for (int i=0; i<mesh.ne * mesh.nve; i++) mesh.t[i] -= 1;      
    
    // determine element type
    if (mesh.nd == 2 && mesh.nve == 4) mesh.elemtype = 1;
    else if (mesh.nd == 3 && mesh.nve == 8) mesh.elemtype = 1;
    else mesh.elemtype = 0;

    mesh.nvf = (mesh.nd == 3) ? (mesh.nd + mesh.elemtype) : mesh.nd;
    mesh.nfe = mesh.nd + (mesh.nd - 1) * mesh.elemtype + 1;                    
}

void readFieldFromBinaryFile(const std::string& filename, vector<double>& xdg, vector<int>& ndims) 
{
    std::string ext = getFileExtension(filename);
    
    if (ext == "bin") {
        std::ifstream in(filename, std::ios::in | std::ios::binary);
        if (!in) error("Unable to open file " + filename);

        readiarrayfromdouble(in, ndims, 3);
        readarray(in, xdg, ndims[0] * ndims[1] * ndims[2]);        
        in.close();   
    } else {
        error("Unsupported mesh file format: " + ext);      
    }
}

void readPartitionFromFile(const std::string& filename, vector<int>& elem2cpu, int ne) 
{
    std::string ext = getFileExtension(filename);
    
    if (ext == "bin") {
        std::ifstream in(filename, std::ios::in | std::ios::binary);
        if (!in) error("Unable to open file " + filename);
        readiarrayfromdouble(in, elem2cpu, ne);        
        in.close();              
    } else if (ext == "txt") {
        std::ifstream infile(filename);
        if (!infile) error("Unable to open file " + filename);
        elem2cpu.resize(ne);
        for (int j = 0; j < ne; ++j) 
            if (!(infile >> elem2cpu[j])) {
                throw std::runtime_error("Error reading point at " + std::to_string(j));
            }
    } else {
        error("Unsupported mesh partition file format: " + ext);      
    }
}


#ifdef HAVE_MPI     
// -----------------------------------------------------------------------------
// Parallel reader: each rank reads *its portion* of p and t from the same file
// -----------------------------------------------------------------------------
void readParMeshFromBinaryFile(const std::string& filename,
                               Mesh& mesh,
                               MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // -----------------------------------------------------------------
    // 1. Rank 0 reads the global header (nd, np_global, nve, ne_global)
    // -----------------------------------------------------------------
    int nd_global = 0, np_global = 0, nve_global = 0, ne_global = 0;

    if (rank == 0) {
        std::ifstream in(filename, std::ios::in | std::ios::binary);
        if (!in) error("Unable to open file " + filename);

        int* ndims = readiarrayfromdouble(in, 4);
        nd_global  = ndims[0];
        np_global  = ndims[1];
        nve_global = ndims[2];
        ne_global  = ndims[3];
        free(ndims);

        in.close();
    }

    // Broadcast global mesh sizes to all ranks
    int header[4] = { nd_global, np_global, nve_global, ne_global };
    MPI_Bcast(header, 4, MPI_INT, 0, comm);
    nd_global  = header[0];
    np_global  = header[1];
    nve_global = header[2];
    ne_global  = header[3];

    // -----------------------------------------------------------------
    // 2. Compute local np and ne and their offsets in the global arrays
    // -----------------------------------------------------------------
    int np_local   = 0;
    int nodeOffset = 0;  // starting node index (0-based) for this rank

    int ne_local   = 0;
    int elemOffset = 0;  // starting element index (0-based) for this rank

    computeLocalRange(np_global, size, rank, np_local, nodeOffset);
    computeLocalRange(ne_global, size, rank, ne_local, elemOffset);

    // -----------------------------------------------------------------
    // 3. Compute byte offsets into the file
    //
    // File layout (as implied by your serial reader):
    //   - 4 doubles: nd, np, nve, ne
    //   - p: np_global * nd_global doubles
    //   - t: ne_global * nve_global doubles (each storing an integer)
    // -----------------------------------------------------------------
    const std::size_t headerDoubles = 4;
    const std::size_t pDoubles      = static_cast<std::size_t>(np_global) * nd_global;
    const std::size_t tDoubles      = static_cast<std::size_t>(ne_global) * nve_global;

    const std::size_t headerBytes = headerDoubles * sizeof(double);
    const std::size_t pBytes      = pDoubles      * sizeof(double);
    const std::size_t tBytes      = tDoubles      * sizeof(double);

    // Our local slices:
    const std::size_t localPDoubles = static_cast<std::size_t>(np_local) * nd_global;
    const std::size_t localTDoubles = static_cast<std::size_t>(ne_local) * nve_global;

    const std::size_t localPBytes = localPDoubles * sizeof(double);
    const std::size_t localTBytes = localTDoubles * sizeof(double);

    // Start offsets (in doubles) for this rank
    const std::size_t nodeStartDouble = static_cast<std::size_t>(nodeOffset) * nd_global;
    const std::size_t elemStartDouble = static_cast<std::size_t>(elemOffset) * nve_global;

    const std::size_t pStartByte = headerBytes + nodeStartDouble * sizeof(double);
    const std::size_t tStartByte = headerBytes + pBytes + elemStartDouble * sizeof(double);

    // -----------------------------------------------------------------
    // 4. Resize local mesh storage and fill metadata
    // -----------------------------------------------------------------
    mesh.nd  = nd_global;   // same on all ranks
    mesh.nve = nve_global;  // same on all ranks
    mesh.np  = np_local;    // local
    mesh.ne  = ne_local;    // local

    mesh.p.resize(static_cast<std::size_t>(mesh.np) * mesh.nd);
    mesh.t.resize(static_cast<std::size_t>(mesh.ne) * mesh.nve);

    if (np_local == 0 && ne_local == 0) {
        // This rank owns no data; nothing more to do
        return;
    }

    // -----------------------------------------------------------------
    // 5. Each rank opens the file and reads its slices of p and t
    // -----------------------------------------------------------------
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (!in) error("Unable to open file " + filename);

    // 5a. Read local p (coordinates)
    if (np_local > 0) {
        in.seekg(static_cast<std::streamoff>(pStartByte), std::ios::beg);
        in.read(reinterpret_cast<char*>(mesh.p.data()),
                static_cast<std::streamsize>(localPBytes));
        if (!in) error("Error reading local p block from " + filename);
    }

    // 5b. Read local t (connectivity), stored as doubles in file
    if (ne_local > 0) {
        in.seekg(static_cast<std::streamoff>(tStartByte), std::ios::beg);

        std::vector<double> tbuf(localTDoubles);
        in.read(reinterpret_cast<char*>(tbuf.data()),
                static_cast<std::streamsize>(localTBytes));
        if (!in) error("Error reading local t block from " + filename);

        // Convert double -> int, consistent with readiarrayfromdouble
        for (std::size_t i = 0; i < localTDoubles; ++i) {
            mesh.t[i] = static_cast<int>(tbuf[i]);
        }
    }

    in.close();
}

// ----------------------------------------------------------
// Generic parallel mesh reader:
//  - ReaderFunc must be: void(const std::string&, Mesh&)
//  - tagBase: MPI tag base (tagBase for p, tagBase+1 for t)
// ----------------------------------------------------------
template <typename ReaderFunc>
void readParMeshGeneric(const std::string& filename,
                        Mesh& mesh,
                        MPI_Comm comm,
                        ReaderFunc reader,
                        int tagBase)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    Mesh globalMesh;   // only used on rank 0

    int nd_global  = 0;
    int np_global  = 0;
    int nve_global = 0;
    int ne_global  = 0;

    // 1. Rank 0 reads full mesh using provided reader
    if (rank == 0) {
        reader(filename, globalMesh);
        nd_global  = globalMesh.nd;
        np_global  = globalMesh.np;
        nve_global = globalMesh.nve;
        ne_global  = globalMesh.ne;
    }

    // Broadcast global metadata
    int header[4] = { nd_global, np_global, nve_global, ne_global };
    MPI_Bcast(header, 4, MPI_INT, 0, comm);
    nd_global  = header[0];
    np_global  = header[1];
    nve_global = header[2];
    ne_global  = header[3];

    // 2. Compute local np, ne and offsets for this rank
    int np_local   = 0;
    int nodeOffset = 0;
    int ne_local   = 0;
    int elemOffset = 0;

    computeLocalRange(np_global, size, rank, np_local, nodeOffset);
    computeLocalRange(ne_global, size, rank, ne_local, elemOffset);

    // 3. Initialize local mesh metadata and allocate storage
    mesh.nd  = nd_global;
    mesh.nve = nve_global;
    mesh.np  = np_local;
    mesh.ne  = ne_local;

    mesh.p.resize(static_cast<std::size_t>(mesh.nd) * mesh.np);
    mesh.t.resize(static_cast<std::size_t>(mesh.nve) * mesh.ne);

    if (np_local == 0 && ne_local == 0) {
        // This rank owns nothing; done
        return;
    }

    // 4. Distribute nodal coordinates (mesh.p), column-major: p(j + i*nd)
    const int tagP = tagBase;
    if (rank == 0) {
        // Rank 0: fill its local block directly
        for (int i = 0; i < np_local; ++i) {
            int gnode = nodeOffset + i;
            for (int j = 0; j < nd_global; ++j) {
                mesh.p[j + i * nd_global] =
                    globalMesh.p[j + gnode * nd_global];
            }
        }

        // Send other ranks their blocks
        for (int r = 1; r < size; ++r) {
            int np_r, nodeOff_r;
            computeLocalRange(np_global, size, r, np_r, nodeOff_r);
            if (np_r == 0) continue;

            std::vector<double> sendP(static_cast<std::size_t>(np_r) * nd_global);
            for (int i = 0; i < np_r; ++i) {
                int gnode = nodeOff_r + i;
                for (int j = 0; j < nd_global; ++j) {
                    sendP[j + i * nd_global] =
                        globalMesh.p[j + gnode * nd_global];
                }
            }

            MPI_Send(sendP.data(),
                     static_cast<int>(sendP.size()),
                     MPI_DOUBLE, r, tagP, comm);
        }
    } else {
        // Non-root ranks receive their p block
        if (np_local > 0) {
            std::vector<double> recvP(static_cast<std::size_t>(np_local) * nd_global);
            MPI_Recv(recvP.data(),
                     static_cast<int>(recvP.size()),
                     MPI_DOUBLE, 0, tagP, comm, MPI_STATUS_IGNORE);

            for (int i = 0; i < np_local; ++i)
                for (int j = 0; j < nd_global; ++j)
                    mesh.p[j + i * nd_global] = recvP[j + i * nd_global];
        }
    }

    // 5. Distribute connectivity (mesh.t), column-major: t(j + i*nve)
    const int tagT = tagBase + 1;
    if (rank == 0) {
        // Rank 0: fill its local block
        for (int i = 0; i < ne_local; ++i) {
            int gele = elemOffset + i;
            for (int j = 0; j < nve_global; ++j) {
                mesh.t[j + i * nve_global] =
                    globalMesh.t[j + gele * nve_global];
            }
        }

        // Send other ranks their blocks
        for (int r = 1; r < size; ++r) {
            int ne_r, elemOff_r;
            computeLocalRange(ne_global, size, r, ne_r, elemOff_r);
            if (ne_r == 0) continue;

            std::vector<int> sendT(static_cast<std::size_t>(ne_r) * nve_global);
            for (int i = 0; i < ne_r; ++i) {
                int gele = elemOff_r + i;
                for (int j = 0; j < nve_global; ++j) {
                    sendT[j + i * nve_global] =
                        globalMesh.t[j + gele * nve_global];
                }
            }

            MPI_Send(sendT.data(),
                     static_cast<int>(sendT.size()),
                     MPI_INT, r, tagT, comm);
        }
    } else {
        // Non-root ranks receive their t block
        if (ne_local > 0) {
            std::vector<int> recvT(static_cast<std::size_t>(ne_local) * nve_global);
            MPI_Recv(recvT.data(),
                     static_cast<int>(recvT.size()),
                     MPI_INT, 0, tagT, comm, MPI_STATUS_IGNORE);

            for (int i = 0; i < ne_local; ++i)
                for (int j = 0; j < nve_global; ++j)
                    mesh.t[j + i * nve_global] = recvT[j + i * nve_global];
        }
    }
}

// Text
void readParMeshFromTextFile(const std::string& filename,
                            Mesh& mesh,
                            MPI_Comm comm)
{
    // tags 500, 501
    readParMeshGeneric(filename, mesh, comm,
                       readMeshFromTextFile,
                       500);
}

// VTK
void readParMeshFromVTKFile(const std::string& filename,
                            Mesh& mesh,
                            MPI_Comm comm)
{
    // tags 300, 301
    readParMeshGeneric(filename, mesh, comm,
                       readMeshFromVTKFile,
                       300);
}

// VTU
void readParMeshFromVTUFile(const std::string& filename,
                            Mesh& mesh,
                            MPI_Comm comm)
{
    // tags 400, 401
    readParMeshGeneric(filename, mesh, comm,
                       readMeshFromVTUFile,
                       400);
}

// Gmsh v2
void readParMeshFromMshV2File(const std::string& filename,
                              Mesh& mesh,
                              MPI_Comm comm)
{
    // tags 100, 101
    readParMeshGeneric(filename, mesh, comm,
                       readMeshFromMshV2File,
                       100);
}

// Gmsh v4
void readParMeshFromMshV4File(const std::string& filename,
                              Mesh& mesh,
                              MPI_Comm comm)
{
    // tags 200, 201
    readParMeshGeneric(filename, mesh, comm,
                       readMeshFromMshV4File,
                       200);
}

void readParMeshFromGmshFile(const std::string& filename, Mesh& mesh, MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    double version = 0.0;
    int fileType   = 0;
    int foundMeshFormat = 0; // bool as int for MPI

    if (rank == 0) {
        std::ifstream infile(filename);
        if (!infile.is_open()) error("Error opening file: " + filename);

        std::string line;
        while (std::getline(infile, line)) {
            if (line == "$MeshFormat") {
                foundMeshFormat = 1;
                break;
            }
        }

        if (!foundMeshFormat) {
            std::cerr << "Error in readParMeshFromGmshFile: $MeshFormat not found in file: "
                      << filename << std::endl;
            std::exit(EXIT_FAILURE);
        }

        std::getline(infile, line); // read version line
        std::istringstream versionStream(line);
        versionStream >> version >> fileType;

        infile.close();
    }

    // Broadcast detection results to all ranks
    MPI_Bcast(&foundMeshFormat, 1, MPI_INT, 0, comm);
    MPI_Bcast(&version,         1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&fileType,        1, MPI_INT,    0, comm);

    if (!foundMeshFormat) {
        // Rank 0 already printed an error; make sure everyone stops.
        std::cerr << "Error in readParMeshFromGmshFile: $MeshFormat not found (broadcast)\n";
        std::exit(EXIT_FAILURE);
    }

    if (fileType != 0) {
        if (rank == 0) {
            std::cerr << "Error in readParMeshFromGmshFile: Only ASCII .msh files are supported (binary detected)"
                      << std::endl;
        }
        std::exit(EXIT_FAILURE);
    }

    if (version >= 4.0) {
        if (rank == 0)
            std::cout << "Detected Gmsh format v" << version << " → using parallel V4 reader\n";
        readParMeshFromMshV4File(filename, mesh, comm);
    } else if (version >= 2.0 && version < 4.0) {
        if (rank == 0)
            std::cout << "Detected Gmsh format v" << version << " → using parallel V2/V3 reader\n";
        readParMeshFromMshV2File(filename, mesh, comm);
    } else {
        if (rank == 0) {
            std::cerr << "Error in readParMeshFromGmshFile: Unsupported Gmsh mesh format version "
                      << version << std::endl;
        }
        std::exit(EXIT_FAILURE);
    }
}

void readParMeshFromFile(const std::string& filename, Mesh& mesh, MPI_Comm comm) 
{
    std::string ext = getFileExtension(filename);

    if (ext == "bin") {
        readMeshFromBinaryFile(filename, mesh, comm);
    } else if (ext == "txt") {
        readMeshFromTextFile(filename, mesh, comm);
    } else if (ext == "msh") {
        readMeshFromGmshFile(filename, mesh, comm);
    } else if (ext == "vtk") {
        readMeshFromVTKFile(filename, mesh, comm);
    } else if (ext == "vtu") {
        readMeshFromVTUFile(filename, mesh, comm);
    } else {
        error("Unsupported mesh file format: " + ext);        
    }
        
    for (int i=0; i<mesh.ne * mesh.nve; i++) mesh.t[i] -= 1;      
    
    // determine element type
    if (mesh.nd == 2 && mesh.nve == 4) mesh.elemtype = 1;
    else if (mesh.nd == 3 && mesh.nve == 8) mesh.elemtype = 1;
    else mesh.elemtype = 0;

    mesh.nvf = (mesh.nd == 3) ? (mesh.nd + mesh.elemtype) : mesh.nd;
    mesh.nfe = mesh.nd + (mesh.nd - 1) * mesh.elemtype + 1;                    
}

void readLocalMeshFromBinaryFile(const std::string& filename,
                                 Mesh& mesh,
                                 const std::vector<idx_t>& epart_local,
                                 MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // -----------------------------------------------------------------
    // 1. Rank 0 reads global header (nd, np_global, nve, ne_global)
    // -----------------------------------------------------------------
    int nd_global = 0, np_global = 0, nve_global = 0, ne_global = 0;

    if (rank == 0) {
        std::ifstream in(filename, std::ios::in | std::ios::binary);
        if (!in) error("Unable to open file " + filename);

        int* ndims = readiarrayfromdouble(in, 4);
        nd_global  = ndims[0];
        np_global  = ndims[1];
        nve_global = ndims[2];
        ne_global  = ndims[3];
        free(ndims);

        in.close();
    }

    // Broadcast global mesh sizes to all ranks
    int header[4] = { nd_global, np_global, nve_global, ne_global };
    MPI_Bcast(header, 4, MPI_INT, 0, comm);
    nd_global  = header[0];
    np_global  = header[1];
    nve_global = header[2];
    ne_global  = header[3];

    // -----------------------------------------------------------------
    // 2. Reconstruct original block distribution used for ParMETIS
    //    (same computeLocalRange as in readParMeshFromBinaryFile)
    // -----------------------------------------------------------------
    int ne_local_orig = 0;
    int elemOffset    = 0;  // starting global element index for this rank
    computeLocalRange(ne_global, size, rank, ne_local_orig, elemOffset);

    // Sanity check: local ParMETIS output must match original local ne
    if (static_cast<int>(epart_local.size()) != ne_local_orig) {
        std::cerr << "readLocalMeshFromBinaryFile: epart_local.size() = "
                  << epart_local.size()
                  << " but original ne_local = " << ne_local_orig
                  << " on rank " << rank << "\n";
        error("Inconsistent local ParMETIS partition size");
    }

    // -----------------------------------------------------------------
    // 3. Build global element-partition array epart_global[e] (size ne_global)
    // -----------------------------------------------------------------
    std::vector<idx_t> epart_global(ne_global);

    // Determine MPI datatype for idx_t
    MPI_Datatype mpi_idx_t = (sizeof(idx_t) == sizeof(int))
                             ? MPI_INT
                             : MPI_LONG_LONG; // adjust if you use 64-bit idx_t differently

    // Allgatherv parameters
    std::vector<int> recvcounts(size), displs(size);
    for (int r = 0; r < size; ++r) {
        int ne_r = 0, off_r = 0;
        computeLocalRange(ne_global, size, r, ne_r, off_r);
        recvcounts[r] = ne_r;
        displs[r]     = off_r;
    }

    MPI_Allgatherv(epart_local.data(),
                   ne_local_orig,
                   mpi_idx_t,
                   epart_global.data(),
                   recvcounts.data(),
                   displs.data(),
                   mpi_idx_t,
                   comm);

    // -----------------------------------------------------------------
    // 4. Determine which global elements belong to *this* rank after ParMETIS
    // -----------------------------------------------------------------
    std::vector<idx_t> ownedElems; // global element indices for this rank
    ownedElems.reserve(ne_global / size + 1);

    for (idx_t e = 0; e < (idx_t)ne_global; ++e) {
        if (epart_global[e] == rank) {
            ownedElems.push_back(e);
        }
    }

    const int ne_local = static_cast<int>(ownedElems.size());

    // If this rank owns no elements, also owns no nodes
    if (ne_local == 0) {
        mesh.nd  = nd_global;
        mesh.nve = nve_global;
        mesh.ne  = 0;
        mesh.np  = 0;
        mesh.p.clear();
        mesh.t.clear();
        return;
    }

    // -----------------------------------------------------------------
    // 5. Compute file layout constants
    //
    // File layout:
    //   - 4 doubles: nd, np, nve, ne
    //   - p: np_global * nd_global doubles
    //   - t: ne_global * nve_global doubles (each storing an integer)
    // -----------------------------------------------------------------
    const std::size_t headerDoubles = 4;
    const std::size_t pDoubles      = static_cast<std::size_t>(np_global) * nd_global;
    const std::size_t tDoubles      = static_cast<std::size_t>(ne_global) * nve_global;

    const std::size_t headerBytes   = headerDoubles * sizeof(double);
    const std::size_t pBytes        = pDoubles      * sizeof(double);
    (void)tDoubles;

    const std::size_t tBlockStartByte = headerBytes + pBytes;

    // -----------------------------------------------------------------
    // 6. Read connectivity for owned elements and build node mask
    // -----------------------------------------------------------------
    mesh.nd  = nd_global;
    mesh.nve = nve_global;
    mesh.ne  = ne_local;

    // Temporarily store global node IDs
    mesh.t.resize(static_cast<std::size_t>(mesh.nve) * mesh.ne);

    std::vector<char> nodeMask(np_global, 0);

    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (!in) error("Unable to open file " + filename);

    std::vector<double> tmpConn(nve_global);

    for (int eL = 0; eL < ne_local; ++eL) {
        idx_t eGlobal = ownedElems[eL];

        std::size_t elemStartDouble = static_cast<std::size_t>(eGlobal) * nve_global;
        std::size_t elemStartByte   = tBlockStartByte + elemStartDouble * sizeof(double);

        in.seekg(static_cast<std::streamoff>(elemStartByte), std::ios::beg);
        in.read(reinterpret_cast<char*>(tmpConn.data()),
                static_cast<std::streamsize>(nve_global * sizeof(double)));
        if (!in) error("Error reading connectivity for element from " + filename);

        for (int j = 0; j < nve_global; ++j) {
            int gnode = static_cast<int>(tmpConn[j]);
            mesh.t[j + eL * nve_global] = gnode; // global node for now

            if (gnode < 0 || gnode >= np_global) {
                error("readLocalMeshFromBinaryFile: node index out of range");
            }
            nodeMask[gnode] = 1;
        }
    }

    // -----------------------------------------------------------------
    // 7. Build list of local nodes and global→local mapping
    // -----------------------------------------------------------------
    std::vector<int> localNodes;
    localNodes.reserve(np_global);

    std::vector<int> g2l(np_global, -1);

    for (int g = 0; g < np_global; ++g) {
        if (nodeMask[g]) {
            int l = static_cast<int>(localNodes.size());
            localNodes.push_back(g);
            g2l[g] = l;
        }
    }

    const int np_local = static_cast<int>(localNodes.size());
    mesh.np = np_local;

    // -----------------------------------------------------------------
    // 8. Read coordinates for local nodes
    // -----------------------------------------------------------------
    mesh.p.resize(static_cast<std::size_t>(mesh.nd) * mesh.np);

    // Clear EOF flags before re-using the stream
    in.clear();

    std::vector<double> coordBuf(nd_global);

    for (int l = 0; l < np_local; ++l) {
        int gnode = localNodes[l];

        std::size_t nodeStartDouble = static_cast<std::size_t>(gnode) * nd_global;
        std::size_t nodeStartByte   = headerBytes + nodeStartDouble * sizeof(double);

        in.seekg(static_cast<std::streamoff>(nodeStartByte), std::ios::beg);
        in.read(reinterpret_cast<char*>(coordBuf.data()),
                static_cast<std::streamsize>(nd_global * sizeof(double)));
        if (!in) error("Error reading coordinates for node from " + filename);

        for (int j = 0; j < nd_global; ++j) {
            mesh.p[j + l * nd_global] = coordBuf[j];
        }
    }

    in.close();

    // -----------------------------------------------------------------
    // 9. Remap connectivity from global node IDs to local node IDs
    // -----------------------------------------------------------------
    for (int eL = 0; eL < ne_local; ++eL) {
        for (int j = 0; j < nve_global; ++j) {
            int gnode = mesh.t[j + eL * nve_global];
            int lnode = g2l[gnode];
            if (lnode < 0) {
                error("readLocalMeshFromBinaryFile: missing local mapping for node");
            }
            mesh.t[j + eL * nve_global] = lnode;
        }
    }
}

#endif