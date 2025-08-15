
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
    int ne;                // number of elements    
    int nf;                // number of elements    
    int elemtype;
    int nvf, nfe;
    int npe, npf;    
    int nbndexpr, nbcm, nprdexpr, nprdcom;        
    
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

