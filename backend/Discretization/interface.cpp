#ifndef __INTERFACE
#define __INTERFACE

int interface_elements(vector<int>& inte, vector<int>& intl, const vector<int>& bf, int nfe, int ne, int coupledboundarycondition)
{
    int n=0;
    for (int i=0; i<ne; i++)
      for (int j=0; j<nfe; j++)
        if (bf[j + nfe*i] == coupledboundarycondition) n += 1;

    inte.resize(n);
    intl.resize(n);
    n = 0; 
    for (int i=0; i<ne; i++)
      for (int j=0; j<nfe; j++) 
        if (bf[j + nfe*i] == coupledboundarycondition) {
          inte[n] = i;
          intl[n] = j;
          n += 1;
        }    

    return n;
}

// All indexing is 0-based and arrays are column-major
void fill_xdg(std::vector<double>& xdg,
              const std::vector<double>& dgnodes,
              const std::vector<int>& perm,
              const std::vector<int>& intl,
              const std::vector<int>& inte,
              int npf, int nd, int nfint, int npe)
{
    // Resize xdg to the correct dimensions
    xdg.resize(static_cast<size_t>(npf) * nfint * nd );

    for (int j = 0; j < nfint; ++j) {
        int intl_j = intl[j];
        int inte_j = inte[j];

        for (int i = 0; i < npf; ++i) {
            int row = perm[i + npf * intl_j]; 
            for (int d = 0; d < nd; ++d) {
                //xdg[i + npf * (d + nd * j)] = dgnodes[row + npe * (d + nd * inte_j)];
                xdg[i + npf * (j + nfint * d)] = dgnodes[row + npe * (d + nd * inte_j)];
            }
        }
    }
}

// Match interface faces between xdg1 and xdg2
// xdg1: (npf × nfint1 × nd)
// xdg2: (npf × nfint2 × nd)
// Outputs:
//   ind1(i,j2): node index in xdg1 corresponding to node i of face j2 in xdg2
//   ind2(i,j2): face index j1 in xdg1 corresponding to face j2 in xdg2
void match_xdg(std::vector<int>& ind1,
               std::vector<int>& ind2,
               const std::vector<double>& xdg1,
               const std::vector<double>& xdg2,
               int npf, int nfint1, int nfint2, int nd,
               double tol = 1e-12)
{
    // Resize outputs: npf × nfint2
    ind1.resize(static_cast<size_t>(npf) * nfint2);
    ind2.resize(static_cast<size_t>(npf) * nfint2);

    // Loop over faces in xdg2
    for (int j2 = 0; j2 < nfint2; ++j2) {

        bool face_found = false;

        // Try to find the matching face in xdg1
        for (int j1 = 0; j1 < nfint1 && !face_found; ++j1) {

            // Compute geometric centers of faces
            std::vector<double> c1(nd, 0.0), c2(nd, 0.0);
            for (int i = 0; i < npf; ++i) {
                for (int d = 0; d < nd; ++d) {
                    c1[d] += xdg1[i + npf * (j1 + nfint1 * d)];
                    c2[d] += xdg2[i + npf * (j2 + nfint2 * d)];
                }
            }
            for (int d = 0; d < nd; ++d) {
                c1[d] /= npf;
                c2[d] /= npf;
            }

            // Compare face centers
            double dist2 = 0.0;
            for (int d = 0; d < nd; ++d)
                dist2 += (c1[d] - c2[d]) * (c1[d] - c2[d]);
            if (dist2 > 1e-8) continue; // different faces

            // Try node-by-node matching
            std::vector<int> temp_ind1(npf);
            bool mismatch = false;

            for (int i2 = 0; i2 < npf; ++i2) {
                double mindist = std::numeric_limits<double>::max();
                int best_i1 = -1;

                for (int i1 = 0; i1 < npf; ++i1) {
                    double dist = 0.0;
                    for (int d = 0; d < nd; ++d) {
                        double diff = xdg1[i1 + npf * (j1 + nfint1 * d)] -
                                      xdg2[i2 + npf * (j2 + nfint2 * d)];
                        dist += diff * diff;
                    }
                    if (dist < mindist) {
                        mindist = dist;
                        best_i1 = i1;
                    }
                }

                if (mindist > tol) {
                    mismatch = true;
                    break;
                }

                temp_ind1[i2] = best_i1;
            }

            if (!mismatch) {
                // Found matching face
                for (int i = 0; i < npf; ++i) {
                    ind1[i + npf * j2] = temp_ind1[i]; // node mapping
                    ind2[i + npf * j2] = j1;           // face index in xdg1
                }
                face_found = true;
            }
        }

        if (!face_found) {
            // No matching face found
            for (int i = 0; i < npf; ++i) {
                ind1[i + npf * j2] = -1;
                ind2[i + npf * j2] = -1;
            }
        }
    }
}

// General nearest-neighbor face matching between xdg1 and xdg2
// xdg1: (npf1 × nfint1 × nd)
// xdg2: (npf2 × nfint2 × nd)
// Output arrays (npf2 × nfint2):
//   ind1(i,j2): index (0..npf1-1) of nearest node in xdg1
//   ind2(i,j2): index (0..nfint1-1) of the face in xdg1 whose node was nearest
void match_xdg_nearest(std::vector<int>& ind1,
                       std::vector<int>& ind2,
                       const std::vector<double>& xdg1,
                       const std::vector<double>& xdg2,
                       int npf1, int nfint1,
                       int npf2, int nfint2,
                       int nd)
{
    ind1.resize(static_cast<size_t>(npf2) * nfint2);
    ind2.resize(static_cast<size_t>(npf2) * nfint2);

    // For each face j2 in xdg2
    for (int j2 = 0; j2 < nfint2; ++j2) {
        // For each point i2 on that face
        for (int i2 = 0; i2 < npf2; ++i2) {

            double best_dist2 = std::numeric_limits<double>::max();
            int best_face = -1;
            int best_node = -1;

            // Coordinates of the query point from xdg2
            double x2[3] = {0.0, 0.0, 0.0};
            for (int d = 0; d < nd; ++d)
                x2[d] = xdg2[i2 + npf2 * (j2 + nfint2 * d)];

            // Search all faces and nodes in xdg1
            for (int j1 = 0; j1 < nfint1; ++j1) {
                for (int i1 = 0; i1 < npf1; ++i1) {
                    double dist2 = 0.0;
                    for (int d = 0; d < nd; ++d) {
                        double diff = xdg1[i1 + npf1 * (j1 + nfint1 * d)] - x2[d];
                        dist2 += diff * diff;
                    }
                    if (dist2 < best_dist2) {
                        best_dist2 = dist2;
                        best_face = j1;
                        best_node = i1;
                    }
                }
            }

            ind1[i2 + npf2 * j2] = best_node;
            ind2[i2 + npf2 * j2] = best_face;
        }
    }
}

#endif        
