/*
 * Domain Decomposition Utilities for Parallel PDE Solvers
 *
 * This header provides structures and functions to support domain decomposition
 * for parallel finite element methods, including LDG and HDG schemes.
 *
 * Structures:
 * -----------
 * struct DMD
 *   - nbsd: List of neighboring processor IDs.
 *   - elemrecv: For each element received, stores [sender, recv_local_idx, sender_global_idx].
 *   - elemsend: For each element sent, stores [receiver, send_local_idx, recv_global_idx].
 *   - elempart: Local element IDs in the partition.
 *   - elem2cpu: Processor ID for each element in the partition.
 *   - elemsendpts: Number of elements sent to each neighbor.
 *   - elemrecvpts: Number of elements received from each neighbor.
 *   - elempartpts: Partition sizes: [interior, interface, exterior].
 *   - intepartpts: Optional, for HDG: [interior, interface1, interface2, exterior].
 *
 * Functions:
 * ----------
 * DMD initializeDMD(const PDE& pde, const Mesh& mesh)
 *   - Initializes a DMD structure for a single-process or coupled interface scenario.
 *
 * void xiny(std::vector<int>& indices, const std::vector<int>& x, const std::vector<int>& y)
 *   - For each element in x, finds its index in y, or -1 if not found.
 *
 * std::vector<int> neighboringelements(const int* e2e, const int* elem, int nfe, int ne, int nelem)
 *   - Returns a sorted, unique list of neighboring elements for a given set of elements.
 *
 * void create_elemsend(std::vector<DMD>& dmd)
 *   - Populates the elemsend field for each DMD partition, based on elemrecv.
 *   - Also computes elemsendpts and elemrecvpts for communication sizing.
 *
 * void build_dmdldg(std::vector<DMD>& dmd, const int* e2e, const int* elem2cpu, int nfe, int ne, const PDE& pde)
 *   - Builds domain decomposition for LDG schemes.
 *   - Partitions elements into interior, interface, exterior, and outer categories.
 *   - Populates communication structures for parallel exchange.
 *
 * void build_dmdhdg(std::vector<DMD>& dmd, const int* e2e, const int* elem2cpu, const int* inte, int nfe, int ne, const PDE& pde)
 *   - Builds domain decomposition for HDG schemes, supporting coupled interfaces.
 *   - Partitions elements into interior, interface1, interface2, and exterior categories.
 *   - Populates communication structures for parallel exchange.
 *
 * Notes:
 * ------
 * - Many functions use STL algorithms for set operations (difference, intersection, union).
 * - Debug output and file writing are available but commented out.
 * - The code assumes existence of PDE and Mesh types, and that mesh connectivity is provided via e2e arrays.
 */

#ifndef __DOMAINDECOMPOSITION
#define __DOMAINDECOMPOSITION

// struct DMD {
//     std::vector<int> nbsd;                  // neighbors    
//     std::vector<std::vector<int>> elemrecv; // each row: [sender, recv_local_idx, sender_global_idx]
//     std::vector<std::vector<int>> elemsend; // each row: [receiver, send_local_idx, recv_global_idx]
//     std::vector<int> elempart;              // local element IDs in the partition
//     std::vector<int> elem2cpu;              // processor ID for each element in the partition
//     std::vector<int> elemsendpts;           // number of elements sent to each neighbor
//     std::vector<int> elemrecvpts;           // number of elements received from each neighbor
//     std::vector<int> elempartpts;           // partition sizes: [interior, interface, exterior]
//     std::vector<int> intepartpts;           // optional: [interior, interface1, interface2, exterior]
//     std::vector<int> nbinfo;                // neighboring information 
//     int numneigh;                           // number of neighbors 
// //     std::vector<int> inte;                  // processor ID for each element in the partition
// //     std::vector<int> intl;                  // processor ID for each element in the partition
// };

DMD initializeDMD(const PDE& pde, const Mesh& mesh)
{              
    DMD dmd;
    
    if (pde.mpiprocs==1) {
      int ne = mesh.ne;
      if (pde.coupledinterface>0) {                
        int n = mesh.inte.size();

        vector<int> tm(ne);
        vector<int> tn(ne-n);
        for (int i=0; i<ne; i++) tm[i] = i;

        std::set_difference(tm.begin(), tm.end(), mesh.inte.begin(), mesh.inte.end(), std::back_inserter(tn));
        for (int i=0; i<ne-n; i++) dmd.elempart[i] = tn[i];
        for (int i=0; i<n; i++) dmd.elempart[ne-n+i] = mesh.inte[i];

        dmd.elempartpts.resize(3);        
        dmd.elempartpts[0] = ne-n;  dmd.elempartpts[1] = n; dmd.elempartpts[2] = 0;
        dmd.intepartpts.resize(4);        
        dmd.intepartpts[0] = ne-n;  dmd.intepartpts[1] = n; dmd.intepartpts[2] = 0; dmd.intepartpts[3] = 0;            
      }
      else {
        dmd.elempart.resize(ne);
        for (int i=0; i<ne; i++) dmd.elempart[i] = i;
        dmd.elempartpts.resize(1);
        dmd.elempartpts[0] = ne;       
        dmd.intepartpts.resize(4);      
        dmd.intepartpts[0] = ne;  dmd.intepartpts[1] = 0; dmd.intepartpts[2] = 0; dmd.intepartpts[3] = 0;            
      }    
    }
    
    return dmd;
}

void xiny(std::vector<int>& indices, const std::vector<int>& x, const std::vector<int>& y) {
    indices.resize(x.size());
    for (int i = 0; i < x.size(); ++i) {
        auto it = std::find(y.begin(), y.end(), x[i]);
        indices[i] = (it != y.end()) ? std::distance(y.begin(), it) : -1;
    }
}

std::vector<int> neighboringelements(const int* e2e, const int* elem, int nfe, int ne, int nelem) {
    std::vector<int> neighbors;
    neighbors.reserve(nelem * nfe);

    for (int k = 0; k < nelem; ++k) {
        int e = elem[k];
        if (e < 0 || e >= ne) continue;
        for (int i = 0; i < nfe; ++i) {
            int nb = e2e[i + e * nfe];
            if (nb >= 0) {
                neighbors.push_back(nb);
            }
        }
    }

    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());

    return neighbors;
}

void create_elemsend(std::vector<DMD>& dmd) 
{
    int nproc = dmd.size();

    for (int k = 0; k < nproc; ++k) {
        dmd[k].elemsend.clear();
    }

    for (int i = 0; i < nproc; ++i) {
        for (int j = 0; j < dmd[i].nbsd.size(); ++j) {
            int k = dmd[i].nbsd[j];

            //std::vector<std::vector<int>> tm;
            std::vector<std::array<int, 3>> tm;
            tm.reserve(dmd[i].elemrecv.size());
            for (const auto& row : dmd[i].elemrecv) {
                if (row[0] == k) tm.push_back(row);
            }

            for (auto& row : tm) row[0] = i;

            std::vector<int> x(tm.size());
            for (int m = 0; m < tm.size(); ++m)
                x[m] = tm[m][2];

            std::vector<int> local_indices;
            local_indices.reserve(x.size());
            xiny(local_indices, x, dmd[k].elempart);

            for (int m = 0; m < tm.size(); ++m)
                tm[m][1] = local_indices[m];

            dmd[k].elemsend.insert(dmd[k].elemsend.end(), tm.begin(), tm.end());
        }
    }

    for (int i = 0; i < nproc; ++i) {
        dmd[i].elemsendpts.resize(dmd[i].nbsd.size());
        dmd[i].elemrecvpts.resize(dmd[i].nbsd.size());

        for (int j = 0; j < dmd[i].nbsd.size(); ++j) {
            int n = dmd[i].nbsd[j];

            // dmd[i].elemsendpts[j] = std::count_if(
            //     dmd[i].elemsend.begin(), dmd[i].elemsend.end(),
            //     [n](const std::vector<int>& row) { return row[0] == n; });
            // 
            // dmd[i].elemrecvpts[j] = std::count_if(
            //     dmd[i].elemrecv.begin(), dmd[i].elemrecv.end(),
            //     [n](const std::vector<int>& row) { return row[0] == n; });
            dmd[i].elemsendpts[j] = std::count_if(
                dmd[i].elemsend.begin(), dmd[i].elemsend.end(),
                [n](const std::array<int, 3>& row) { return row[0] == n; });
            
            dmd[i].elemrecvpts[j] = std::count_if(
                dmd[i].elemrecv.begin(), dmd[i].elemrecv.end(),
                [n](const std::array<int, 3>& row) { return row[0] == n; });
        }
    }
}

void build_dmdldg(std::vector<DMD>& dmd, const int* e2e, const int* elem2cpu, int nfe, int ne, const PDE& pde) 
{
    int nproc = static_cast<int>(dmd.size());

    for (int i = 0; i < nproc; ++i) {
        std::vector<int> intelem;
        intelem.reserve(ne / nproc);
        for (int e = 0; e < ne; ++e) {
            if (elem2cpu[e] == i) intelem.push_back(e);
        }

        std::vector<int> elem = neighboringelements(e2e, intelem.data(), nfe, ne, intelem.size());
                
        std::vector<int> extelem;
        extelem.reserve(elem.size());
        std::set_difference(elem.begin(), elem.end(), intelem.begin(), intelem.end(), std::back_inserter(extelem));

        std::vector<int> elem2 = neighboringelements(e2e, extelem.data(), nfe, ne, extelem.size());
        
        std::vector<int> bndelem;
        bndelem.reserve(elem2.size());
        std::set_intersection(elem2.begin(), elem2.end(), intelem.begin(), intelem.end(), std::back_inserter(bndelem));
        
        std::vector<int> intex;
        std::set_union(intelem.begin(), intelem.end(), extelem.begin(), extelem.end(), std::back_inserter(intex));
        std::vector<int> outelem;
        outelem.reserve(elem2.size());
        std::set_difference(elem2.begin(), elem2.end(), intex.begin(), intex.end(), std::back_inserter(outelem));
        
        std::vector<int> part1, part2, part3, part4;
        part1.reserve(intelem.size());
        part2 = bndelem;
        part3 = extelem;
        part4 = outelem;

        std::set_difference(intelem.begin(), intelem.end(), bndelem.begin(), bndelem.end(), std::back_inserter(part1));

        dmd[i].elempart.clear();
        dmd[i].elempart.reserve(part1.size() + part2.size() + part3.size() + part4.size());
        dmd[i].elempart.insert(dmd[i].elempart.end(), part1.begin(), part1.end());
        dmd[i].elempart.insert(dmd[i].elempart.end(), part2.begin(), part2.end());
        dmd[i].elempart.insert(dmd[i].elempart.end(), part3.begin(), part3.end());
        dmd[i].elempart.insert(dmd[i].elempart.end(), part4.begin(), part4.end());
        
        dmd[i].elempartpts = {static_cast<int>(part1.size()), static_cast<int>(part2.size()), static_cast<int>(part3.size()), static_cast<int>(part4.size())};

        dmd[i].elem2cpu.resize(dmd[i].elempart.size());
        for (int j = 0; j < dmd[i].elempart.size(); ++j)
            dmd[i].elem2cpu[j] = elem2cpu[dmd[i].elempart[j]];

        std::vector<int> recvelem = extelem;
        recvelem.insert(recvelem.end(), outelem.begin(), outelem.end());
        
        std::vector<int> ind;
        xiny(ind, recvelem, dmd[i].elempart);

        dmd[i].elemrecv.clear();
        dmd[i].elemrecv.reserve(recvelem.size());
        int offset = part1.size() + part2.size();
        for (int j = 0; j < recvelem.size(); ++j) {
            if (ind[j] >= 0) {
                dmd[i].elemrecv.push_back({dmd[i].elem2cpu[ind[j]], static_cast<int>(offset + j), recvelem[j]});
            }
        }

        std::sort(dmd[i].elemrecv.begin(), dmd[i].elemrecv.end());

        dmd[i].nbsd.clear();
        dmd[i].nbsd.reserve(dmd[i].elemrecv.size());
        for (const auto& row : dmd[i].elemrecv)
            dmd[i].nbsd.push_back(row[0]);

        std::sort(dmd[i].nbsd.begin(), dmd[i].nbsd.end());
        dmd[i].nbsd.erase(std::unique(dmd[i].nbsd.begin(), dmd[i].nbsd.end()), dmd[i].nbsd.end());
                
//         if (pde.debugmode==1) {
//           int n = dmd[i].elemrecv.size();
//           vector<int> elemrecv(n);
//           for (int j=0; j<n; j++) {
//             elemrecv[j] = dmd[i].elemrecv[j][1];
//             //elemrecv[j + n*0] = dmd[i].elemrecv[j][0];
//             //elemrecv[j + n*1] = dmd[i].elemrecv[j][1];
//             //elemrecv[j + n*2] = dmd[i].elemrecv[j][2];
//           }
//             
//           writearray2file(pde.datapath + "/intelem" + std::to_string(i) + ".bin", intelem.data(), intelem.size());
//           writearray2file(pde.datapath + "/bndelem" + std::to_string(i) + ".bin", bndelem.data(), bndelem.size());
//           writearray2file(pde.datapath + "/extelem" + std::to_string(i) + ".bin", extelem.data(), extelem.size());
//           writearray2file(pde.datapath + "/outelem" + std::to_string(i) + ".bin", outelem.data(), outelem.size());
//           writearray2file(pde.datapath + "/elempart" + std::to_string(i) + ".bin", dmd[i].elempart.data(), dmd[i].elempart.size());
//           writearray2file(pde.datapath + "/elempartpts" + std::to_string(i) + ".bin", dmd[i].elempartpts.data(), dmd[i].elempartpts.size());
//           writearray2file(pde.datapath + "/elem2cpu" + std::to_string(i) + ".bin", dmd[i].elem2cpu.data(), dmd[i].elem2cpu.size());
//           writearray2file(pde.datapath + "/nbsd" + std::to_string(i) + ".bin", dmd[i].nbsd.data(), dmd[i].nbsd.size());
//           writearray2file(pde.datapath + "/elemrecv" + std::to_string(i) + ".bin", elemrecv.data(), elemrecv.size());
//         }
    }    
    
    create_elemsend(dmd);
    
//     for (int i = 0; i < nproc; ++i) {
//         if (pde.debugmode==1) {
//           int n = dmd[i].elemsend.size();
//           vector<int> elemsend(n);
//           for (int j=0; j<n; j++) elemsend[j] = dmd[i].elemsend[j][1];
//                       
//           writearray2file(pde.datapath + "/elemsend" + std::to_string(i) + ".bin", elemsend.data(), elemsend.size());
//           writearray2file(pde.datapath + "/elemsendpts" + std::to_string(i) + ".bin", dmd[i].elemsendpts.data(), dmd[i].elemsendpts.size());
//           writearray2file(pde.datapath + "/elemrecvpts" + std::to_string(i) + ".bin", dmd[i].elemrecvpts.data(), dmd[i].elemrecvpts.size());         
//         }
//     }
    
    std::cout << "Finished build_dmdldg.\n";
}

void build_dmdhdg(std::vector<DMD>& dmd, const int* e2e, const int* elem2cpu, const int* inte, int nfe, int ne, const PDE& pde) 
{
    int coupledinterface = pde.coupledinterface;
    int nproc = static_cast<int>(dmd.size());

    for (int i = 0; i < nproc; ++i) {
        std::vector<int> intelem;
        intelem.reserve(ne / nproc);
        for (int e = 0; e < ne; ++e) {
            if (elem2cpu[e] == i) intelem.push_back(e);
        }

        std::vector<int> elem = neighboringelements(e2e, intelem.data(), nfe, ne, intelem.size());
        std::vector<int> extelem;
        extelem.reserve(elem.size());
        std::set_difference(elem.begin(), elem.end(), intelem.begin(), intelem.end(), std::back_inserter(extelem));

        elem = neighboringelements(e2e, extelem.data(), nfe, ne, extelem.size());
        std::vector<int> bndelem;
        bndelem.reserve(elem.size());
        std::set_intersection(elem.begin(), elem.end(), intelem.begin(), intelem.end(), std::back_inserter(bndelem));

        std::vector<int> bndelem1, bndelem2;
        if (coupledinterface > 0) {
            bndelem1.reserve(bndelem.size());
            bndelem2.reserve(bndelem.size());
            std::set_intersection(bndelem.begin(), bndelem.end(), inte, inte + ne, std::back_inserter(bndelem1));
            std::set_difference(bndelem.begin(), bndelem.end(), bndelem1.begin(), bndelem1.end(), std::back_inserter(bndelem2));
            bndelem = bndelem1;
            bndelem.insert(bndelem.end(), bndelem2.begin(), bndelem2.end());
        }

        std::vector<int> part1, part2;
        part1.reserve(intelem.size());
        part2 = bndelem;
        std::set_difference(intelem.begin(), intelem.end(), bndelem.begin(), bndelem.end(), std::back_inserter(part1));
        std::vector<int>& part3 = extelem;

        dmd[i].elempart.clear();
        dmd[i].elempart.reserve(part1.size() + part2.size() + part3.size());
        dmd[i].elempart.insert(dmd[i].elempart.end(), part1.begin(), part1.end());
        dmd[i].elempart.insert(dmd[i].elempart.end(), part2.begin(), part2.end());
        dmd[i].elempart.insert(dmd[i].elempart.end(), part3.begin(), part3.end());

        dmd[i].elempartpts = {static_cast<int>(part1.size()), static_cast<int>(part2.size()), static_cast<int>(part3.size())};

        if (coupledinterface > 0) {
            dmd[i].intepartpts = {static_cast<int>(part1.size()), static_cast<int>(bndelem1.size()), static_cast<int>(bndelem2.size()), static_cast<int>(part3.size())};
        }

        dmd[i].elem2cpu.resize(dmd[i].elempart.size());
        for (int j = 0; j < dmd[i].elempart.size(); ++j)
            dmd[i].elem2cpu[j] = elem2cpu[dmd[i].elempart[j]];

        std::vector<int> recvelem = extelem;
        std::vector<int> ind;
        xiny(ind, recvelem, dmd[i].elempart);

        dmd[i].elemrecv.clear();
        dmd[i].elemrecv.reserve(recvelem.size());
        int offset = part1.size() + part2.size();
        for (int j = 0; j < recvelem.size(); ++j) {
            if (ind[j] >= 0) {
                dmd[i].elemrecv.push_back({dmd[i].elem2cpu[ind[j]], static_cast<int>(offset + j), recvelem[j]});
            }
        }

        std::sort(dmd[i].elemrecv.begin(), dmd[i].elemrecv.end());

        dmd[i].nbsd.clear();
        dmd[i].nbsd.reserve(dmd[i].elemrecv.size());
        for (const auto& row : dmd[i].elemrecv)
            dmd[i].nbsd.push_back(row[0]);

        std::sort(dmd[i].nbsd.begin(), dmd[i].nbsd.end());
        dmd[i].nbsd.erase(std::unique(dmd[i].nbsd.begin(), dmd[i].nbsd.end()), dmd[i].nbsd.end());
        
//         if (pde.debugmode==1) {
//           int n = dmd[i].elemrecv.size();
//           vector<int> elemrecv(n);
//           for (int j=0; j<n; j++) {
//             elemrecv[j] = dmd[i].elemrecv[j][1];
//             //elemrecv[j + n*0] = dmd[i].elemrecv[j][0];
//             //elemrecv[j + n*1] = dmd[i].elemrecv[j][1];
//             //elemrecv[j + n*2] = dmd[i].elemrecv[j][2];
//           }
//             
//           writearray2file(pde.datapath + "/intelem" + std::to_string(i) + ".bin", intelem.data(), intelem.size());
//           writearray2file(pde.datapath + "/bndelem" + std::to_string(i) + ".bin", bndelem.data(), bndelem.size());
//           writearray2file(pde.datapath + "/extelem" + std::to_string(i) + ".bin", extelem.data(), extelem.size());
//           writearray2file(pde.datapath + "/elempart" + std::to_string(i) + ".bin", dmd[i].elempart.data(), dmd[i].elempart.size());
//           writearray2file(pde.datapath + "/elempartpts" + std::to_string(i) + ".bin", dmd[i].elempartpts.data(), dmd[i].elempartpts.size());
//           writearray2file(pde.datapath + "/elem2cpu" + std::to_string(i) + ".bin", dmd[i].elem2cpu.data(), dmd[i].elem2cpu.size());
//           writearray2file(pde.datapath + "/nbsd" + std::to_string(i) + ".bin", dmd[i].nbsd.data(), dmd[i].nbsd.size());
//           writearray2file(pde.datapath + "/elemrecv" + std::to_string(i) + ".bin", elemrecv.data(), elemrecv.size());
//         }
    }
    
    create_elemsend(dmd);
    
//     for (int i = 0; i < nproc; ++i) {
//         if (pde.debugmode==1) {
//           int n = dmd[i].elemsend.size();
//           vector<int> elemsend(n);
//           for (int j=0; j<n; j++) elemsend[j] = dmd[i].elemsend[j][1];
//                       
//           writearray2file(pde.datapath + "/elemsend" + std::to_string(i) + ".bin", elemsend.data(), elemsend.size());
//           writearray2file(pde.datapath + "/elemsendpts" + std::to_string(i) + ".bin", dmd[i].elemsendpts.data(), dmd[i].elemsendpts.size());
//           writearray2file(pde.datapath + "/elemrecvpts" + std::to_string(i) + ".bin", dmd[i].elemrecvpts.data(), dmd[i].elemrecvpts.size());         
//         }
//     }
    
    std::cout << "Finished build_dmdhdg.\n";
}

#endif

