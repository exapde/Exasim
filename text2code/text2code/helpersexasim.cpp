/*
    helpers.cpp

    Utility functions and macros for file I/O, array manipulation, and LAPACK bindings.

    LAPACK Bindings:
    - DGEMM, DGETRF, DGETRI: For matrix operations via Fortran LAPACK routines.

    Macros:
    - CPUFREE(x): Safely frees memory and sets pointer to nullptr.
    - error(errmsg): Prints error message with file and line, then exits.

    File and Path Utilities:
    - trim_dir(s): Returns parent directory of a given path.
    - ensure_dir(dir): Ensures a directory exists, creates if needed.
    - make_path(str1, str2): Joins two paths, normalizes result.
    - trimToSubstringAtFirstOccurence(fullPath, keyword): Trims string at first occurrence of keyword.
    - trimToSubstringAtLastOccurence(fullPath, keyword): Trims string at last occurrence of keyword.

    Array Printing:
    - print2darray(a, m, n): Prints m x n double array in scientific format.
    - print2iarray(a, m, n): Prints m x n int array.

    Array File I/O:
    - writearray2file(filename, a, N): Writes array to binary file.
    - readarray(in, a, N): Reads array from binary stream into vector.
    - readiarrayfromdouble(in, a, N): Reads N doubles from file, rounds and stores as ints.
    - readfile(filename, ndims, data, m): Reads dimensions and data from binary file (int/double overloads).

    Array Manipulation:
    - xiny(out, A, B, m, n, dim, tol): Finds indices of A in B (row-wise, column-major).
    - xiny2(out, A, B, m, n, dim, tol): Finds indices of A in B (row-wise, alternate layout).
    - cumsum_int(in, out, n): Computes cumulative sum of int array.
    - unique_ints(arr, n): Sorts array and removes duplicates, returns new size.
    - extract_subset(b, a, ind, k): Extracts subset from array a using indices ind.
    - find(a, b, m, n, k, opts): Counts occurrences or matches in array with options.
    - find(indices, a, b, m, n, k, opts): Finds indices of matches in array with options.
    - simple_bubble_sort(b, ind, a, n): Sorts array a, returns sorted array and indices.
    - unique_count(b, c, a, n): Counts unique values and their occurrences.
    - select_columns(a_new, a, ind, m, k): Selects columns from array a using indices ind (int/double overloads).
    - permute_columns(a, ind, m, k): Permutes columns of array a according to ind.

    Notes:
    - Uses C++ STL (vector, string, filesystem).
    - Some functions use raw pointers and manual memory management.
    - Designed for scientific computing and data manipulation tasks.
*/
#ifndef __HELPERS
#define __HELPERS


// std::string trimToSubstringAtFirstOccurence(const std::string& fullPath, const std::string& keyword) {
//     std::size_t pos = fullPath.find(keyword);  // Use find to get the first occurrence
//     if (pos != std::string::npos) {
//         return fullPath.substr(0, pos + keyword.length());
//     }
//     else {      
//       return "";
//     }
// }
// 
// std::string trimToSubstringAtFirstOccurence(const std::filesystem::path& fullPath, const std::string& keyword) {
//     const std::string s = fullPath.generic_string();
//     std::size_t pos = s.find(keyword);  // Use find to get the first occurrence
//     if (pos != std::string::npos) {
//         return s.substr(0, pos + keyword.length());
//     }
//     else {      
//       return "";
//     }
// }
// 
// std::string trimToSubstringAtLastOccurence(const std::string& fullPath, const std::string& keyword) {
//     std::size_t pos = fullPath.rfind(keyword);  // Use rfind to get the last occurrence
//     if (pos != std::string::npos) {
//         return fullPath.substr(0, pos + keyword.length());
//     }
//     else {      
//       return "";
//     }
// }
// 
// std::string trimToSubstringAtLastOccurence(const std::filesystem::path& fullPath,
//                                            const std::string& keyword)
// {
//     // generic_string uses forward slashes on all platforms (nice for substring ops)
//     const std::string s = fullPath.generic_string();
//     const auto pos = s.rfind(keyword);
//     if (pos != std::string::npos)
//         return s.substr(0, pos + keyword.size());
//     return {};
// }
// 
// void print2darray(const double* a, int m, int n)
// {
//     //cout.precision(4);
//     for (int i=0; i<m; i++) {
//         for (int j=0; j<n; j++)
//             cout << scientific << a[j*m+i] << "   ";
//         cout << endl;
//     }
//     cout << endl;
// }
// 
// void print2iarray(const int* a, int m, int n)
// {
//     for (int i=0; i<m; i++) {
//         for (int j=0; j<n; j++)
//             cout << a[j*m+i] << "   ";
//         cout << endl;
//     }
//     cout << endl;
// }
// 
// template <typename T> void writearray2file(string filename, T *a, int N)
// {
//     if (N>0) {
//         // Open file to read
//         ofstream out(filename.c_str(), ios::out | ios::binary);
// 
//         if (!out) {
//             error("Unable to open file " + filename);
//         }
// 
//         out.write( reinterpret_cast<char*>( &a[0] ), sizeof(T) * N );
// 
//         out.close();
//     }
// }

template <typename T> 
void readarray(ifstream &in, vector<T> &a, int N)
{    
    if (N>0) {        
        a.resize(N);
        in.read( reinterpret_cast<char*>( a.data() ), sizeof(T)*N );        
    }    
}

void readiarrayfromdouble(ifstream &in, vector<int> &a, int N) {
    if (N > 0) {
        a.resize(N);
        double read;
        for (int i = 0; i < N; i++) {
            in.read(reinterpret_cast<char*>(&read), sizeof(read));
            a[i] = static_cast<int>(round(read));
        }
    }
}

void readfile(string filename, vector<int> &ndims, vector<int> &data, int m) 
{
    ifstream in(filename, ios::in | ios::binary);
    if (!in) error("Unable to open file " + filename);
    
    int n = 1;
    readiarrayfromdouble(in, ndims, m);        
    for (int i = 0; i<m; i++) n = n * ndims[i];
    readiarrayfromdouble(in, data, n);
    
    in.close();
}    

void readfile(string filename, vector<int> &ndims, vector<double> &data, int m) 
{
    ifstream in(filename, ios::in | ios::binary);
    if (!in) error("Unable to open file " + filename);
    
    int n = 1;
    readiarrayfromdouble(in, ndims, m);       
    for (int i = 0; i<m; i++) n = n * ndims[i];
    readarray(in, data, n);
    
    in.close();
}    

// template <typename T>
// void readarray(ifstream &in, T **a, int N) {
//     if (N > 0) {
//         *a = (T*) malloc(sizeof(T) * N);
//         in.read(reinterpret_cast<char*>(*a), sizeof(T) * N);
//     }
// }
// 
// int* readiarrayfromdouble(ifstream &in, int N) {
//     int *a = nullptr;
//     if (N > 0) {
//         a = (int*) malloc(sizeof(int) * N);
//         double read;
//         for (int i = 0; i < N; i++) {
//             in.read(reinterpret_cast<char*>(&read), sizeof(read));
//             a[i] = static_cast<int>(round(read));
//         }
//     }
//     return a;
// }

template<typename T>
void xiny(int* out, const T* A, const T* B, int m, int n, int dim, double tol = 1e-12) {
    for (int i = 0; i < m; ++i) {
        out[i] = -1;
        for (int j = 0; j < n; ++j) {
            bool match = true;
            for (int d = 0; d < dim; ++d) {
                if (std::abs(A[i * dim + d] - B[j * dim + d]) > tol) {
                    match = false;
                    break;
                }
            }
            if (match) {
                out[i] = j;
                break;
            }
        }
    }
}

// template<typename T>
// void xiny2(int* out, const T* A, const T* B, int m, int n, int dim, double tol = 1e-12) {
//     for (int i = 0; i < m; ++i) {
//         out[i] = -1;
//         for (int j = 0; j < n; ++j) {
//             bool match = true;
//             for (int d = 0; d < dim; ++d) {
//                 if (std::abs(A[i + n*d] - B[j + m*d]) > tol) {
//                     match = false;
//                     break;
//                 }
//             }
//             if (match) {
//                 out[i] = j;
//                 break;
//             }
//         }
//     }
// }

void cumsum_int(const int* in, int* out, int n) {
    if (n <= 0) return;
    out[0] = in[0];
    for (int i = 1; i < n; ++i) {
        out[i] = out[i - 1] + in[i];
    }
}

// Sort and remove duplicates
// int unique_ints(int* arr, int n) {
//     // Simple insertion sort
//     for (int i = 1; i < n; ++i) {
//         int key = arr[i];
//         int j = i - 1;
//         while (j >= 0 && arr[j] > key) {
//             arr[j + 1] = arr[j];
//             j--;
//         }
//         arr[j + 1] = key;
//     }
//     // Remove duplicates
//     int m = 0;
//     for (int i = 0; i < n; ++i) {
//         if (m == 0 || arr[i] != arr[m - 1]) {
//             arr[m++] = arr[i];
//         }
//     }
//     return m;
// }

// int unique_ints(int* arr, int n) {
//     if (n <= 1) return n;
//     std::sort(arr, arr + n);                       // introsort (quick+heap+insertion)
//     int* last = std::unique(arr, arr + n);         // compacts uniques to front
//     return int(last - arr);                        // new length m
// }

void extract_subset(int* b, const int* a, const int* ind, size_t k) 
{
    for (size_t i = 0; i < k; ++i) {
        b[i] = a[ind[i]];
    }
}

// int find(const int* a, int b, int m, int n, int k, int opts) 
// {
//     int count = 0;
//     if (opts==0) {
//       for (int i = 0; i < n; ++i) {
//           if (a[k + i*m] == b) count++;                        
//       }
//     }
//     else if (opts==1) {
//       for (int i = 0; i < n; ++i) {
//           if (a[k + i*m] <= b) count++;
//       }
//     }
//     else if (opts==2) {
//       for (int i = 0; i < n; ++i) {
//           if (a[k + i*m] >= b) count++;
//       }
//     }    
//     return count;
// }

// int find(int* indices, const int* a, int b, int m, int n, int k, int opts) 
// {
//     int count = 0;
//     if (opts==0) {
//       for (int i = 0; i < n; ++i) {
//           if (a[k + i*m] == b) {
//               indices[count++] = i;
//           }
//       }
//     }
//     else if (opts==1) {
//       for (int i = 0; i < n; ++i) {
//           if (a[k + i*m] <= b) {
//               indices[count++] = i;
//           }
//       }
//     }
//     else if (opts==2) {
//       for (int i = 0; i < n; ++i) {
//           if (a[k + i*m] >= b) {
//               indices[count++] = i;
//           }
//       }
//     }
// 
//     return count;
// }

// void simple_bubble_sort(int* b, int* ind, const int* a, int n) 
// {
//     // Initialize b and ind
//     for (int i = 0; i < n; ++i) {
//         b[i] = a[i];
//         ind[i] = i;
//     }
//     // Bubble sort
//     for (int i = 0; i < n-1; ++i) {
//         for (int j = 0; j < n-i-1; ++j) {
//             if (b[j] > b[j+1]) {
//                 // Swap values
//                 int tmp = b[j];
//                 b[j] = b[j+1];
//                 b[j+1] = tmp;
//                 // Swap indices
//                 int tmp_idx = ind[j];
//                 ind[j] = ind[j+1];
//                 ind[j+1] = tmp_idx;
//             }
//         }
//     }
// }
// 
// int unique_count(int *b, int *c, const int *a, int n)
// {
//     if (n == 0) return 0;
// 
//     int uniq = 0;          /* index in b/c */
//     double current = a[0];
//     int     count   = 1;
// 
//     for (int i = 1; i < n; ++i) {
//         if (a[i] == current) {
//             ++count;          /* same value – just bump counter */
//         } else {
//             /* flush previous run */
//             b[uniq] = current;
//             c[uniq] = count;
//             ++uniq;
// 
//             /* start new run */
//             current = a[i];
//             count   = 1;
//         }
//     }
//     /* flush final run */
//     b[uniq] = current;
//     c[uniq] = count;
//     ++uniq;
// 
//     return uniq;              /* number of unique values */
// }
// 
// void select_columns(int* a_new, const int* a, const int* ind, int m, int k) 
// {
//     for (int j = 0; j < k; ++j) {
//         int col = ind[j];
//         for (int i = 0; i < m; ++i) {
//             a_new[i + j * m] = a[i + col * m];
//         }
//     }
// }

void select_columns(double* a_new, const double* a, const int* ind, int m, int k) 
{
    for (int j = 0; j < k; ++j) {
        int col = ind[j];
        for (int i = 0; i < m; ++i) {
            a_new[i + j * m] = a[i + col * m];
        }
    }
}

// void permute_columns(int* a, const int* ind, int m, int k) 
// {    
//   if ((k > 0) && (m > 0)) {
//     int* a_new = (int*)malloc(m * k * sizeof(int));
// 
//     for (int j = 0; j < k; ++j) {
//         int col = ind[j];
//         for (int i = 0; i < m; ++i) {
//             a_new[i + j * m] = a[i + col * m];
//         }
//     }
// 
//     for (int j = 0; j < k; ++j)
//       for (int i = 0; i < m; ++i) 
//         a[i + j * m] = a_new[i + j * m];
// 
//     free(a_new);
//   }
// }

#endif

