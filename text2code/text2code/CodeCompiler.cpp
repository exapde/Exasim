
/*
    CodeCompiler.cpp

    This file provides functions for generating, compiling, and building dynamic libraries for C++ code
    based on a parsed model specification. It is designed to support multiple platforms (Windows, macOS, Linux)
    and toolchains (GCC, Clang, MSVC), as well as optional CUDA and HIP backends via Kokkos.

    Main Components:

    1. generateCppCode(PDE& pde)
        - Parses a model file and generates corresponding C++ source files using SymEngine and custom code generators.
        - Handles framework-specific code generation for CUDA and HIP.
        - Generates empty source files for optional outputs if not specified.
        - Outputs generated files to the model directory.

    2. Compiler Detection and Toolchain Utilities
        - Functions to detect available C++ compilers and their kind (GCC, Clang, MSVC).
        - Utilities for quoting paths, running commands silently, and checking compiler support for flags.
        - Toolchain struct encapsulates compiler executable, kind, and relevant flags.

    3. executeCppCode(ParsedSpec& spec)
        - Compiles the generated C++ source file into an executable using the detected toolchain.
        - Links against SymEngine library.
        - Runs the resulting executable to perform code generation tasks.

    4. buildDynamicLibraries(ParsedSpec& spec)
        - Builds shared libraries for the generated model code using Kokkos backends (serial, CUDA, HIP).
        - Handles platform-specific shared library extensions and compiler flags.
        - Compiles and links against Kokkos core and container libraries.
        - Supports detection and usage of nvcc_wrapper and hipcc_wrapper for CUDA and HIP compilation.

    5. Utility Functions
        - getSharedLibExtension(): Returns the appropriate shared library extension for the current platform.
        - with_exe_if_windows(): Appends ".exe" to executable names on Windows.
        - join(): Joins vector of strings with a separator.

    Notes:
        - Error handling is performed via error() calls when compilation or execution fails.
        - The code is designed to be robust across platforms and compilers, with careful quoting and flag selection.
        - Some legacy code for building dynamic libraries is commented out for reference.

*/

void generateCppCode(ParsedSpec spec)
{
    std::cout << "\nGenerating the C++ code for this model file ("<< spec.modelfile << ") ... \n\n";

    // ParsedSpec spec = TextParser::parseFile(make_path(pde.datapath, pde.modelfile));        
    // spec.exasimpath = pde.exasimpath;
    // spec.modelpath = make_path(pde.exasimpath, "/backend/Model/");
    // spec.symenginepath = make_path(pde.exasimpath, "/text2code/symengine");
    
    // Generate SymEngine code
    CodeGenerator gen(spec);
    gen.generateCode2Cpp(make_path(spec.modelpath, "Code2Cpp.cpp"));
    gen.generateSymbolicFunctionsHpp(make_path(spec.modelpath, "SymbolicFunctions.hpp"));
    gen.generateSymbolicFunctionsCpp(make_path(spec.modelpath, "SymbolicFunctions.cpp"));        
    gen.generateSymbolicScalarsVectorsHpp(make_path(spec.modelpath, "SymbolicScalarsVectors.hpp"));
    gen.generateSymbolicScalarsVectorsCpp(make_path(spec.modelpath, "SymbolicScalarsVectors.cpp"));   

    std::string framework = spec.framework;
    if ((framework == "Cuda") || (framework == "CUDA") || (framework == "cuda")) {
      gen.generateCudaHipHpp(make_path(spec.modelpath, "CudaHip.hpp"));
    } 
    else if ((framework == "Hip") || (framework == "HIP") || (framework == "hip")) {
      gen.generateCudaHipHpp(make_path(spec.modelpath, "CudaHip.hpp"));
    } 

    if (spec.isoutput[6]==false) gen.generateEmptySourcewCpp(spec.modelpath);
    if (spec.isoutput[7]==false) gen.generateEmptyOutputCpp(spec.modelpath);
    if (spec.isoutput[8]==false) gen.generateEmptyMonitorCpp(spec.modelpath);
    if (spec.isoutput[9]==false) gen.generateEmptyInituCpp(spec.modelpath);
    if (spec.isoutput[10]==false) gen.generateEmptyInitqCpp(spec.modelpath);
    if (spec.isoutput[11]==false) gen.generateEmptyInitudgCpp(spec.modelpath);
    if (spec.isoutput[12]==false) gen.generateEmptyInitwdgCpp(spec.modelpath);
    if (spec.isoutput[13]==false) gen.generateEmptyInitodgCpp(spec.modelpath);
    if (spec.isoutput[14]==false) gen.generateEmptyAvfieldCpp(spec.modelpath);
    if (spec.isoutput[15]==false) gen.generateEmptyFintCpp(spec.modelpath);
    if (spec.isoutput[16]==false) gen.generateEmptyEoSCpp(spec.modelpath);
    if (spec.isoutput[17]==false) gen.generateEmptyVisScalarsCpp(spec.modelpath);
    if (spec.isoutput[18]==false) gen.generateEmptyVisVectorsCpp(spec.modelpath);
    if (spec.isoutput[19]==false) gen.generateEmptyVisTensorsCpp(spec.modelpath);
    if (spec.isoutput[20]==false) gen.generateEmptyQoIvolumeCpp(spec.modelpath);
    if (spec.isoutput[21]==false) gen.generateEmptyQoIboundaryCpp(spec.modelpath);
    gen.generateEmptyFhatCpp(spec.modelpath);
    gen.generateEmptyUhatCpp(spec.modelpath);
    gen.generateEmptyStabCpp(spec.modelpath);
    gen.generateLibPDEModelHpp(spec.modelpath);
    gen.generateLibPDEModelCpp(spec.modelpath);

    std::cout << "C++ source files are generated and located in this directory " << spec.modelpath << "\n";
    std::cout << "Finished generating code.\n";
    
    //return spec;
}

namespace fs = std::filesystem;

static std::string quote(const std::string& s) {
#ifdef _WIN32
    // crude but effective quoting for paths with spaces
    if (s.find_first_of(" \t\"") == std::string::npos) return s;
    return "\"" + s + "\"";
#else
    if (s.find_first_of(" \t\"'\\$()") == std::string::npos) return s;
    return "'" + s + "'";
#endif
}

static int run_silently(const std::string& cmd) {
#ifdef _WIN32
    std::string full = cmd + " >nul 2>&1";
#else
    std::string full = cmd + " >/dev/null 2>&1";
#endif
    return std::system(full.c_str());
}

enum class CompilerKind { GCC, Clang, MSVC, Unknown };

static bool probe_compiler(const std::string& exe) {
#ifdef _WIN32
    // cl/clang-cl don’t accept --version; /? is fine
    if (exe.find("cl") != std::string::npos)
        return run_silently(quote(exe) + " /?") == 0;
#endif
    // Most compilers accept --version or -v
    int ok = run_silently(quote(exe) + " --version");
    if (ok != 0) ok = run_silently(quote(exe) + " -v");
    return ok == 0;
}

static std::string first_nonempty(const char* s) { return s ? std::string(s) : std::string(); }

static std::string detect_cxx_exe() {
    // 1) respect CXX env
    if (auto env = first_nonempty(std::getenv("CXX")); !env.empty()) {
        if (probe_compiler(env)) return env;
    }
#ifdef _WIN32
    // 2) common Windows candidates
    const char* candidates[] = {"g++", "clang++", "c++", "cl", "clang-cl", nullptr};
#else
    // 2) common Unix/mac candidates (prefer c++ wrapper)
    const char* candidates[] = {"g++", "clang++", "c++", nullptr};
#endif
    for (const char** p = candidates; *p; ++p) {
        if (probe_compiler(*p)) return *p;
    }
    throw std::runtime_error("No working C++ compiler found in PATH. Set CXX or install clang++/g++.");
}

static CompilerKind detect_kind(const std::string& exe) {
#ifdef _WIN32
    if (exe.find("cl") != std::string::npos) return CompilerKind::MSVC;
    if (exe.find("clang-cl") != std::string::npos) return CompilerKind::MSVC; // via clang driver
#endif
    // Cheap heuristic via -v output
    // (We don’t capture here; kind is mostly for MSVC vs not)
    if (exe.find("clang") != std::string::npos) return CompilerKind::Clang;
    if (exe.find("g++")   != std::string::npos || exe == "c++") return CompilerKind::GCC;
    return CompilerKind::Unknown;
}

static bool compiler_supports_flag(const std::string& exe, const std::string& flag) {
    fs::path tmpdir = fs::temp_directory_path();
    fs::path src = tmpdir / "cmk_check.cpp";
    fs::path obj = tmpdir / "cmk_check.o";
    std::ofstream(src.string()) << "int main(){return 0;}\n";
#ifdef _WIN32
    std::string out = quote((tmpdir / "cmk_check.obj").string());
    std::string cmd = quote(exe) + " " + flag + " /nologo /c " + quote(src.string()) + " /Fo" + out;
#else
    std::string cmd = quote(exe) + " -x c++ " + flag + " -c " + quote(src.string()) +
                      " -o " + quote(obj.string());
#endif
    int rc = run_silently(cmd);
    std::error_code ec; fs::remove(src, ec); fs::remove(obj, ec);
    return rc == 0;
}

struct Toolchain {
    std::string cxx;
    CompilerKind kind{};
    std::vector<std::string> cxx17_flag;
    std::vector<std::string> warn_squelch;
};

static Toolchain detect_toolchain() {
    Toolchain tc;
    tc.cxx  = detect_cxx_exe();
    tc.kind = detect_kind(tc.cxx);

    // C++17 flag
    if (tc.kind == CompilerKind::MSVC) {
        tc.cxx17_flag = {"/std:c++17", "/EHsc"};
    } else {
        tc.cxx17_flag = {"-std=c++17"};
    }

    // Only add warning-squelch if supported (Clang has this, GCC doesn’t)
    if (compiler_supports_flag(tc.cxx, "-Wno-inconsistent-missing-override"))
        tc.warn_squelch.push_back("-Wno-inconsistent-missing-override");

    return tc;
}

static std::string join(const std::vector<std::string>& v, const std::string& sep = " ") {
    std::string out;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) out += sep;
        out += v[i];
    }
    return out;
}

constexpr std::string_view to_string(CompilerKind k) {
    switch (k) {
        case CompilerKind::GCC:    return "GCC";
        case CompilerKind::Clang:  return "Clang";
        case CompilerKind::MSVC:   return "MSVC";
        default:                   return "Unknown";
    }
}

static std::string with_exe_if_windows(std::string s) {
#ifdef _WIN32
    fs::path p = s;
    if (p.extension().empty()) p.replace_extension(".exe");
    return p.string();
#else
    return s;
#endif
}

int executeCppCode(ParsedSpec& spec) 
{
    std::string sourcefile = make_path(spec.modelpath , "Code2Cpp.cpp");                    
    std::string symengine_include = make_path(spec.symenginepath, "include");    
    std::string exefile = make_path(spec.modelpath, "code2cpp");

    Toolchain tc = detect_toolchain();

    std::cout<<"Compiling " + sourcefile + "\n";
    
    // Construct compile command using text2code_path
    std::stringstream cmd;
    
    if (tc.kind == CompilerKind::MSVC) {    
      std::string symengine_lib = make_path(spec.symenginepath, "lib/symengine.lib");
      cmd << tc.cxx << " /std:c++17 /EHsc /W0 "
          << "/I" << quote(spec.symenginepath) << " "
          << "/I" << quote(symengine_include) << " "
          << quote(sourcefile) << " "
          << quote(symengine_lib) << " "
          << "/Fe:" << quote(exefile);
    } else {
      std::string symengine_lib = make_path(spec.symenginepath, "lib/libsymengine.a");
      cmd << tc.cxx << " -std=c++17 -w "          
          << "-I" << quote(spec.symenginepath) << " "
          << "-I" << quote(symengine_include) << " "       
          << quote(sourcefile) << " "        
          << quote(symengine_lib) << " -o "    
          << quote(exefile);
    }
                           
    std::cout<<cmd.str().c_str()<<std::endl;
    
    int status = std::system(cmd.str().c_str());
    if (status != 0)  
        error("Compiling Code2Cpp.cpp failed! Please compile text2code with -DUSE_CMAKE=ON and use cmake to run text2code.");        
    else
        std::cout<<"Compiling " + sourcefile + " successfully.\n";
                                  
    std::cout << "Running " << exefile << " ...\n";
    std::string run = with_exe_if_windows(exefile);
    status = std::system(quote(run).c_str());
    if (status != 0) 
        error("code2cpp execution failed! Please compile text2code with -DUSE_CMAKE=ON and use cmake to run text2code.");            
    else
        std::cout<<"Running " + exefile + " successfully.\n";
    
    return 0;
}

std::string getSharedLibExtension() {
#if defined(_WIN32)
    return ".dll";
#elif defined(__APPLE__)
    return ".dylib";
#else // Linux and other POSIX systems
    return ".so";
#endif
}

int buildDynamicLibraries(ParsedSpec& spec) 
{
    const std::string kokkos_serial_path = make_path(spec.exasimpath, "kokkos/buildserial");
    const std::string kokkos_cuda_path   = make_path(spec.exasimpath, "kokkos/buildcuda");
    const std::string kokkos_hip_path    = make_path(spec.exasimpath, "kokkos/buildhip");

    Toolchain tc = detect_toolchain();
    auto ext = getSharedLibExtension();

    auto out_name = [&](const std::string& stem) {
        return make_path(spec.modelpath, stem) + ext;   // build full name, THEN quote
    };

    auto inc = [&](const std::string& dir) {
        return (tc.kind == CompilerKind::MSVC)
               ? ("/I" + quote(dir))
               : ("-I" + quote(dir));
    };

    auto shlib_flags = [&] {
        if (tc.kind == CompilerKind::MSVC) return std::string{"/LD /nologo /std:c++17 /EHsc /W0 "};
#ifdef __APPLE__
        return std::string{"-dynamiclib -std=c++17 -fPIC "};  // macOS
#else
        return std::string{"-shared -std=c++17 -fPIC "};      // Linux/MinGW
#endif
    }();

    int status = 0;

    // -------------------- SERIAL --------------------
    if (fs::exists(kokkos_serial_path) && fs::is_directory(kokkos_serial_path)) {
        const std::string inc_dir = make_path(kokkos_serial_path, "include");
        const std::string src      = make_path(spec.modelpath, "libpdemodel.cpp");
        const std::string out      = out_name("libpdemodelserial");

        std::cout << "Compiling libpdemodelserial\n";
        std::stringstream cmd;

        if (tc.kind == CompilerKind::MSVC) {
            const std::string lib_core = make_path(kokkos_serial_path, "lib/kokkoscore.lib");
            const std::string lib_cont = make_path(kokkos_serial_path, "lib/kokkoscontainers.lib");
            const std::string lib_simd = make_path(kokkos_serial_path, "lib/kokkossimd.lib");
          
            cmd << tc.cxx << " " << shlib_flags
                << inc(spec.modelpath) << " " << inc(inc_dir) << " "
                << quote(src) << " "
                << "/Fe:" << quote(out) << " "
                << "/link " << quote(lib_core) << " " << quote(lib_cont) << " " << quote(lib_simd);
        } else {
            const std::string lib_core = make_path(kokkos_serial_path, "lib/libkokkoscore.a");
            const std::string lib_cont = make_path(kokkos_serial_path, "lib/libkokkoscontainers.a");
            const std::string lib_simd = make_path(kokkos_serial_path, "lib/libkokkossimd.a");
          
            cmd << tc.cxx << " " << shlib_flags
                << inc(spec.modelpath) << " " << inc(inc_dir) << " "  // includes BEFORE source
                << quote(src) << " "
                << quote(lib_core) << " " << quote(lib_cont) << " " << quote(lib_simd) << " "
                << "-o " << quote(out);
        }

        std::cout << cmd.str() << "\n";
        status = std::system(cmd.str().c_str());
        if (status != 0) { error("Compiling libpdemodelserial failed!"); return 1; }
        std::cout << "libpdemodelserial built: " << out << "\n";
    }

    // -------------------- CUDA (sketch) --------------------
    if (fs::exists(kokkos_cuda_path) && fs::is_directory(kokkos_cuda_path)) {
        // Prefer Kokkos nvcc_wrapper if available:
        std::string cuda_cxx = tc.cxx;
        std::string nvccw = make_path(kokkos_cuda_path, "bin/nvcc_wrapper");
        if (fs::exists(nvccw)) cuda_cxx = nvccw;

        const std::string inc_dir = make_path(kokkos_cuda_path, "include");
        const std::string lib_dir = make_path(kokkos_cuda_path, "lib");
        const std::string lib_core = make_path(kokkos_cuda_path, "lib/libkokkoscore.a");
        const std::string lib_cont = make_path(kokkos_cuda_path, "lib/libkokkoscontainers.a");
        const std::string src      = make_path(spec.modelpath, "libpdemodel.cpp");
        const std::string out      = out_name("libpdemodelcuda");

        std::cout << "Compiling libpdemodelcuda\n";
        std::stringstream cmd;
        cmd << quote(cuda_cxx) << " -shared -std=c++17 -Xcompiler -fPIC --expt-extended-lambda --expt-relaxed-constexpr "
            << inc(spec.modelpath) << " " << inc(inc_dir) << " "
            << quote(src) << " -L" << quote(lib_dir) << " -lkokkoscore -lkokkoscontainers -Wl,-rpath,"
            << quote(lib_dir) << " " << "-o " << quote(out);

        std::cout << cmd.str() << "\n";
        status = std::system(cmd.str().c_str());
        if (status != 0) { error("Compiling libpdemodelcuda failed!"); return 1; }
        std::cout << "libpdemodelcuda built: " << out << "\n";
    }

    // -------------------- HIP (sketch) --------------------
    if (fs::exists(kokkos_hip_path) && fs::is_directory(kokkos_hip_path)) {
        std::string hip_cxx = tc.cxx;
        std::string amdccw = make_path(kokkos_cuda_path, "bin/hipcc_wrapper");
        if (fs::exists(amdccw)) hip_cxx = amdccw;
        else {        
          std::string hipcc = "hipcc";
          if (run_silently(quote(hipcc) + " --version") == 0) hip_cxx = hipcc;
        }
        
        const std::string inc_dir = make_path(kokkos_hip_path, "include");
        const std::string lib_core = make_path(kokkos_hip_path, "lib/libkokkoscore.a");
        const std::string lib_cont = make_path(kokkos_hip_path, "lib/libkokkoscontainers.a");
        const std::string src      = make_path(spec.modelpath, "libpdemodel.cpp");
        const std::string out      = out_name("libpdemodelhip");

        std::cout << "Compiling libpdemodelhip\n";
        std::stringstream cmd;
        cmd << quote(hip_cxx) << " -shared -std=c++17 -fPIC "
            << inc(spec.modelpath) << " " << inc(inc_dir) << " "
            << quote(src) << " "
            << quote(lib_core) << " " << quote(lib_cont) << " "
            << "-o " << quote(out);

        std::cout << cmd.str() << "\n";
        status = std::system(cmd.str().c_str());
        if (status != 0) { error("Compiling libpdemodelhip failed!"); return 1; }
        std::cout << "libpdemodelhip built: " << out << "\n";
    }

    return 0;
}








// int buildDynamicLibraries(ParsedSpec& spec) 
// {
//     std::string kokkos_path = spec.exasimpath + "/kokkos";          
//     std::string kokkos_serial_path = spec.exasimpath + "/kokkos/buildserial";
//     std::string kokkos_cuda_path = spec.exasimpath + "/kokkos/buildcuda";
//     std::string kokkos_hip_path = spec.exasimpath + "/kokkos/buildhip";
// 
//     Toolchain tc = detect_toolchain();
//     int status = 0;
//     
//     if (std::filesystem::exists(kokkos_serial_path) && std::filesystem::is_directory(kokkos_serial_path)) {   
//         std::cout<<"Compiling libpdemodelserial\n";
//         std::stringstream cmp;
//         cmp << tc.cxx << " -fPIC -shared -std=c++17 "        
//             << quote(make_path(spec.modelpath, "libpdemodel.cpp"))<<" "   
//             << "-I" << quote(make_path(kokkos_serial_path, "include"))<<" "       
//             << quote(make_path(kokkos_serial_path, "lib/libkokkoscore.a"))<<" "    
//             << quote(make_path(kokkos_serial_path, "lib/libkokkoscontainers.a"))<<" "        
//             << quote(make_path(kokkos_serial_path, "lib/libkokkossimd.a"))<<" -o "        
//             << quote(make_path(spec.modelpath, "libpdemodelserial"))<< getSharedLibExtension();                                   
// 
//         std::cout<<cmp.str().c_str()<<std::endl;
//         
//         status = std::system(cmp.str().c_str());
//         if (status != 0) 
//           error("Compiling libpdemodelserial failed! Please compile text2code with -DUSE_CMAKE=ON and use cmake to run text2code.");        
//         else
//           std::cout<<"Compiling libpdemodelserial successfully.\n";
//     }                   
// 
//     if (std::filesystem::exists(kokkos_cuda_path) && std::filesystem::is_directory(kokkos_cuda_path)) {
//         std::cout<<"Compiling libpdemodelcuda\n";
//         std::stringstream cmp;
//         cmp << tc.cxx << " -fPIC -shared -std=c++17 "     
//             << quote(make_path(spec.modelpath, "libpdemodel.cpp"))<<" "      
//             << "-I" << quote(make_path(kokkos_cuda_path, "include"))<<" "     
//             << quote(make_path(kokkos_cuda_path, "lib/libkokkoscore.a"))<<" " 
//             << quote(make_path(kokkos_cuda_path, "lib/libkokkoscontainers.a"))<<" -o "    
//             << quote(make_path(spec.modelpath, "libpdemodelcuda"))<< getSharedLibExtension();                                    
// 
//         std::cout<<cmp.str().c_str()<<std::endl;
//         
//         status = std::system(cmp.str().c_str());
//         if (status != 0) 
//           error("Compiling libpdemodelcuda failed! Please compile text2code with -DUSE_CMAKE=ON and use cmake to run text2code.");        
//         else
//           std::cout<<"Compiling libpdemodelcuda successfully.\n";
//     }                   
// 
//     if (std::filesystem::exists(kokkos_hip_path) && std::filesystem::is_directory(kokkos_hip_path)) {
//         std::cout<<"Compiling libpdemodelhip\n";
//         std::stringstream cmp;
//         cmp << tc.cxx << " -fPIC -shared -std=c++17 "        
//             << quote(make_path(spec.modelpath, "libpdemodel.cpp"))<<" "   
//             << "-I" << quote(make_path(kokkos_hip_path, "include"))<<" "     
//             << quote(make_path(kokkos_hip_path, "lib/libkokkoscore.a"))<<" " 
//             << quote(make_path(kokkos_hip_path, "lib/libkokkoscontainers.a"))<<" -o "        
//             << quote(make_path(spec.modelpath, "libpdemodelhip"))<< getSharedLibExtension();                              
// 
//         std::cout<<cmp.str().c_str()<<std::endl;
//         
//         status = std::system(cmp.str().c_str());
//         if (status != 0) 
//           error("Compiling libpdemodelhip failed! Please compile text2code with -DUSE_CMAKE=ON and use cmake to run text2code.");        
//         else
//           std::cout<<"Compiling libpdemodelhip successfully.\n";
//     }                             
//     
//     return 0;
// }
// 
