#include <iostream>

#include "ExasimSolver.hpp"
#include "../Model/BuiltIn/builtinlibprovider.cpp"

namespace {

const char* SelectExasimDriverProviderName()
{
    return "BuiltInLibrary";
}

void PrintModelProvider(const int modelnumber, const int builtinmodelID)
{
    int rank = 0;
#ifdef HAVE_MPI
    MPI_Comm_rank(EXASIM_COMM_WORLD, &rank);
#endif

    if (rank == 0) {
        std::cout << "Model " << modelnumber
                  << ": provider = " << SelectExasimDriverProviderName()
                  << ", builtinmodelID = " << builtinmodelID << std::endl;
    }
}

int ConfigureModelDefinitions(ExasimSolver& solver)
{
    const ExasimDriverABI& abi = getBuiltInLibraryExasimDriverABI();

    for (int i = 0; i < solver.NumModelDefinitions(); i++) {
        const int builtinmodelID = solver.BuiltinModelID(i);
        const int err = solver.SetModelDefinition(i, builtinmodelID, abi);
        if (err) return err;
        PrintModelProvider(i, builtinmodelID);
    }

    return 0;
}

} // namespace

int main(int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Comm comm = MPI_COMM_WORLD;
#else
    MPI_Comm comm = MPI_COMM_NULL;
#endif

    ExasimSolver solver;

    int err = solver.InitializeEnvironment(argc, argv, comm);
    if (err) return err;

    err = solver.ParseInputs(argc, argv);
    if (err) {
        solver.Finalize();
        return err;
    }

    err = ConfigureModelDefinitions(solver);
    if (err) {
        solver.Finalize();
        return err;
    }

    err = solver.InitializeModels();
    if (err) {
        solver.Finalize();
        return err;
    }

    err = solver.Solve();
    if (err) {
        solver.Finalize();
        return err;
    }

    return solver.Finalize();
}
