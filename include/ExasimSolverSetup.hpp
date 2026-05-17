/**
 * @file ExasimSolverSetup.hpp
 * @brief High-level helper API for selecting and attaching Exasim model providers.
 *
 * This header centralizes the provider-selection logic that was previously
 * duplicated in executable entry points. It selects the provider ABI from the
 * active compile-time mode and configures every model definition owned by an
 * ExasimSolver instance.
 *
 * The low-level ABI getter functions declared here are defined by the provider
 * modules for FrontendGenerated, Text2codeGenerated, UserDefined, and BuiltIn
 * models. This header intentionally declares them but does not include the
 * provider implementation translation units, so it remains install-safe.
 */

#pragma once

#include <iostream>
#include <vector>

#include "ExasimSolver.hpp"

// Low-level provider ABI entry points supplied by the provider modules.
#ifdef HAVE_BUILTINMODEL
const ExasimDriverABI& getBuiltInExasimDriverABI();
#endif
const ExasimDriverABI& getKokkosKernelExasimDriverABI();
const ExasimDriverABI& getText2codeGeneratedExasimDriverABI();
const ExasimDriverABI& getUserDefinedExasimDriverABI();
const ExasimDriverABI& getFrontendGeneratedExasimDriverABI();

inline const ExasimDriverABI& SelectExasimDriverABI(const int builtinmodelID)
{
#ifdef HAVE_BUILTINMODEL
    if (builtinmodelID > 0)
        return getBuiltInExasimDriverABI();
#endif

#ifdef HAVE_KOKKOSKERNEL
    return getKokkosKernelExasimDriverABI();
#elif defined(HAVE_TEXT2CODE)
    return getText2codeGeneratedExasimDriverABI();
#elif defined(HAVE_UDERDEFINED) || defined(HAVE_USERDEFINED)
    return getUserDefinedExasimDriverABI();
#else
    return getFrontendGeneratedExasimDriverABI();
#endif
}

inline const char* SelectExasimDriverProviderName(const int builtinmodelID)
{
#ifdef HAVE_BUILTINMODEL
    if (builtinmodelID > 0)
        return "BuiltIn";
#endif

#ifdef HAVE_KOKKOSKERNEL
    return "KokkosKernel";
#elif defined(HAVE_TEXT2CODE)
    return "Text2codeGenerated";
#elif defined(HAVE_UDERDEFINED) || defined(HAVE_USERDEFINED)
    return "UserDefined";
#else
    return "FrontendGenerated";
#endif
}

inline void PrintModelProvider(const int modelnumber, const int builtinmodelID)
{
    int rank = 0;
#ifdef HAVE_MPI
    MPI_Comm_rank(EXASIM_COMM_WORLD, &rank);
#endif

    if (rank == 0) {
        std::cout << "Model " << modelnumber
                  << ": provider = "
                  << SelectExasimDriverProviderName(builtinmodelID)
                  << ", builtinmodelID = " << builtinmodelID << std::endl;
    }
}

inline int ConfigureModelDefinitions(ExasimSolver& solver)
{
    const int numModelDefinitions = solver.NumModelDefinitions();

    for (int i = 0; i < numModelDefinitions; i++) {
        const int builtinmodelID = solver.BuiltinModelID(i);
        const int err = solver.SetModelDefinition(
            i, builtinmodelID, SelectExasimDriverABI(builtinmodelID));
        if (err)
            return err;
        PrintModelProvider(i, builtinmodelID);
    }

    return 0;
}

inline int ConfigureModelDefinitions(ExasimSolver& solver,
                                     const std::vector<int>& vecBuiltinModelID)
{
    const int numModelDefinitions = solver.NumModelDefinitions();

    if (static_cast<int>(vecBuiltinModelID.size()) != numModelDefinitions)
        return 1;

    for (int i = 0; i < numModelDefinitions; i++) {
        int builtinmodelID = vecBuiltinModelID[i];

        if (builtinmodelID <= 0)
            builtinmodelID = solver.BuiltinModelID(i);

        const int err = solver.SetModelDefinition(
            i, builtinmodelID, SelectExasimDriverABI(builtinmodelID));
        if (err)
            return err;
        PrintModelProvider(i, builtinmodelID);
    }

    return 0;
}

inline int RunExasimSolver(ExasimSolver& solver, int argc, char** argv,
                           MPI_Comm comm)
{
    int err = solver.InitializeEnvironment(argc, argv, comm);
    if (err)
        return err;

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
    if (err)
        return err;

    err = solver.Solve();
    if (err) {
        solver.Finalize();
        return err;
    }

    return solver.Finalize();
}

inline int InitializeExasimSolver(ExasimSolver& solver, int argc, char** argv,
                                  MPI_Comm comm)
{
    int err = solver.InitializeEnvironment(argc, argv, comm);
    if (err)
        return err;

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
    if (err)
        return err;

    return 0;
}

inline int InitializeExasimSolver(ExasimSolver& solver, int argc, char** argv,
                                  MPI_Comm comm,
                                  const std::vector<int>& vecBuiltinModelID)
{
    int err = solver.InitializeEnvironment(argc, argv, comm);
    if (err)
        return err;

    err = solver.ParseInputs(argc, argv);
    if (err) {
        solver.Finalize();
        return err;
    }

    err = ConfigureModelDefinitions(solver, vecBuiltinModelID);
    if (err) {
        solver.Finalize();
        return err;
    }

    err = solver.InitializeModels();
    if (err)
        return err;

    return 0;
}

inline int RunExasimSolver(ExasimSolver& solver, int argc, char** argv,
                           MPI_Comm comm,
                           const std::vector<int>& vecBuiltinModelID)
{
    int err = solver.InitializeEnvironment(argc, argv, comm);
    if (err)
        return err;

    err = solver.ParseInputs(argc, argv);
    if (err) {
        solver.Finalize();
        return err;
    }

    err = ConfigureModelDefinitions(solver, vecBuiltinModelID);
    if (err) {
        solver.Finalize();
        return err;
    }

    err = solver.InitializeModels();
    if (err)
        return err;

    err = solver.Solve();
    if (err) {
        solver.Finalize();
        return err;
    }

    return solver.Finalize();
}
