/*
* arrayProduct.cpp - example in MATLAB External Interfaces
*
* Multiplies an input scalar (np)
* times a MxN matrix (inMatrix)
* and outputs a MxN matrix (outMatrix)
*
* Usage : from MATLAB
*         >> outMatrix = arrayProduct(np, inMatrix)
*
* This is a C++ MEX-file for MATLAB.
* Copyright 2017 The MathWorks, Inc.
*
*/

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "mutation++.h"
#include "XYZ_lib.h"
#include "MatlabDataArray.hpp"


class MexFunction : public matlab::mex::Function {
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        using namespace matlab::data;
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        ArrayFactory factory;

        checkArguments(outputs, inputs);
        int np = inputs[0][0];
        int nd = 2;
        double vareps = inputs[1][0];
        const TypedArray<double> inArray = std::move(inputs[1]);

        Mutation::MixtureOptions opts("air_5");
        opts.setStateModel("ChemNonEq1T");
        opts.setThermodynamicDatabase("NASA-9");
        Mutation::Mixture mix(opts);
        int ns = mix.nSpecies();

        double T = 300.0;
        double rhoe;
        
        Array omega_i_array = factory.createArray<double>({ np,ns });
        Array domega_i_drho_j_array = factory.createArray<double>({ np,ns,ns });
        Array domega_i_dT_array = factory.createArray<double>({ np,ns });

        Array dT_drhoi_array = factory.createArray<double>({ np,ns });
        Array dT_drhoe_array = factory.createArray<double>({ np,1 });

        // create working arrays 
        // TODO: there must be a better way than this. 
        double wdot[5];
        double wdot_eps[5]; // used for dw_dt
        double rho_i[5];
        double dwdot_drhoi[25];
        double dTdU[6];
        double Uwork[5];
        // double vareps = 1.0e-6;
        double denom;
        for (int i = 0; i < np; i++)
        {
            for (int j = 0; j < ns; j++)
            {
                rho_i[j] = inputs[2][i][j];
            }
            rhoe = inputs[2][i][ns];
            mix.setState(rho_i, &rhoe, 0); 

            // Evaluate source term w_i
            mix.netProductionRates(wdot);

            // Evaluate dw_i/drho_j
            mix.jacobianRho(dwdot_drhoi);

            // Evaluate dT/drho_i 
            denom = mix.density() * mix.mixtureFrozenCvMass();
            mix.getEnergiesMass(Uwork);
            for (int j = 0; j<ns; j++)
            {
                dTdU[j] = -Uwork[j] / denom;
            }
            // Evaluate dT/drhoe;
            dTdU[ns] = 1.0 / denom;

            // Pack w_i, dw/drho_i, dT/drho_i, dT/drhoe into output arrays
            dT_drhoe_array[i] = dTdU[ns];
            for (int j = 0; j < ns; j++)
            {
                omega_i_array[i][j] = wdot[j];
                dT_drhoi_array[i][j] = dTdU[j];
                for (int k = 0; k < ns; k++)
                {
                    domega_i_drho_j_array[i][j][k] = dwdot_drhoi[j*ns + k];
                }
            }

            // For dw_i/dt, use finite differences
            double Teps = mix.T() + vareps;
            mix.setState(rho_i, &Teps, 1);
            mix.netProductionRates(wdot_eps);
            for (int j = 0; j < ns; j++)
            {
                domega_i_dT_array[i][j] = (wdot_eps[j] - wdot[j])/vareps;
            }
        }

        outputs[0] = omega_i_array;
        outputs[1] = domega_i_drho_j_array;
        outputs[2] = domega_i_dT_array;
        outputs[3] = dT_drhoi_array;
        outputs[4] = dT_drhoe_array;
        // we also need dT/dre;
    }

// Okay for some reason having a hard time defining member functions...
    // void dT_dUstate(double* dTdU, double* Uwork, int ns, int ndim, Mutation::Mixture mix)
    // {
    // // Uwork is a temporary array for storing things
    // // Ustate = (rho_i, rhoe, rhoev)
    // // Works entirely in dimensional variables, must be nondimensionalized after the fact. 
    // double denom = mix.density() * mix.mixtureFrozenCvMass();
    // mix.getEnergiesMass(Uwork);
    // for (int j = 0; j<ns; j++)
    // {
    //   dTdU[j] = -Uwork[j] / denom;
    // }
    // dTdU[ns] = 1.0 / denom;
    // }


    void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        matlab::data::ArrayFactory factory;

        if (inputs.size() != 3) {
            matlabPtr->feval(u"error", 
                0, std::vector<matlab::data::Array>({ factory.createScalar("Three inputs required") }));
        }

        if (inputs[0].getNumberOfElements() != 1) {
            matlabPtr->feval(u"error", 
                0, std::vector<matlab::data::Array>({ factory.createScalar("First input (np) must be a scalar") }));
        }

        if (inputs[1].getNumberOfElements() != 1) {
            matlabPtr->feval(u"error", 
                0, std::vector<matlab::data::Array>({ factory.createScalar("Second input (nd) must be a scalar") }));
        }
        
        if (inputs[2].getType() != matlab::data::ArrayType::DOUBLE ||
            inputs[2].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE) {
            matlabPtr->feval(u"error", 
                0, std::vector<matlab::data::Array>({ factory.createScalar("Input matrix must be type double") }));
        }

        if (inputs[2].getDimensions().size() != 2) {
            matlabPtr->feval(u"error", 
                0, std::vector<matlab::data::Array>({ factory.createScalar("Input must be m-by-n dimension") }));
        }
    }
};