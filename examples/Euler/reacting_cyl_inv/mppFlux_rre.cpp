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
        int nd = inputs[1][0];
        const TypedArray<double> inArray = std::move(inputs[1]);

        // std::cout << "Initializing mixture..." <<std::endl;

        Mutation::MixtureOptions opts("air_5");
        opts.setStateModel("ChemNonEq1T");
        opts.setThermodynamicDatabase("NASA-9");
        Mutation::Mixture mix(opts);
        // mix.addComposition("N:0.8, O:0.2", true);
        int ns = mix.nSpecies();

        // std::cout << "Done initializing..." <<std::endl;

        // Setup the default composition
        double T = 300.0;
        double RU = 8.314471468617452;
        // double P = Mutation::ONEATM;
        // mix.setState(&T,&P,1);

        // std::cout << "The mixture temp is: " <<mix.T() <<std::endl;
        
        Array P_array = factory.createArray<double>({ np,1 });
        Array T_array = factory.createArray<double>({ np,1 });
        Array dT_drhoi_array = factory.createArray<double>({ np,ns });
        Array dT_drhoe_array = factory.createArray<double>({ np,1 });
        // Array omega_i_array = factory.createArray<double>({ np,ns });
        // Array X_array = factory.createArray<double>({ np,ns });
        // Array Y_array = factory.createArray<double>({ np,ns });

        // create working arrays 
        // TODO: there must be a better way than this. 
        double wdot[5];
        double X[5];
        double Y[5];
        double rho_i[5];
        double rhoe;
        double dTdU[6];
        double Uwork[5];
        double denom;

        for (int i = 0; i < np; i++)
        {
            for (int j = 0; j < ns; j++)
            {
                rho_i[j] = inputs[2][i][j];
            }
            rhoe = inputs[2][i][ns];
            mix.setState(rho_i, &rhoe, 0);

            P_array[i] = mix.P();
            T_array[i] = mix.T();

            // Evaluate dT/drho_i 
            denom = mix.density() * mix.mixtureFrozenCvMass();
            mix.getEnergiesMass(Uwork);
            for (int j = 0; j<ns; j++)
            {
                dTdU[j] = -Uwork[j] / denom;
            }
            // Evaluate dT/drhoe;
            dTdU[ns] = 1.0 / denom;

            // dT/drho_i, dT/drhoe into output arrays
            dT_drhoe_array[i] = dTdU[ns];
            for (int j = 0; j < ns; j++)
            {
                dT_drhoi_array[i][j] = dTdU[j];
            }

        }

        outputs[0] = P_array;
        outputs[1] = T_array;
        outputs[2] = dT_drhoi_array;
        outputs[3] = dT_drhoe_array;
    }

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