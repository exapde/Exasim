void nondimensionalizeConsVars(double* Ucons, double* uinf, int nspecies, int nd)
{
  // Modifies Ucons vector in place to have nondimensional quantities
  double rho_scale = uinf[0];
  double u_scale = uinf[1];
  double rhoe_scale = uinf[2];
  for (int i=0; i<nspecies; i++)
  {
    Ucons[i] = Ucons[i] / rho_scale;
  }
  for (int i=0; i<nd; i++)
  {
    Ucons[nspecies+i] = Ucons[nspecies+i] / (rho_scale * u_scale);
  }
  Ucons[nspecies+nd] = Ucons[nspecies+nd] / (rhoe_scale);
}

void dimensionalizeConsVars(double* Ucons, double* uinf, int nspecies, int nd)
{
  // Modifies Ucons vector in place to have dimensional quantities
  double rho_scale = uinf[0];
  double u_scale = uinf[1];
  double rhoe_scale = uinf[2];
  for (int i=0; i<nspecies; i++)
  {
    Ucons[i] = Ucons[i] * rho_scale;
  }
  for (int i=0; i<nd; i++)
  {
    Ucons[nspecies+i] = Ucons[nspecies+i] * (rho_scale * u_scale);
  }
  Ucons[nspecies+nd] = Ucons[nspecies+nd] * (rhoe_scale);
}

void nondimensionalizeStateVars(double* Ustate, double* uinf, int nspecies)
{
  // Modifies Ucons vector in place to have nondimensional quantities
  double rho_scale = uinf[0];
  double u_scale = uinf[1];
  double rhoe_scale = uinf[2];
  for (int i=0; i<nspecies; i++)
  {
    Ustate[i] = Ustate[i] / rho_scale;
  }
  Ustate[nspecies] = Ustate[nspecies] / (rhoe_scale);
  // TODO: allow multiple energy equations 
}

void dimensionalizeStateVars(double* Ustate, double* uinf, int nspecies)
{
  // Modifies Ucons vector in place to have nondimensional quantities
  double rho_scale = uinf[0];
  double u_scale = uinf[1];
  double rhoe_scale = uinf[2];
  for (int i=0; i<nspecies; i++)
  {
    Ustate[i] = Ustate[i] * rho_scale;
  }
  Ustate[nspecies] = Ustate[nspecies] * (rhoe_scale);
  // TODO: allow multiple energy equations 
}

void conservativeToState(double* Ucons, double* Ustate, double* uinf, int nspecies)
{   
    // Maps from conservative fluid variables (r, ru, rE, rEv) to state variables (r, re, rev)
    //Ucons[nspecies+2]: state variables rho_i, rhou, rhoE (dimensional)
    //Ustate[nspecies+1]: inputs for mutation: rho_i, rhoe (dimensional)
    double rho = 0.0;
    for (int i=0; i<nspecies; i++)
    {
      Ustate[i] = Ucons[i]; //rho_i
      rho = rho + Ucons[i];
    }
    double rhou = Ucons[nspecies];
    double rhoE = Ucons[nspecies+1];
    double u = rhou/rho;
    Ustate[nspecies] = rhoE - 0.5 * rho * u * u; //rhoe
}

void stateToConsVars(double* Ucons, double u, double* Ustate, double* uinf, int nspecies)
{   
    // Maps from state variables (r, re) and velocity (u) to cons vars (r, ru, rE)
    //Ucons[nspecies+2]: state variables rho_i, rhou, rhoE (dimensional)
    //Ustate[nspecies+1]: inputs for mutation: rho_i, rhoe (dimensional)
    // TODO: currently for 1D
    double rho = 0.0;
    for (int i=0; i<nspecies; i++)
    {
      Ucons[i] = Ustate[i]; //rho_i
      rho = rho + Ucons[i];
    }
    double rhou = rho * u; 
    double rhoe = Ustate[nspecies];
    Ucons[nspecies] = rhou; // TODO: should be Ucons[nspecies + 1:nd]
    Ucons[nspecies+1] = rhoe + 0.5 * rho * u * u;
}


// void uinflow(double* Ucons, double p, double *param, double* uinf, Mutation::Mixture *mix)
// { // Subsonic inflow: use pressure from solution
//   //                  and given Tinf, Yinf, uinf to solve for rho, rhou, rhoE
//     int nspecies = 5;
//     // Scaling information
//     double rho_inf = uinf[0];
// 		double uv_inf = uinf[1];
// 		double rhoe_inf = uinf[2];

//     // Known values from prescribed inflow
// 		double T_inf = param[3*nspecies+ 5];
//     double Yinf[5];
//     double uv_inflow = param[3*nspecies+4]; //inflow; here same as scaling but not necessarily true
//     for (int ispecies = 0; ispecies<nspecies; ispecies++)
//     {
//       Yinf[ispecies] = param[2*nspecies + 4 + ispecies];
//     }

//     // Known values from solution (nondim pressure)
//     double pdim = p * rhoe_inf;
//     double pt_arr[2] = {pdim, T_inf};

//     // Update state with Y_i, P, T
//     mix->setState(Yinf, pt_arr, 2);

//     // Get conservative quantities  
//     double rho_i[5]; // memory for densities
//     mix->densities(rho_i);
//     double rho = mix->density();
//     double rhoe = mix->mixtureEnergyMass() * rho;
//     double Ustate[6] = {rho_i[0],rho_i[1],rho_i[2],rho_i[3],rho_i[4],rhoe};

//     stateToConsVars(Ucons, uv_inflow, Ustate, uinf, nspecies);
//     nondimensionalizeConsVars(Ucons, uinf, nspecies, 1);
// }

// void uoutflow(double* Ucons, double* param, double* uinf, Mutation::Mixture *mix)
// {
//     // Subsonic outflow: use rho_i, rhou from solution
//     //                   given Pout, solve for energy rhoE
//     int nspecies = 5;
//     // scaling information
//     double rho_inf = uinf[0];
// 		double uv_inf = uinf[1];
// 		double rhoe_inf = uinf[2];
//     double RU = Mutation::RU;

//     double Ustate[6];
//     double rho_i[5]; // If kept around, these arrays should be initialized as 
//                     //  some member of temp (maybe tmpapp array so it's separate from tmpg?) 
//     double rho = 0.0;

//     dimensionalizeConsVars(Ucons, uinf, nspecies, 1);
//     conservativeToState(Ucons, Ustate, uinf, nspecies);

//     // Pressure outlet
//     double P_out = param[3*nspecies+6];

//     // Use pressure and rho from solution to get Temperature
//     double denom = 0.0;
//     for (int i = 0; i < nspecies; i++)
//     {
//       rho_i[i] = Ustate[i];
//       rho = rho + rho_i[i];
//       // Ucons[i] = rho_i[i];
//       denom = denom + RU * Ustate[i]/mix->speciesMw(i);
//     }
//     double T_out = P_out / denom;
//     double u_out = Ucons[nspecies]/rho;
//     // double u_out = 370.0;
//     // printf("Tout: %f\n", T_out);
//     // Set state with desnities and temperature
//     mix->setState(rho_i, &T_out, 1);

//     // Update energy
//     double rhoe = mix->mixtureEnergyMass() * rho;
//     Ustate[nspecies] = rhoe;

//     // printf("rhou before update: %f\n", Ucons[5]);
//     stateToConsVars(Ucons, u_out, Ustate, uinf, nspecies);
//     //  printf("rhou after update: %f\n", Ucons[5]);  
//     nondimensionalizeConsVars(Ucons, uinf, nspecies, 1);
   
//     // Here I can check that Ustate is the same mix.densities? 
//     // and that rho u doesn't change ? ? ? 
//     // double rhotest[5];
//     // mix->densities(rhotest);
//     // for (int i = 0; i< nspecies; i++)
//     // {
//       // printf("rho mix %i: %f\n", i, rhotest[i]/rho_inf);
//       // printf("rho state %i: %f\n", i, Ucons[i]);
//     // }
// }

void uinflow(double* Ucons, double p, double T_inf, double uv_inflow, double* Yinf, double* scales, Mutation::Mixture *mix)
{ // Subsonic inflow: use pressure from solution
  //                  and given Tinf, Yinf, uinf to solve for rho, rhou, rhoE
    int nspecies = 5;
    // Scaling information
    double rho_inf = scales[0];
		double uv_inf = scales[1];
		double rhoe_inf = scales[2];

    // Known values from prescribed inflow
		// double T_inf = datainflow[0];
    // double Yinf[5];
    // double uv_inflow = datainflow[1]; //inflow; here same as scaling but not necessarily true
    // for (int ispecies = 0; ispecies<nspecies; ispecies++)
    // {
    //   Yinf[ispecies] = param[2*nspecies + 4 + ispecies];
    // }
    // double *Yinf = datainflow[2];

    // Known values from solution (nondim pressure)
    double pdim = p * rhoe_inf;
    double pt_arr[2] = {pdim, T_inf};

    // Update state with Y_i, P, T
    mix->setState(Yinf, pt_arr, 2);

    // Get conservative quantities  
    double rho_i[5]; // memory for densities
    mix->densities(rho_i);
    double rho = mix->density();
    double rhoe = mix->mixtureEnergyMass() * rho;
    double Ustate[6] = {rho_i[0],rho_i[1],rho_i[2],rho_i[3],rho_i[4],rhoe};

    stateToConsVars(Ucons, uv_inflow, Ustate, scales, nspecies);
    nondimensionalizeConsVars(Ucons, scales, nspecies, 1);
}

void uoutflow(double* Ucons, double outP, double* scales, Mutation::Mixture *mix)
{
    // Subsonic outflow: use rho_i, rhou from solution
    //                   given Pout, solve for energy rhoE
    int nspecies = 5;
    // scaling information
    double rho_inf = scales[0];
		double uv_inf = scales[1];
		double rhoe_inf = scales[2];
    double RU = Mutation::RU;

    double Ustate[6];
    double rho_i[5]; // If kept around, these arrays should be initialized as 
                    //  some member of temp (maybe tmpapp array so it's separate from tmpg?) 
    double rho = 0.0;

    dimensionalizeConsVars(Ucons, scales, nspecies, 1);
    conservativeToState(Ucons, Ustate, scales, nspecies);

    // Pressure outlet
    double P_out = outP;

    // Use pressure and rho from solution to get Temperature
    double denom = 0.0;
    for (int i = 0; i < nspecies; i++)
    {
      rho_i[i] = Ustate[i];
      rho = rho + rho_i[i];
      // Ucons[i] = rho_i[i];
      denom = denom + RU * Ustate[i]/mix->speciesMw(i);
    }
    double T_out = P_out / denom;
    // double T_out = 300.0;
    // double u_out = Ucons[nspecies]/rho;
    double u_out = 0.0;
    // printf("u_out: %f\n", u_out);
    // double u_out = 370.0;
    // printf("Tout: %f\n", T_out);
    // Set state with desnities and temperature
    mix->setState(rho_i, &T_out, 1);

    // Update energy
    double rhoe = mix->mixtureEnergyMass() * rho;
    Ustate[nspecies] = rhoe;

    // printf("rhou before update: %f\n", Ucons[5]);
    stateToConsVars(Ucons, u_out, Ustate, scales, nspecies);
    //  printf("rhou after update: %f\n", Ucons[5]);  
    nondimensionalizeConsVars(Ucons, scales, nspecies, 1);
}

// void getdTdU(double* dTdU, double* Ucons, double* Uwork, Mutation::Mixture *mix, int nspecies, int ndim)
// {
//   // Uwork is a temporary array for storing things
//     double rho = 0.0;
//     for (int i=0; i<nspecies; i++)
//     {
//       Uwork[i] = Ucons[i]; //rho_i
//       rho = rho + Ucons[i];
//     }
//     for (int j = 0; j<ndim; j++)
//     {
//       Uwork[nspecies+j] = Ucons[nspecies+j] / rho; // u_i
//     }
//     Uwork[nspecies+ndim] = mix->mixtureEnergyMass() // OKAY SHOULD I DO THIS ALL IN DIMENSIONAL OR NONDIM COORDS...what did I do last time? 
//     double rhou = Ucons[nspecies];
//     double rhoE = Ucons[nspecies+1];

// }

// void getdPdU()
// {

// }

// double getdTdr(double* dTdr_i, double* Ucons, Mutation::Mixture *mix, int nspecies, int ndim)
// {
//     double rho = 0.0;
//     for (int i=0; i<nspecies; i++)
//     {
//       Ustate[i] = Ucons[i]; //rho_i
//       rho = rho + Ucons[i];
//     }
//     double rhou = Ucons[nspecies];
//     double rhoE = Ucons[nspecies+1];
//     // double e = rhoe/rho;
//     // double denom = (rho*mixture->mixtureEquilibriumCvMass());
//     // return -e/denom;
// }

// double getdTdr(Mutation::Mixture *mix)
// {
//     // double e = mixture->mixtureEnergyMass();
//     // double rho = mixture->density();
//     // double denom = (rho*mixture->mixtureEquilibriumCvMass());
//     // return -e/denom;
// }

// double dTdru(double* rho_i, double* rhou)

// double getdTdre(double rho, Mutation::Mixture *mix)
// {
//     // double Cv = mixture->mixtureEquilibriumCvMass();
//     // double denom = rho * Cv;
//     // return 1.0/denom;
// }

// double getdTdre(Mutation::Mixture *mix)
// {
// //     double Cv = mixture->mixtureEquilibriumCvMass();
// //     double rho = mixture->density();
// //     double denom = rho * Cv;
// //     return 1.0/denom;
// }

// double getP()
// {

// }

// double getKineticSource()
// {

// }