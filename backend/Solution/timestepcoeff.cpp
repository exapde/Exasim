/*
    timestepcoeff.cpp

    This file provides functions to compute coefficients for time-stepping schemes
    used in numerical integration, specifically Diagonally Implicit Runge-Kutta (DIRK)
    and Backward Differentiation Formula (BDF) methods.

    Functions:

    - void DIRKcoeff(dstype * c, dstype * d, dstype * t, Int nstage, Int torder)
        Computes the coefficients for DIRK schemes based on the number of stages (nstage)
        and the temporal order (torder). The coefficients are stored in arrays c, d, and t.
        Supported DIRK configurations:
            - DIRK(1,1), DIRK(1,2), DIRK(2,2), DIRK(2,3), DIRK(3,3), DIRK(3,4)
        If an unsupported configuration is requested, an error is raised.

    - void BDFcoeff(dstype * c, dstype * t, Int torder)
        Computes the coefficients for BDF schemes based on the temporal order (torder).
        The coefficients are stored in arrays c and t.
        Supported BDF orders: 1, 2, 3
        If an unsupported order is requested, an error is raised.

    - void TimestepCoefficents(commonstruct common)
        Determines the temporal scheme (DIRK or BDF) from the 'common' structure and
        computes the corresponding coefficients by calling DIRKcoeff or BDFcoeff.
        For BDF schemes, the number of stages is set to 1.

    Notes:
    - The functions assume that the input arrays have sufficient size for the requested scheme.
    - Error handling is performed for unsupported configurations.
*/
#ifndef __TIMESTEPCOEFF
#define __TIMESTEPCOEFF

void DIRKcoeff(dstype * c, dstype * d, dstype * t, Int nstage, Int torder) 
{
    if (nstage == 1 && torder == 1) {   //DIRK(1,1)
        d[0] = 1.0;
        c[0] = 1.0;
        t[0] = 1.0;
    }
    else if (nstage == 1 && torder == 2) {   //DIRK(1,2)
        d[0] = 2.0;
        c[0] = 2.0;
        t[0] = 0.5;
    }
    else if (nstage == 2 && torder == 2) {   //DIRK(2,2)           
        d[0] = 3.414213562373096;
        d[1] = -8.242640687119289;
        d[2] = 0.0;
        d[3] = 3.414213562373096;
        c[0] = 0.0;
        c[1] = 1.0;
        t[0] = 0.292893218813452;
        t[1] = 1.0;
    }
    else if (nstage == 2 && torder == 3) {   //DIRK(2,3)
        d[0] = 1.267949192431123;
        d[1] = 0.928203230275509;
        d[2] = 0.0;
        d[3] = 1.267949192431123;
        c[0] = 1.098076211353316;
        c[1] = 0.633974596215561;
        t[0] = 0.788675134594813;
        t[1] = 0.211324865405187;
    }
    else if (nstage == 3 && torder == 3) {   //DIRK(3,3)
        d[0] = 2.294280360323567;
        d[1] = -1.484721005721436;
        d[2] = -8.556127801733982;
        d[3] = 0.0;
        d[4] = 2.294280360323568;
        d[5] = 3.391748836909791;
        d[6] = 0.0;
        d[7] = 0.0;
        d[8] = 2.294280360323567;
        c[0] = 0.0;
        c[1] = 0.0;
        c[2] = 1.0;
        t[0] = 0.435866521500000;
        t[1] = 0.717933260750000;
        t[2] = 1.000000000000000;
    }
    else if (nstage == 3 && torder == 4) {   //DIRK(3,4)
        d[0] = 0.935822227524088;
        d[1] = 0.497940606760015;
        d[2] = -0.345865915800969;
        d[3] = 0.0;
        d[4] = 0.935822227524088;
        d[5] = 2.867525668568206;
        d[6] = 0.0;
        d[7] = 0.0;
        d[8] = 0.935822227524088;
        c[0] = 0.445622407287714;
        c[1] = 1.064177772475912;
        c[2] = 0.120614758428183;
        t[0] = 1.068579021301629;
        t[1] = 0.500000000000000;
        t[2] = -0.068579021301629;
    }
    else {        
        error("DIRK is not implemented yet \n");
    }                
}

void BDFcoeff(dstype * c, dstype * t, Int torder) 
{
    if (torder == 1) {   //BDF1        
        c[0] = 1.0;
        c[1] = -1.0;        
        t[0] = 1.0;
    }
    else if (torder == 2) {   //BDF2        
        c[0] = 1.5;
        c[1] = -2.0;
        c[2] = 0.5;                
        t[0] = 1.0;
    }
    else if (torder == 3) {   //BDF3        
        c[0] = 11.0/6.0;
        c[1] = -3.0;
        c[2] = 1.5;
        c[3] = -1.0/3.0;                
        t[0] = 1.0;
    }
    else {
        printf("BDF%d  not implemented yet.", torder);
        error("\n");
    }                    
}

void TimestepCoefficents(commonstruct common) 
{    
    if (common.temporalScheme==0)  // DIRK 
    {
        DIRKcoeff(common.DIRKcoeff_c, common.DIRKcoeff_d, common.DIRKcoeff_t, common.tstages, common.torder);     
    }
    else // BDF 
    {
        BDFcoeff(common.BDFcoeff_c, common.BDFcoeff_t, common.torder); 
        common.tstages = 1;
    }        
}

#endif
