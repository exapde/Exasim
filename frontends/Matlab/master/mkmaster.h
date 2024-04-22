#ifndef __MKMASTER_H
#define __MKMASTER_H

#include "wrapper.h"

struct masterstruct {
    Int nd;
    Int porder;
    Int pgauss;
    Int npv;
    Int npf;
    Int ngv;
    Int ngf;
    Int elemtype;
    Int nodetype;
    
    vector<double> plocvl;
    vector<double> gpvl;
    vector<double> gwvl;
    vector<double> plocfc;
    vector<double> gpfc;
    vector<double> gwfc;
    vector<double> shapvl;
    vector<double> shapvt;
    vector<double> shapvg;
    vector<double> shapvgdotshapvl;
    vector<double> shapfc;
    vector<double> shapft;
    vector<double> shapfg;
    vector<double> shapfgdotshapfc;
    vector<double> shapmv;
    vector<double> shapmf;
    vector<Int> perm;
    
    vector<double> shapnv;
    vector<double> shapnvt;
    
    vector<double> projLowP;
    
    
    // NEW FIELDS FOR NEW MATRIX ASSEMBLY:
    Int pgaussR;
    Int pgaussJ;
    Int pgaussQ;
    
    Int quadTypeR;
    Int quadTypeJ;
    Int quadTypeQ;
    
    Int nqvR;
    Int nqvQ;
    Int nqvJ;
    vector<double> gpvlR;
    vector<double> gpvlJ;
    vector<double> gpvlQ;
    vector<double> gwvlR;
    vector<double> gwvlJ;
    vector<double> gwvlQ;
    vector<double> shapvlR;
    vector<double> shapvlJ;
    vector<double> shapvlQ;
    vector<double> shapvtR;
    vector<double> shapvtJ;
    vector<double> shapvtQ;
    vector<double> shapvgR;
    vector<double> shapvgJ;
    vector<double> shapvgQ;
    vector<double> shapvgdotshapvlR;
    vector<double> shapvgdotshapvlJ;
    vector<double> shapvgdotshapvlQ;
    
    Int nqfR;
    Int nqfQ;
    Int nqfJ;
    vector<double> gpfcR;
    vector<double> gpfcJ;
    vector<double> gpfcQ;
    vector<double> gwfcR;
    vector<double> gwfcJ;
    vector<double> gwfcQ;
    vector<double> shapfcR;
    vector<double> shapfcJ;
    vector<double> shapfcQ;
    vector<double> shapftR;
    vector<double> shapftJ;
    vector<double> shapftQ;
    vector<double> shapfgR;
    vector<double> shapfgJ;
    vector<double> shapfgQ;
    vector<double> shapfgdotshapfcR;
    vector<double> shapfgdotshapfcJ;
    vector<double> shapfgdotshapfcQ;
};

void mkmaster(masterstruct &master, Int porder, Int *pgauss, Int *quadType, Int nd, Int elemtype, Int nodetype);

#endif
