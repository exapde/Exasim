#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>

using namespace std;

#include "mkmaster.h"
typedef int Int;

void writedarray(ofstream &out, vector<double> &a, Int N)
{
    out.write( reinterpret_cast<char*>( &a[0] ), sizeof(double) * N );
}

template <typename T> string NumberToString ( T Number )
{
    ostringstream ss;
    ss << Number;
    return ss.str();
}

int main(int argc, char** argv) 
{
    Int quadType[] = {0,0,0};
    Int nodetype = 0;
    masterstruct master;
    
    for (Int porder=0; porder <= 6; porder++)
        for (Int quadOrder=1; quadOrder <= 3; quadOrder++)
            for (Int elemtype=0; elemtype <= 1; elemtype++)
                for (Int nd=1; nd <= 3; nd++) {
                    if (quadOrder*porder <= 16) {
                        Int pgauss[] = {quadOrder*porder,quadOrder*porder,quadOrder*porder};
                        
                        mkmaster(master, porder, &pgauss[0], &quadType[0], nd, elemtype, nodetype);
                        
                        string filename = "P" + NumberToString(porder) + "D" + NumberToString(nd) +"E" + NumberToString(elemtype) + "Q" + NumberToString(quadOrder) + ".bin";
                        ofstream out(filename.c_str(), ios::out | ios::binary);
                        
                        if (!out)
                            exit(-1);

                        writedarray(out, master.plocvl, (Int) master.plocvl.size());
                        writedarray(out, master.plocfc, (Int) master.plocfc.size());
                        writedarray(out, master.gpvlR, (Int) master.gpvlR.size());
                        writedarray(out, master.gwvlR, (Int) master.gwvlR.size());
                        writedarray(out, master.gpfcR, (Int) master.gpfcR.size());
                        writedarray(out, master.gwfcR, (Int) master.gwfcR.size());
                        writedarray(out, master.shapvlR, (Int) master.shapvlR.size());
                        writedarray(out, master.shapfcR, (Int) master.shapfcR.size());

                        out.close();
                    }
                }
    return 0;
}
