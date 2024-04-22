#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>

using namespace std;

#include "mkmaster.h"
typedef int Int;

int main(int argc, char** argv) 
{
    Int porder = 2;
    Int pgauss[] = {3*porder,3*porder,3*porder};
    Int quadType[] = {0,0,0};
    Int nd = 2;
    Int elemtype = 0;
    Int nodetype = 0;

    masterstruct master;

    mkmaster(master, porder, &pgauss[0], &quadType[0], nd, elemtype, nodetype);

    return 0;
}
