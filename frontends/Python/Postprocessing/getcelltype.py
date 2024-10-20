import os

def getcelltype(nd,elemtype):

    if nd==2:
        if elemtype==0:
            cell_t = 5;
        else:
            cell_t = 9;
    elif nd==3:
        if elemtype==0:
            cell_t = 10;
        else:
            cell_t = 12;
    else:
        sys.exit("Number of mesh spatial dimensions should be 2 or 3");

    return cell_t
