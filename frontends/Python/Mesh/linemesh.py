from numpy import linspace, array

def linemesh(m):
    p = array([linspace(0,1,m)])
    t = array([[k,k+1] for k in range(0,m-1)])
    t = t.transpose()
    return p,t