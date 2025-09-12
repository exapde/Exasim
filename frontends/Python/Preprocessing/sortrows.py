from numpy import *

def sortrows(a):

  n = a.shape[1];

  if (n==2):
    ind = argsort(a.view('int64,int64'), order=['f0','f1'], axis=0);
  elif (n==3):
    ind = argsort(a.view('int64,int64,int64'), order=['f0','f1','f2'], axis=0);
  elif (n==4):
    ind = argsort(a.view('int64,int64,int64,int64'), order=['f0','f1','f2','f3'], axis=0);
  elif (n==5):
    ind = argsort(a.view('int64,int64,int64,int64,int64'), order=['f0','f1','f2','f3','f4'], axis=0);
  elif (n==6):
    ind = argsort(a.view('int64,int64,int64,int64,int64,int64'), order=['f0','f1','f2','f3','f4','f5'], axis=0);
  elif (n==7):
    ind = argsort(a.view('int64,int64,int64,int64,int64,int64,int64'), order=['f0','f1','f2','f3','f4','f5','f6'], axis=0);
  elif (n==8):
    ind = argsort(a.view('int64,int64,int64,int64,int64,int64,int64,int64'), order=['f0','f1','f2','f3','f4','f5','f6','f7'], axis=0);
  elif (n==9):
    ind = argsort(a.view('int64,int64,int64,int64,int64,int64,int64,int64,int64'), order=['f0','f1','f2','f3','f4','f5','f6','f7','f8'], axis=0);
  elif (n==10):
    ind = argsort(a.view('int64,int64,int64,int64,int64,int64,int64,int64,int64,int64'), order=['f0','f1','f2','f3','f4','f5','f6','f7','f8','f9'], axis=0);
  else:
      ind = argsort(a.view('int64'), axis=0);
      #ind = 0;

  ind = ind.flatten('F');
  b = a[ind,:];

  return b, ind;
