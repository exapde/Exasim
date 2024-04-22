from numpy import *
from sys import platform, exit
from findexec import findexec
import shutil
import csv
import os

def partition(t2f,ne,np,metis):

    # Generate a temporary file to be used in METIS
    print("Writing input files for METIS...");
    t2f = t2f.T;
    t2f = t2f.astype(int);
    with open('temp.txt','w') as f_handle:
        savetxt(f_handle, array([ne]).astype(int), delimiter=' ', fmt='%d')
    with open('temp.txt','a') as f_handle:
        savetxt(f_handle, t2f, delimiter=' ', fmt='%d')

    metisstatus0 = shutil.which(metis);
    metisstatus1 = shutil.which("mpmetis");
    metisstatus2 = shutil.which("/usr/bin/mpmetis");
    metisstatus3 = shutil.which("/usr/local/bin/mpmetis");
    metisstatus4 = shutil.which("/opt/local/bin/mpmetis");

    if metisstatus0 != None:
        metis = metis;
    elif metisstatus1 != None:
        metis = "mpmetis"
    elif metisstatus2 != None:
        metis = "/usr/bin/mpmetis";
    elif metisstatus3 != None:
        metis = "/usr/local/bin/mpmetis";
    elif metisstatus4 != None:
        metis = "/opt/local/bin/mpmetis";
    else:
        exit("Exasim search in /usr/bin, /usr/local/bin, and /opt/local/bin and could not find Metis. Please see the documentation to install it. After installation, please set its path to app.metis");

    print("Calling METIS and reading output files...");
    # call mpmetis
    nps = str(np);
    mystr = metis + " temp.txt " + nps;
    os.system(mystr);

    # get mesh partitioning data
    estr = "temp.txt.epart." + nps;
    #epart = genfromtxt(estr, delimiter=',')
    epart = genfromtxt(estr);

    # get node partitioning data
    nstr = "temp.txt.npart." + nps;
    #npart = genfromtxt(nstr, delimiter=',')
    npart = genfromtxt(nstr);

    # remove files
    os.remove("temp.txt");
    os.remove(estr);
    os.remove(nstr);

    return epart, npart



    # # current directory
    # cdir = os.getcwd();
    # ii = cdir.find("Exasim");
    #
    # if platform == "linux" or platform == "linux2":
    #     #cd([cdir "/metis/linux"]);
    #     metisdir = cdir[0:ii] + "Exasim/metis/linux";
    # elif platform == "darwin":
    #     #cd([cdir "/metis/mac"]);
    #     metisdir = cdir[0:ii] + "Exasim/metis/mac";
    # elif platform == "win32":
    #     metisdir = cdir[0:ii] + "Exasim/metis/windows";
    #
    # os.chdir(metisdir);
    #
    # # Generate a temporary file to be used in METIS
    # print("Writing input files for METIS...");
    #
    # t2f = t2f.T;
    # t2f = t2f.astype(int);
    # with open('temp.txt','w') as f_handle:
    #     savetxt(f_handle, array([ne]).astype(int), delimiter=' ', fmt='%d')
    # with open('temp.txt','a') as f_handle:
    #     savetxt(f_handle, t2f, delimiter=' ', fmt='%d')
    # # with open('temp.txt', 'w', newline='') as out:
    # #     writer = csv.writer(out, delimiter=',')
    # #     writer.writerow(ne)
    # # with open('temp.txt', 'a', newline='') as out:
    # #     writer = csv.writer(out, delimiter=',')
    # #     writer.writerow(t2f.transpose())
    #
    # print("Calling METIS and reading output files...");
    # # call mpmetis
    # nps = str(np);
    # mystr = "./mpmetis temp.txt " + nps;
    # os.system(mystr);
    #
    # # get mesh partitioning data
    # estr = "temp.txt.epart." + nps;
    # epart = genfromtxt(estr, delimiter=',')
    #
    # # get node partitioning data
    # nstr = "temp.txt.npart." + nps;
    # npart = genfromtxt(nstr, delimiter=',')
    #
    # # remove files
    # os.remove("temp.txt");
    # os.remove(estr);
    # os.remove(nstr);
    #
    # # move back to current directory
    # os.chdir(cdir);
    #
    # return epart, npart
