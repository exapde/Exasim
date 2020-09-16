from sys import platform, exit
import shutil
from time import sleep
import subprocess

def findexec(filename, version):

    status = shutil.which(filename);
    if status != None:
        return filename;

    status = shutil.which("/usr/bin/" + filename);
    if status != None:
        filename = "/usr/bin/" + filename;
        return filename;

    status = shutil.which("/usr/local/bin/" + filename);
    if status != None:
        filename = "/usr/local/bin/" + filename;
        return filename;

    status = shutil.which("/opt/local/bin/" + filename);
    if status != None:
        filename = "/opt/local/bin/" + filename;
        return filename;

    if status == None:
        print("Exasim could not find " + filename + " in /usr/bin, /usr/local/bin, /opt/local/bin\n");
        if platform == "darwin":
            print("Exasim tries to find " + filename + " in /Applications. It may take a while.\n");
            result = subprocess.run(['find', '/Applications', '-name', filename, '-type', 'f'], stdout=subprocess.PIPE);
            a = result.stdout.decode('utf-8');
            if len(a)>0:
                ii = a.find("/MacOS/" + filename);
                newfilename = a[0:ii] + "/MacOS/" + filename;
                print("Exasim found " + filename + " at " + newfilename + "\n");
                print("Please open initializepde.py in the folder Exasim/" + version + "/Python/Preprocessing\n");
                qstr1 = """ " """ + filename + """ " """;
                qstr1 = qstr1.replace(" ", "");
                qstr2 = """ " """ + newfilename + """ " """;
                qstr2 = qstr2.replace(" ", "");
                print("Replace pde." + filename + " =" + qstr1 + " with pde." + filename + " =" + qstr2 + "\n");
                print("Doing so will prevent Exasim from searching " + filename + " again to save time.\n");
                print("Read the above instructions carefully and press any key to continue...\n");
                sleep(20);
                filename = newfilename;
                return filename;
        # elseif platform == "linux" or platform == "linux2":
        # elseif platform == "win32":

            mystr = "Exasim could not find " + filename + " on your system.\n";
            mystr = mystr + "Please see the documentation to install " + filename + ".\n"
            mystr = mystr + "After installing " + filename + ", open initializepde.py in the folder Exasim/" + version + "/Python/Preprocessing\n";
            qstr1 = """ " """ + filename + """ " """;
            qstr1 = qstr1.replace(" ", "");
            newfilename = "/path/to/executable/" + filename;
            qstr2 = """ " """ + newfilename + """ " """;
            qstr2 = qstr2.replace(" ", "");
            mystr = mystr + "and replace pde." + filename + " = " + qstr1 + " with pde." + filename + " = " + qstr2 + "\n";
            print(mystr);
            exit("Exit Exasim");
