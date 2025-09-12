import os
import subprocess

def compilepdemodel(pde):
    # Save current working directory
    cdir = os.getcwd()

    # Construct path to the build directory
    cmp = os.path.join(pde['exasimpath'], "backend", "Model", "build")

    # Change to build directory
    os.chdir(cmp)

    # Run cmake and make
    subprocess.run(["cmake", ".."], check=True)
    subprocess.run(["make"], check=True)

    # Return to original directory
    os.chdir(cdir)

    return cmp
