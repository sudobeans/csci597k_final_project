import subprocess
import os
from pathlib import Path

output = "fakes"

files = subprocess.check_output("find . -type f", shell=True)

for path in files.decode().splitlines():
    # destination file and directory
    directory = os.path.dirname(path)
    directory = os.path.join(output, directory)
    file = os.path.basename(path)
    
 
    # Make path
    Path(directory).mkdir(exist_ok=True, parents=True)

    # Make file
    open(os.path.join(directory, file), "w").close()
