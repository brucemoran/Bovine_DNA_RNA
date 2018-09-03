#! /bin/bash

##wrapper to build container, and index references if required
##in a system using 'modules' to load tools, need to specify that singularity is in the user PATH, as sudo cannot se it otherwise
sudo -E env "PATH=$PATH" singularity build --writable <name_the_container_file.simg> Singularity

##genome indexes are not necessarily required by users as they may already exist on their local system; if genome indexes are required for BWA (~4GB) and STAR (~30GB), run the Nextflow scripts in /data
