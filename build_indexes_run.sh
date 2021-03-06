#! /bin/bash

##wrapper to build container, option to run reference download and indexing, and run full NextFlow pipeline for analysis

##in a system using 'modules' to load tools, need to specify that singularity is in the user PATH, as sudo cannot see it otherwise
sudo -E env "PATH=$PATH" singularity build --writable bovine_DNA_RNA.simg Singularity

##genome indexes are not necessarily required by users as they may already exist on their local system; if genome indexes are required for BWA (~4GB) and STAR (~30GB), run the Nextflow scripts into a local reference direcetory which is then also supplied to the NextFlow main script
nextflow run umd3.1.create_ref_indexes.simg.nf \
            --refDir refs \
            -with-singularity bovine_DNA_RNA.simg \
            -c "bovine_DNA_RNA.nextflow.simg.config" \
            -with-report "ref.report.html" \
            -with-timeline "ref.timeline.html"

nextflow clean -f

##run NExtFlow to download data, align, call variants
