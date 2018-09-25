Bootstrap:docker
From:centos:centos7.4.1708

%help
    Container for Bovine DNA and RNA tools

%labels
    MAINTAINER Colin Sauze

%environment
    #environment variables defined in post section using
    ##SINGULARITY_ENVIRONMENT variable

%post
    #essential utilities
    yum -y install git wget bzip2 unzip which emacs

    #language and libraries
    yum -y install java-1.8.0-openjdk-devel gcc gcc-c++ glibc-devel make ncurses-devel zlib-devel libbzip2-devel bzip2-devel xz-devel perl-DBI lapack-devel atlas-devel freetype freetype-devel libpng-devel readline-devel pcre pcre-devel libcurl libcurl-devel

    #libclas and libatlas aren't put in the right places
    ln -s /usr/lib64/atlas/libtatlas.so /usr/lib64/libatlas.so
    ln -s /usr/lib64/atlas/libsatlas.so /usr/lib64/libcblas.so

    yum -y install https://centos7.iuscommunity.org/ius-release.rpm
    yum -y install python36u python36u-pip python36u-devel
    pip3.6 install --upgrade pip
    pip3.6 install --upgrade setuptools

    mkdir -p /usr/local/src
    cd /usr/local/src

    ##define env vars via S..._E... env var when in post
    ##see: https://www.sylabs.io/guides/2.5/user-guide/environment_and_metadata.html
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib' >> $SINGULARITY_ENVIRONMENT

    ##R
    #required libs
    yum -y install readline readline-devel pcre pcre-devel libcurl libcurl-devel

    #source
    wget https://cran.rstudio.com/src/base/R-3/R-3.5.1.tar.gz
    tar xf R-3.5.1.tar.gz
    cd R-3.5.1
    ./configure --with-x=no --prefix=/usr/local/
    make
    make install
    cd /usr/local/src

    #packages
    R --slave -e 'install.packages("BiocManager", repos="https://cloud.r-project.org/")'
    R --slave -e 'lapply(c("ggplot2","tidyverse","plyr","dplyr","biomaRt","reshape2","roots","Biobase","rgexf","fgsea","gtools","Rplinkseq","Rserve"),function(f){library("BiocManager"); BiocManager::install(f,update=TRUE,ask=FALSE)})'

    #multiqc
    pip3.6 install multiqc

    # ensembl vep
    # note built without mysql support, error message:
    # WARNING: DBD::mysql module not found. VEP can only run in offline (--offline) mode without DBD::mysql installed
    ##used --offline by preference
    yum install -y perl-CPAN perl-IO-Socket-SSL perl-Archive-Any perl-YAML perl-CPAN-Meta perl-Digest-MD5 perl-App-cpanminus perl-local-lib

    ##setting more that LANG locale is an issue for several tools
    ##https://github.com/CentOS/sig-cloud-instance-images/issues/71
    localedef -i en_US -f UTF-8 en_US.UTF-8
    echo -e "LANGUAGE="C"\nLC_ALL="C"" >> /etc/locale.conf
    echo 'export LANG=en_US.UTF-8' >> $SINGULARITY_ENVIRONMENT
    echo 'export LANGUAGE=C' >> $SINGULARITY_ENVIRONMENT
    echo 'export LC_ALL=C' >> $SINGULARITY_ENVIRONMENT

    #ensembl VEP
    #The VEP installer needs Ensembl which needs bioperl, Ensembl didn't have a cpan or rpm package, installing like its 1990!
    ##thank goodness it is such a wonderful tool=D
    cd /usr/local/lib
    wget ftp://ftp.ensembl.org/pub/ensembl-api.tar.gz
    wget https://cpan.metacpan.org/authors/id/C/CJ/CJFIELDS/BioPerl-1.6.924.tar.gz
    tar xf ensembl-api.tar.gz
    tar xf BioPerl-1.6.924.tar.gz


    echo 'export PERL5LIB=${PERL5LIB}:/usr/local/lib/ensembl/modules:/usr/local/lib/ensembl-io/modules:/usr/local/lib/ensembl-variation/modules:/usr/local/lib/ensembl-funcgen/modules:/usr/local/lib/BioPerl-1.6.924' >>$SINGULARITY_ENVIRONMENT

    #but the bioperl we just downloaded is missing the DB:HTS module, so use cpan to install it
    #this prompts us to install a load of extra bits of bioperl, tell it not to
    #yes n | cpan install Bio::DB::HTS
    #cpanm install --force Bio::DB::HTS #version mismatch in the tests requires forcing, actual program ok
    #cpanm DBI Set::IntervalTree PerlIO::gzip

    #cpan Module::Build

    cd /usr/local/src
    wget https://github.com/Ensembl/ensembl-vep/archive/release/92.5.tar.gz
    tar xf 92.5.tar.gz
    cd ensembl-vep-release-92.5

    echo "y" > commands
    echo "y" >> commands #agree to install cache files
    echo "10" >> commands #choose bos_taurus_merged_vep_92_UMD3.1.tar.gz
    echo "n" >> commands #overwrite cached files, is this needed -> no
    echo "n" >> commands #don't install FASTA files -> can use and specify in NextFlow scripts
    echo "n" >> commands #don't install plugins -> may need to revise

    mkdir /usr/local/src/ensembl-vep-release-92.5/cache
    perl ./INSTALL.pl --CACHEDIR /usr/local/src/ensembl-vep-release-92.5/cache < commands
    ln -s /usr/local/src/ensembl-vep-release-92.5/vep /usr/local/bin/vep
    cd /usr/local/src


    #samtools
    wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
    tar xf samtools-1.8.tar.bz2
    cd samtools-1.8
    ./configure --prefix=/usr/local/
    make
    make install
    cd /usr/local/src

    #bcftools
    wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2
    tar xf bcftools-1.8.tar.bz2
    cd bcftools-1.8
    ./configure --prefix=/usr/local
    make
    make install
    cd /usr/local/src

    #htslib
    wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2
    tar xf htslib-1.8.tar.bz2
    cd htslib-1.8
    ./configure --prefix=/usr/local
    make
    make install
    cd /usr/local/src

    #cramtools
    wget https://github.com/enasequence/cramtools/archive/v3.0.tar.gz
    tar xf v3.0.tar.gz
    cd cramtools-3.0/
    chmod a+x cramtools-3.0.jar
    mv cramtools-3.0.jar /usr/local/bin
    echo -e "#! /bin/bash\nexec java -jar /data/genome/reference/hg19/refGene/cramtools-3.0/cramtools-3.0.jar "$@"" > /usr/local/bin/cramtools
    chmod a+x /usr/local/bin/cramtools

    #bedtools
    wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz
    tar xf bedtools-2.27.1.tar.gz
    cd bedtools2
    make
    make install
    cd /usr/local/src

    #picard
    wget https://github.com/broadinstitute/picard/releases/download/2.18.9/picard.jar -O /usr/local/lib/picard.jar
    chmod a+x /usr/local/lib/picard.jar
    echo -e "#! /bin/bash\njavamem=""\nif [[ \$1 =~ "-Xmx" ]];then javamem=\$1; shift 1; fi\nexec java \$javamem -jar /usr/local/lib/picard.jar "\$@"" > /usr/local/bin/picard-tools
    chmod a+x /usr/local/bin/picard-tools

    #BWA
    wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
    tar xf bwa-0.7.17.tar.bz2
    cd bwa-0.7.17
    make
    mv bwa /usr/local/bin
    mv *.pl /usr/local/bin
    mv libbwa.a /usr/local/lib
    cd /usr/local/src

    #STAR
    wget https://github.com/alexdobin/STAR/archive/2.6.0c.tar.gz
    tar xf 2.6.0c.tar.gz
    cd STAR-2.6.0c/source
    make STAR
    mv STAR /usr/local/bin
    #make install
    cd /usr/local/src

    #kallisto
    wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz
    tar xf kallisto_linux-v0.44.0.tar.gz
    mv kallisto_linux-v0.44.0/kallisto /usr/local/bin
    cd /usr/local/src

    #fastqc
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
    unzip fastqc_v0.11.7.zip
    chmod a+x /usr/local/src/FastQC/fastqc
    ln -s /usr/local/src/FastQC/fastqc /usr/local/bin/fastqc
    cd /usr/local/src

    #BBMap
    wget "https://downloads.sourceforge.net/project/bbmap/BBMap_38.11.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbbmap%2Ffiles%2FBBMap_38.11.tar.gz%2Fdownload&ts=1531223392" -O BBMap_38.11.tar.gz
    tar xf BBMap_38.11.tar.gz
    cd bbmap
    ln -s /usr/local/src/bbmap/*.sh /usr/local/bin/
    cd /usr/local/src

    #sam2bam
    wget https://github.com/t-ogasawara/sam-to-bam/raw/master/build.sh
    bash build.sh

    #gatk
    wget https://github.com/broadinstitute/gatk/releases/download/4.0.6.0/gatk-4.0.6.0.zip
    unzip gatk-4.0.6.0.zip
    cd gatk-4.0.6.0
    mv gatk-package-4.0.6.0-local.jar /usr/local/lib
    echo 'export GATK_LOCAL_JAR=/usr/local/lib/gatk-package-4.0.6.0-local.jar' >> $SINGULARITY_ENVIRONMENT
    mv gatk /usr/local/bin
    mv gatk-completion.sh /usr/local/bin
    cd /usr/local/src

    #SnpEff
    wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip
    cd snpEff
    mv snpEff.jar /usr/local/lib
    chmod a+x /usr/local/lib/snpEff.jar
    echo -e "#! /bin/bash\njavamem=""\nif [[ \$1 =~ "-Xmx" ]];then javamem=\$1; shift 1; fi\nexec java \$javamem -jar /usr/local/lib/snpEff.jar "\$@"" > /usr/local/bin/snpEff
    chmod a+x /usr/local/bin/snpEff

    #plink-ng
    wget https://github.com/chrchang/plink-ng/archive/b0cec5e.tar.gz
    tar xf b0cec5e.tar.gz
    cd plink-ng-b0cec5e/2.0/build_dynamic
    make
    mv plink2 /usr/local/bin
    mv pgen_compress /usr/local/bin
    cd /usr/local/src

    #gtfToGenePred
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
    chmod a+x gtfToGenePred
    mv gtfToGenePred /usr/local/bin/

    #pigz
    wget https://github.com/madler/pigz/archive/v2.4.tar.gz
    tar xf v2.4.tar.gz
    cd pigz-2.4
    make
    mv pigz /usr/local/bin
    mv unpigz /usr/local/bin
    cd /usr/local/src

    #NextFlow
    curl -s https://get.nextflow.io | bash
    mv nextflow /usr/local/bin/
    chmod 777 /usr/local/bin/nextflow

    #SRA toolkit
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
    tar xf sratoolkit.current-centos_linux64.tar.gz
    cd $(ls -d sratoolkit* | grep -v .gz)
    mv ./bin/fasterq-dum* /usr/local/bin/
    mv ./bin/prefetc* /usr/local/bin/
    cd /usr/local/src

%runscript
    #if need to run stuff, put here
