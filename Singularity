Bootstrap:docker
From:centos:centos7.4.1708

%help
    Container for Bovine DNA and RNA tools

%labels
    MAINTAINER Colin Sauze

%environment
    #define environment variables here

%post  
    #essential utilities
     yum -y install git wget bzip2 unzip
# 
#     #language and libraries
     yum -y install java-1.8.0-openjdk-devel gcc gcc-c++ glibc-devel make ncurses-devel zlib-devel libbzip2-devel bzip2-devel xz-devel
#     #yum -y install https://centos7.iuscommunity.org/ius-release.rpm
#     yum -y install python36u python36u-pip python36u-devel
# 

     mkdir -p /usr/local/src
     cd /usr/local/src


#     #samtools
#     wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
#     tar xvf samtools-1.8.tar.bz2
#     cd samtools-1.8
#     ./configure --prefix=/usr/local/
#     make
#     make install
#     cd ..
# 
#     #bcftools
#     wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2
#     tar xvf bcftools-1.8.tar.bz2
#     cd bcftools-1.8
#     ./configure --prefix=/usr/local
#     make
#     make install
#     cd ..
# 
#     #htslib
#     wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2
#     tar xvf htslib-1.8.tar.bz2
#     cd htslib-1.8
#     ./configure --prefix=/usr/local
#     make
#     make install
#     cd ..
# 
#     #bedtools
#     wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz
#     tar xvf bedtools-2.27.1.tar.gz
#     cd bedtools2
#     make
#     make install
#     cd ..
    
#    #picard
#    wget https://github.com/broadinstitute/picard/releases/download/2.18.9/picard.jar -O /usr/local/lib/picard.jar
   
#     #BWA
#     wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
#     tar xvf bwa-0.7.17.tar.bz2
#     cd bwa-0.7.17
#     make
#     mv bwa /usr/local/bin
#     mv *.pl /usr/local/bin
#     mv libbwa.a /usr/local/lib
#     cd ..
     
#     #STAR
#     wget https://github.com/alexdobin/STAR/archive/2.6.0c.tar.gz
#     tar xvf 2.6.0c.tar.gz
#     cd STAR-2.6.0c/source
#     make STAR
#     mv STAR /usr/local/bin
#     #make install
#     cd ..
#     
#     #pigz
#     wget https://github.com/madler/pigz/archive/v2.4.tar.gz
#     tar xvf v2.4.tar.gz
#     cd pigz-2.4
#     make
#     ls -lrt
#     mv pigz /usr/local/bin
#     mv unpigz /usr/local/bin
#     cd ..
#     
    #kallisto
#    wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz
#   tar xvf kallisto_linux-v0.44.0.tar.gz
#   mv kallisto_linux-v0.44.0/kallisto /usr/local/bin
    
#     #fastqc
#     wget https://github.com/s-andrews/FastQC/archive/v0.11.7.tar.gz
#     tar xvf v0.11.7.tar.gz
#     cd FastQC-0.11.7
#     mv *.jar /usr/local/lib
#     mv fastqc /usr/local/bin
#     chmod 755 /usr/local/bin/fastqc
#     cd ..
    
#    #multiqc
#    pip3.6 install multiqc
    
#     #BBMap
#     wget "https://downloads.sourceforge.net/project/bbmap/BBMap_38.11.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbbmap%2Ffiles%2FBBMap_38.11.tar.gz%2Fdownload&ts=1531223392" -O BBMap_38.11.tar.gz
#    
#     tar xvf BBMap_38.11.tar.gz
#     #its all shell scripts, don't really know where to put them
#     
#     #gatk
#     wget https://github.com/broadinstitute/gatk/releases/download/4.0.6.0/gatk-4.0.6.0.zip
#     unzip gatk-4.0.6.0.zip
#     cd gatk-4.0.6.0
#     mv gatk-package-4.0.6.0-local.jar /usr/local/lib
#     mv gatk /usr/local/bin
#     mv gatk-completion.sh /usr/local/bin
#     cd ..
#     #don't think we need anything from the scripts subdir
#     
    #SnpEff
    wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip
    #doesn't really make sense to put this elsewhere in /usr/local, so just leaving it in /usr/local/src
    
    

    
    
    
%runscript
    #set locale so multiqc doesn't complain
    export LANG=en_US.UTF-8
    
    #running stuff

