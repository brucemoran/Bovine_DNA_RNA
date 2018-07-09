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
#     yum -y install git wget bzip2 
# 
#     #language and libraries
#     yum -y install java-1.8.0-openjdk-devel gcc gcc-c++ glibc-devel make ncurses-devel zlib-devel libbzip2-devel bzip2-devel xz-devel
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
     
    #STAR
    wget https://github.com/alexdobin/STAR/archive/2.6.0c.tar.gz
    tar xvf 2.6.0c.tar.gz
    cd STAR-2.6.0c/source
    make STAR
    mv STAR /usr/local/bin
    #make install
    cd ..


%runscript
    #running stuff

