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
    yum -y install git wget bzip2 

    #language and libraries
    yum -y install java-1.8.0-openjdk-devel gcc gcc-c++ glibc-devel make ncurses-devel zlib-devel libbzip2-devel bzip2-devel xz-devel
    #yum -y install https://centos7.iuscommunity.org/ius-release.rpm
    yum -y install python36u python36u-pip python36u-devel

    #samtools
    mkdir -p /usr/local/src
    cd /usr/local/src
    wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
    tar xvf samtools-1.8.tar.bz2
    cd samtools-1.8
    ./configure --prefix=/usr/local/
    make
    make install
    cd ..

    #bcftools
    wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2
    tar xvf bcftools-1.8.tar.bz2
    cd bcftools-1.8
    ./configure --prefix=/usr/local
    make
    make install
    cd ..

    #htslib
    wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2
    tar xvf htslib-1.8.tar.bz2
    cd htslib-1.8
    ./configure --prefix=/usr/local
    make
    make install
    cd ..
#

%runscript
    #running stuff

