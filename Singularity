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
    yum -y install git wget bzip2 unzip which

    #language and libraries
    yum -y install java-1.8.0-openjdk-devel gcc gcc-c++ glibc-devel make ncurses-devel zlib-devel libbzip2-devel bzip2-devel xz-devel perl-DBI lapack-devel atlas-devel
    #libclas and libatlas aren't put in the right places
    ln -s /usr/lib64/atlas/libtatlas.so /usr/lib64/libatlas.so
    ln -s /usr/lib64/atlas/libsatlas.so /usr/lib64/libcblas.so

    yum -y install https://centos7.iuscommunity.org/ius-release.rpm
    yum -y install python36u python36u-pip python36u-devel

    mkdir -p /usr/local/src
    cd /usr/local/src
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib


    #NextFlow
    curl -s https://get.nextflow.io | bash
    mv nextflow /usr/local/bin/

    #for required data files
    mkdir /data
    cd /data

    #NextFlow reference indexing/manipulation script, config
    wget https://raw.githubusercontent.com/brucemoran/Bovine_DNA_RNA/master/umd3.1.create_ref_indexes.simg.nf
    wget https://raw.githubusercontent.com/brucemoran/Bovine_DNA_RNA/master/bovine_DNA_RNA.nextflow.simg.config

    #only needed to speed things up in aber
    #export http_proxy="http://wwwcache.aber.ac.uk:8080"

    #1000 Bulls
    wget http://www.1000bullgenomes.com/doco/1000bulls_v6_annotated_snps.tab.gz
    wget http://www.1000bullgenomes.com/doco/1000bulls_v6_annotated_indels.tab.gz

    #DNA genome fasta (toplevel no masking)
    wget ftp://ftp.ensembl.org/pub/release-92/fasta/bos_taurus/dna/Bos_taurus.UMD3.1.dna.toplevel.fa.gz

    #cDNA genome fasta
    wget ftp://ftp.ensembl.org/pub/release-92/fasta/bos_taurus/cdna/Bos_taurus.UMD3.1.cdna.all.fa.gz

    #CDS fasta
    wget ftp://ftp.ensembl.org/pub/release-92/fasta/bos_taurus/cds/Bos_taurus.UMD3.1.cds.all.fa.gz

    #ncRNA fasta
    wget ftp://ftp.ensembl.org/pub/release-92/fasta/bos_taurus/ncrna/Bos_taurus.UMD3.1.ncrna.fa.gz

    #GTF
    wget ftp://ftp.ensembl.org/pub/release-92/gtf/bos_taurus/Bos_taurus.UMD3.1.92.abinitio.gtf.gz
    wget ftp://ftp.ensembl.org/pub/release-92/gtf/bos_taurus/Bos_taurus.UMD3.1.92.chr.gtf.gz
    wget ftp://ftp.ensembl.org/pub/release-92/gtf/bos_taurus/Bos_taurus.UMD3.1.92.gtf.gz

    #variants
    wget ftp://ftp.ensembl.org/pub/release-92/variation/vcf/bos_taurus/bos_taurus_incl_consequences.vcf.gz

    #exome
    wget https://raw.githubusercontent.com/brucemoran/Bovine_DNA_RNA/master/130604_Btau_UMD3_Exome_BM_EZ_HX1.bed.gz

    #run NextFlow reference script
    #first set ulimit for allowing coredump of LevelDB kills
    ulimit -c unlimited
    /usr/local/bin/nextflow run umd3.1.create_ref_indexes.simg.nf \
      --dataDir /data \
      --fa Bos_taurus.UMD3.1.dna.toplevel.fa.gz \
      --gtf Bos_taurus.UMD3.1.92.gtf.gz \
      --vcf bos_taurus_incl_consequences.vcf.gz \
      --bed 130604_Btau_UMD3_Exome_BM_EZ_HX1.bed.gz \
      -c "bovine_DNA_RNA.nextflow.simg.config" \
      -with-report "ref.report.html" \
      -with-timeline "ref.timeline.html"

%runscript
    #set locale so multiqc doesn't complain
    export LANG=en_US.UTF-8

    #running stuff
