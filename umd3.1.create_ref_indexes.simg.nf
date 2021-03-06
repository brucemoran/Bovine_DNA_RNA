#!/usr/bin/env nextflow

/* indexes required for exome, RNAseq
* BWA, STAR, exome interval list
*/

params.help = ""

if (params.help) {
  log.info ''
  log.info '---------------------------------'
  log.info 'NEXTFLOW EXOME, RNASEQ REFERENCES'
  log.info '---------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run umd3.1.create_ref_indexes.simg.nf \
            --refDir refs \
            -c "umd3.1.RNAseq_STAR.exome_BWA.GATK4-HC.simg.nextflow.config" \
            -with-report "ref.report.html" \
            -with-timeline "ref.timeline.html"'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '--refDir    STRING    where reference data will be store'
  exit 1
}

/* 0.0: Download ref data
*/

process downloadFa {

  output:
  file('*.fa.gz') into fa

  script:
  """
  wget -O ./Bos_taurus.UMD3.1.dna.toplevel.fa.gz \
    ftp://ftp.ensembl.org/pub/release-92/fasta/bos_taurus/dna/Bos_taurus.UMD3.1.dna.toplevel.fa.gz
  """
}

process downloadGtf {


  output:
  file('*.gtf.gz') into gtf

  script:
  """
  wget -O ./Bos_taurus.UMD3.1.92.gtf.gz \
    ftp://ftp.ensembl.org/pub/release-92/gtf/bos_taurus/Bos_taurus.UMD3.1.92.gtf.gz
  """
}

process downloadBed {

  output:
  file('*.bed.gz') into bed

  script:
  """
  wget -O ./130604_Btau_UMD3_Exome_BM_EZ_HX1.bed.gz \
    https://raw.githubusercontent.com/brucemoran/Bovine_DNA_RNA/master/130604_Btau_UMD3_Exome_BM_EZ_HX1.bed.gz
  """
}

process downloadVcf {

  publishDir "$params.refDir", mode: "copy", pattern: "*corr.sd.v*"

  output:
  file('*') into complete0_1
  file(dict) from vcf_dict

  script:
  """
  wget -O ./bos_taurus_incl_consequences.vcf.gz \
    ftp://ftp.ensembl.org/pub/release-92/variation/vcf/bos_taurus/bos_taurus_incl_consequences.vcf.gz
  ##need to remove empty ALT
  gunzip -c bos_taurus_incl_consequences.vcf.gz | \
  gawk 'BEGIN{FS="\t"; OFS="\t"}{if (NF>1 && \$5=="") {\$5="."; print \$0} else print \$0}' > bos_taurus_incl_consequences.corr.vcf

  picard-tools UpdateVcfSequenceDictionary \
    I=bos_taurus_incl_consequences.corr.vcf \
    O=bos_taurus_incl_consequences.corr.sd.vcf \
    SD=$dict

  bgzip bos_taurus_incl_consequences.corr.sd.vcf
  tabix bos_taurus_incl_consequences.corr.sd.vcf.gz
  """
}
complete0_1.subscribe { println "Completed VCF download, tabix indexing" }

/* 1.0: Channels from files to unpigz
* sorts fasta on chrnames, then reorders
*/
process sortfa {

  publishDir "$params.refDir", mode: "copy", pattern: "*[.dict,*.sort.fa,*.sort.fa.fai,*.bed,*.gtf]"

  input:
  file(fagz) from fa
  file(gtfgz) from gtf
  file(bedgz) from bed

  output:
  set file('*.sort.fa'), file('*.sort.fa.fai') into (bwa_fasta, star_fasta)
  file('*.gtf') into (star_gtf, refflat_gtf, rrna_gtf)
  file('*.bed') into exome_bed
  file('*.dict') into (fasta_dict, rrna_dict, vcf_dict)

  script:
  """
  BED=\$(echo $bedgz | sed 's/.gz\$//')
  unpigz -c $bedgz | sort -Vk1,1 -k2,2n | uniq | perl -ane 'if(scalar(@F) == 3){print \$_;}' > \$BED

  unpigz -kf $fagz
  unpigz -kf $gtfgz

  FA=\$(echo $fagz | sed 's/.gz\$//')
  samtools faidx \$FA
  SORTFA=\$(echo \$FA | sed 's/.fa/.sort.fa/')
  grep ">" \$FA | sort -V > sort.chrs

  perl -ane 'chomp; @s=split(/:/);print "\$s[3]:\$s[4]-\$s[5]\\t\$_\\n";' sort.chrs | \
  while read LINE; do
    REG=\$(echo "\$LINE" | cut -f 1)
    NAM=\$(echo "\$LINE" | cut -f 2)
    samtools faidx \$FA \$REG > tmp.fa
    sed -i "1s/.*/\$NAM/" tmp.fa
    cat tmp.fa >> \$SORTFA
  done

  samtools faidx \$SORTFA
  picard-tools CreateSequenceDictionary R=\$SORTFA
  """
}

/* 1.1: Index BWA
*
*/
process bwaidx {

  publishDir "$params.refDir", mode: "copy", pattern: "*"

  input:
  set file(fasta), file(faidx) from bwa_fasta

  output:
  file('*') into complete1_1

  script:
  """
  bwa index $fasta
  """
}
complete1_1.subscribe { println "Completed BWA indexing" }

/* 1.2: Index STAR
*
*/
process staridx {

  publishDir "$params.refDir", mode: "copy", pattern: "*"

  input:
  set file(fasta), file(faidx) from star_fasta
  file(gtf) from star_gtf

  output:
  file('*') into complete1_2

  script:
  """
  mkdir STAR_${params.STARsjdbOverhang}
  STAR --runMode genomeGenerate \
       --genomeDir STAR_${params.STARsjdbOverhang} \
       --genomeFastaFiles $fasta \
       --sjdbGTFfile $gtf \
       --sjdbOverhang ${params.STARsjdbOverhang}
  """
}
complete1_2.subscribe { println "Completed STAR indexing" }

/* 1.3: BED intervalList
*
*/
process intlist {

  publishDir "$params.refDir", mode: "copy", pattern: "*.interval_list"

  input:
  file(fastadict) from fasta_dict
  file(exomebed) from exome_bed

  output:
  file('*') into complete1_3

  script:
  """
  INTLIST=\$(echo $fastadict | sed 's/dict/interval_list/')
  cat $fastadict > \$INTLIST
  perl -ane 'chomp; print \$_ . "\\t+\\tinterval_\$.\\n";' $exomebed >> \$INTLIST
  """
}
complete1_3.subscribe { println "Completed exome interval_list" }

/* 1.4: refFlat for RNAseq
*
*/
process refFlat {

  publishDir "$params.refDir", mode: "copy", pattern: "*.refFlat"

  input:
  file(gtf) from refflat_gtf

  output:
  file('*') into complete1_4

  script:
  """
  REFFLAT=\$(echo $gtf | sed 's/gtf/refFlat/')
   gtfToGenePred \
    -genePredExt \
    -geneNameAsName2 \
    -ignoreGroupsWithoutExons \
    $gtf \
    /dev/stdout | \
    awk 'BEGIN { OFS="\t"} {print \$12, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10}' \
    > \$REFFLAT
  """
}
complete1_4.subscribe { println "Completed refFlat" }

/* 1.5: rRNA for RNAseq
*
*/
process rRNA {

  publishDir "$params.refDir", mode: "copy", pattern: "*.rRNA.interval_list"

  input:
  file(dict) from rrna_dict
  file(gtf) from rrna_gtf

  output:
  file('*') into complete1_5

  script:
  """
  RRNA=\$(echo $gtf | sed 's/gtf/rRNA.interval_list/')
  cat $dict > \$RRNA
  grep "gene_biotype \\"rRNA\\"" $gtf | \
  perl -ane 'if(\$F[2] eq "transcript"){
    \$gn=\$F[13]; \$gn=~s/\"//g; \$gn=~s/;//g; print "\$F[0]\\t\$F[3]\\t\$F[4]\\t\$F[6]\\t\$gn\\n";
  }' | sort -V >> \$RRNA
  """
}
complete1_5.subscribe { println "Completed rRNA interval_list" }
