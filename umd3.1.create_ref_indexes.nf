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
  log.info 'nextflow run umd3.1.create_ref_indexes.nf \
            --dataDir data \
            --fa Bos_taurus.UMD3.1.dna.toplevel.fa.gz \
            --gtf Bos_taurus.UMD3.1.92.gtf.gz \
            --vcf bos_taurus_incl_consequences.vcf.gz \
            --bed 130604_Btau_UMD3_Exome_BM_EZ_HX1.bed.gz \
            -c "bovine_DNA_RNA.nextflow.config" \
            -with-report "ref.report.html" \
            -with-timeline "ref.timeline.html"'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '--dataDir    STRING    where data is found, and written'
  exit 1
}

/* 1.0: Channels from files to unpigz
* sorts fasta on chrnames, then reorders
*/
FA = Channel.fromPath("$params.dataDir/$params.fa", type: "file")
GTF = Channel.fromPath("$params.dataDir/$params.gtf", type: "file")
BED = Channel.fromPath("$params.dataDir/$params.bed", type: "file")

process sortfa {

  publishDir "$params.dataDir/", mode: "copy", pattern: "*"

  input:
  file(fagz) from FA
  file(gtfgz) from GTF
  file(bedgz) from BED

  output:
  set file('*.sort.fa'), file('*.sort.fa.fai') into (bwa_fasta, star_fasta)
  file('*.gtf') into (star_gtf, refflat_gtf, rrna_gtf)
  file('*.bed') into exome_bed
  file('*.dict') into (fasta_dict, rrna_dict)

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

  publishDir "$params.dataDir", mode: "copy", pattern: "*"

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

  publishDir "$params.dataDir", mode: "copy", pattern: "*"

  input:
  set file(fasta), file(faidx) from star_fasta
  file(gtf) from star_gtf

  output:
  file('*') into complete1_2

  script:
  """
  mkdir STAR_99
  STAR --runMode genomeGenerate \
       --genomeDir STAR_99 \
       --genomeFastaFiles $fasta \
       --sjdbGTFfile $gtf \
       --sjdbOverhang 99
  """
}
complete1_2.subscribe { println "Completed STAR indexing" }

/* 1.3: BED intervalList
*
*/
process intlist {

  publishDir "$params.dataDir", mode: "copy", pattern: "*"

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

  publishDir "$params.dataDir", mode: "copy", pattern: "*"

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

  publishDir "$params.dataDir", mode: "copy", pattern: "*"

  input:
  file(dict) from rrna_dict
  file(gtf) from rrna_gtf

  output:
  file('*') into complete1_5

  script:
  """
  RRNA=\$(echo $gtf | sed 's/gtf/rRNA.interval_list/')
  cat $dict > \$RRNA
  grep "gene_biotype \"rRNA\"" $gtf | \
  perl -ane 'if(\$F[2] eq "transcript"){
    \$gn=\$F[13]; \$gn=~s/\"//g; print "\$F[0]\t\$F[3]\t\$F[4]\t\$F[6]\t\$gn\n";
  }' >> \$RRNA
  """
}
complete1_6.subscribe { println "Completed rRNA interval_list" }
