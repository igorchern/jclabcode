#!/bin/bash
## bowtie2 paired-end hg38 alignment script to run with SLURM batch job
## ARGUMENTS: 
## 1) = fq file 1 (extension .fq ot .fq.gz)
## 2) = fq file 2 (extension .fq ot .fq.gz)
## 3) = OUTPUT directory where pipeline files will be stored
## 4) trim (optional) in format: trim=[l|r]<read length after trim> where <l> means trim from 5' and <r> from 3' end correspondingly
## 5) genome coverage file (bigwig) (optional) in format: bigwig=yes

## SBATCH headers omitted here as they are provided by pipeline configuration database.  
##########################################################
BOWTIE=/home/hpcuser/local/bin/bowtie2
GREF=/mnt/scratchc/jclab/Ref/Homo_sapiens_GRC_h38/Bowtie2
GSIZE=/mnt/scratchc/jclab/Ref/Homo_sapiens_GRC_h38/hg38.chrom.size 
SAMTOOLS=/home/hpcuser/local/bin/samtools
BEDTOOLS=/home/hpcuser/local/bin/bedtools
UBIN=/home/hpcuser/local/bin
TRIM=
TMP=/tmp/hpcuser
UTILS=/home/hpcuser/utils
GENCOV=no
JAVA=/usr/bin/java
. /home/hpcuser/pipeline/slurm/bin/utils

## at least 3 args should be provided (specified above)
[[ "$1" && "$2" && "$3" ]] || exit -1

OUTD=$3

## Using TMP local patrition:
mkdir $TMP 2>/dev/null
cd $TMP
## copy files over local network to tmp working directory
if ! rsync $1 . 
  then
  echo FAILED > $OUTD/$BN.status
  exit 1
fi
if ! rsync $2 . 
  then
  echo FAILED > $OUTD/$BN.status
  exit 2
fi
F1=$1
F1=${F1##*/}
F2=$2
F2=${F2##*/}

UN=$(printf '%s\n%s\n' $F1 $F2 | sed -E 's,_[12]([.:]),\1,;s,\.fq.*$,,' | sort -u)
BN=`basename $UN`
BAM=pipeline_$$.$BN.bam
BGR=pipeline_$$.bgr
BIGWIG=

if [ $# -gt 3 ]               # test if there are some other arguments 
  then 
    for VAR in $@
      do 
       if [[ $VAR == trim=* ]] 
        then
          S=`echo $VAR | sed 's,trim=,,'`
          P=${S:0:1} 
          R=${S:1}
          echo Trimming reads to $R bp : fqtrim.${P}.sh ... 

          $UBIN/fqtrim.${P}.sh $F1 $R > ${F1%.gz}.tmp$$
          rm -f $F1
          mv ${F1%.gz}.tmp$$ ${F1%.gz}
          F1=${F1%.gz}

          $UBIN/fqtrim.${P}.sh $F2 $R > ${F2%.gz}.tmp$$
          rm -f $F2
          mv ${F2%.gz}.tmp$$ ${F2%.gz}
          F2=${F2%.gz}

        elif [[ $VAR == bigwig=* ]]
          then
            BIGWIG=`echo $VAR | sed 's,bigwig=,,'`
        fi
     done
fi

echo RUNNING > $OUTD/$BN.status
echo cmd: "$BOWTIE --rg-id REF --rg RN:hg38 -p 4 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x $GREF/hg38 -1 $F1 -2 $F2 $SAMTOOLS view -1 -S -o $BAM -" 
$BOWTIE --rg-id REF --rg RN:hg38 -p 4 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x $GREF/hg38 -1 $F1 -2 $F2 | $SAMTOOLS view -1 -S -o $BAM -
if [ ! -s $BAM ]
then
    rm -f *.bam
    echo FAILED > $OUTD/$BN.status
    exit 3
fi

$SAMTOOLS sort -o $BN.srt.bam $BAM 

#########################################
## Generating genome alignment coverage
#########################################
if [[ $BIGWIG == yes ]]
  then
  $BEDTOOLS genomecov -split -bg -g $GSIZE -ibam $BN.srt.bam | sort -T . -k1,1 -k2,2n > $BGR
  $UBIN/bedGraphToBigWig $BGR $GSIZE $OUTD/$BN.bw 
## cleaning up
  rm $BGR
fi
#################
## copy bam file to output directory
cp $BN.srt.bam $OUTD/

## Cleaning up
rm -f $BAM $F1 $F2 $BN.srt.bam
echo DONE > $OUTD/$BN.status
##################################END##################################
