#!/bin/bash
# This is peaks-generating routine for MACS2 hg peak calling. 
# It uses SGE scheduler for multiple parallel calling. This is a part
# of analysis pipeline.
# Requires <contrasts.txt> file to be present in working directory 
# which is a simple tab/space-delimited text file with two columns where 
# the first column is path to TEST (Chip) sam/bam file and the second 
# is path to corresponding input file. 
# MACS Q or P value can be supplied as argument.
# It assumes <macs2> path is in PATH env variable.
# Naming output file: Test-X-Control
#
#
UTIL=/home/ppl/utils
[ -e contrasts.txt ] || exit 1
REF=hs
ARGQ=
ARGP=
while getopts ':q:p:' arg
    do
      case $arg in
      q)
        ARGQ=$OPTARG
        ;;
      p)
        ARGP=$OPTARG
      ;;
      :)
        echo 'Invalid argument'
        exit 2
      ;;
  esac
done
shift $((OPTIND-1))

while read -r X CX
  do 
    S=`basename $X .bam`
    C=`basename $CX .bam`
    if [[ ! "$ARGQ" && ! "$ARGP" ]]
      then
        Q=0.0001
    else
      Q=$ARGQ
      P=$ARGP
    fi
    echo making $S-$C ...
    [[ -e $X && -e $CX ]] || { echo "Either TEST or INPUT is missing"; break; }
    if [ "$Q" ]
      then
     echo "macs2 callpeak -t $X -c $CX --tempdir /storage/data01/tmp -f BAM -g $REF -n ${S}-X-${C} -q $Q -m 5 50 --nomodel" | qsub -q all.q -cwd -V -e /home/ppl/log/mkpeaks.e -o /home/ppl/log/mkpeaks.o
    elif [ "$P" ]
      then
    echo "macs2 callpeak -t $X -c $CX --tempdir /storage/data01/tmp -f BAM -g $REF -n ${S}-X-${C} -p $P -m 5 50 --nomodel" | qsub -q all.q -cwd -V -e /home/ppl/log/mkpeaks.e -o /home/ppl/log/mkpeaks.o
    fi
done < contrasts.txt
##################################END#################################
