#!/bin/bash

curdir=${PWD}
refdir="$curdir/demo/input"
mkdir $refdir

echo $refdir
cd $refdir

# curl -s "http://mutation.sph.umich.edu/hg19/GRCh37_RefSeq_sorted.bed" >  "$refdir/GRCh37_RefSeq_sorted.bed"

wget -P $refdir -r --no-parent -nH --cut-dirs=2 --reject="index.html*" \
  http://mutation.sph.umich.edu/share/doomsayer_demo/

cd $curdir
echo $curdir
