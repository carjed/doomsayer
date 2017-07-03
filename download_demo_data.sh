#!/bin/bash

curdir=${PWD}
refdir="$curdir/demo/input"
mkdir $refdir

echo "Entering $refdir"
cd $refdir

wget -P $refdir -r --no-parent -nH --cut-dirs=2 --reject="index.html*" \
  http://mutation.sph.umich.edu/share/doomsayer_demo/

cd $curdir
echo "Returning to $curdir"
