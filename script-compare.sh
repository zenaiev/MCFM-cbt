#!/bin/bash
cleandir=/home/zenaiev/sven-msbarm/210219/MCFM-6.8-addon-clean
for file in `find src`; do 
  if [ ! -f $file ]; then continue; fi; 
  newf="$cleandir/$file"; 
  #echo "*** $file ***"; diff $file $newf; 
  if [ `diff $file $newf | wc -l` -ge 1 ]; then
    echo "cp $newf $file"
  fi
done
