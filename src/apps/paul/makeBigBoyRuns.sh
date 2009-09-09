#!/bin/bash

execFile="bigboylb.cc"
cp $execFile $execFile.bak

for lbmode in 2 1 0; do
  sed "s|LBMODE=[0-9]|LBMODE=$lbmode|g" $execFile > $execFile.new
  mv $execFile.new $execFile
  make bigboylb
  cp bigboylb bigboylb_lb$lbmode
done

