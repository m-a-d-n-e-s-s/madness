#!/bin/bash

execFile="testshqi.cc"
cp $execFile $execFile.bak

for lbmode in 2 1; do
  sed "s|LBMODE=[0-9]|LBMODE=$lbmode|g" $execFile > $execFile.new
  mv $execFile.new $execFile
  make testshqi
  cp testshqi testshqi_lb$lbmode
done

