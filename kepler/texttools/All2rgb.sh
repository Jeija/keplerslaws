#!/bin/sh
for file in Untertitel*tex
do
   if test -f $file; then
      Untertitel2rgb.sh $file
   fi
done
for file in Text*tex
do
   if test -f $file; then
      Text2rgb.sh $file
   fi
done
for file in Crop*tex
do
   if test -f $file; then
      Crop2rgb.sh $file
   fi
done
