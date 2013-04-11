#!/bin/sh
base=`basename $1 .tex`
latex $base
dvips -ta3 -o $base.ps $base
ppmmake rgbi:0/0/0 800 600 > $base.ppm_
ps2ppm.sh < $base.ps | pnmcrop > $base.ppm__
pnmpaste $base.ppm__ 0 0 $base.ppm_ | pnmcrop > $base.ppm
convert $base.ppm sgi:$base.rgb
rm $base.ppm
rm $base.ps
rm $base.ppm__
rm $base.ppm_ 
rm $base.aux
rm $base.log
rm $base.dvi
