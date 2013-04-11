/* Library of useful astronomical routines for C/C++ programming

   Copyright (C) 1999 Tobias Kramer

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public
   License along with this library; see the file COPYING.LIB.  If not,
   write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  
*/

/*****************************************************************************/
/* Module: PPOS_LO1.CPP                                                      */
/* Version 1.0                                                               */
/* Last modified: February 22, 1993                                          */
/*****************************************************************************/
/* Abkuerzungen (Literatur):                                                 */
/*                                                                           */
/* AmdPC : Montenbruck,Pfleger: Astronomie mit dem Personal Computer         */
/*****************************************************************************/

#include <math.h>
#include "ppos_lo.h"

/*****************************************************************************/
/* Sonne-, Mond- und Planetenreihenentwicklung nach                          */
/*                                                                           */
/* Low-Precision Formulae for Planetary Positions                            */
/* T.C.van Flandern, K.F.Pulkkinen                                           */
/* The Astrophysical Journal Supplement Series, 41:391-411, 1979 November    */
/*                                                                           */
/* Genauigkeit (1680-2280):                                                  */
/* Laenge und Breite besser als 1' (Pluto 15')                               */
/*                                                                           */
/* Eingabe:                                                                  */
/* jd : Julianisches (Ephemeriden-) Datum                                    */
/*                                                                           */
/* Return:                                                                   */
/* Geometrische ekliptikale Koordinaten fuer das Aequinoktium des Datums     */
/*****************************************************************************/

sphaer son( double jd )

{ double t,d,lm,ls,gs,g2,g4,g5;
  sphaer s;

  d= jd-2451545.0;
  t= d/36525.0+1.0;
  lm= M_2PI*frac(0.606434+0.03660110129*d); /*  1 */
  ls= M_2PI*frac(0.779072+0.00273790931*d); /*  7 */
  gs= M_2PI*frac(0.993126+0.00273777850*d); /*  8 */
  g2= M_2PI*frac(0.140023+0.00445036173*d); /* 13 */
  g4= M_2PI*frac(0.053856+0.00145561327*d); /* 16 */
  g5= M_2PI*frac(0.056531+0.00023080893*d); /* 19 */
/*                      1=lm   7=ls   8=gs  13=g2  16=g4  19=g5       */
s.lambda= 6910.0*    sin(                  gs                     )+
	 72.0*    sin(              2.0*gs                     )+
	-17.0*t*  sin(                  gs                     )+
	 -7.0*    cos(                  gs                  -g5)+
	  6.0*    sin(    lm    -ls                            )+
	  5.0*    sin(              4.0*gs       -8.0*g4+3.0*g5)+
	 -5.0*    cos(              2.0*gs-2.0*g2              )+
	 -4.0*    sin(                  gs    -g2              )+
	  4.0*    cos(              4.0*gs       -8.0*g4+3.0*g5)+
	  3.0*    sin(              2.0*gs-2.0*g2              )+
	 -3.0*    sin(                                       g5)+
	 -3.0*    sin(              2.0*gs              -2.0*g5);
s.lambda=mod(ls+s.lambda/648000.0*PI,M_2PI);
s.beta=0.0;
s.r=1.00014                                                   +
     -0.01675*    cos(                  gs                     )+
     -0.00014*    cos(              2.0*gs                     );
  return(s);
}

sphaer mer( double jd )

{ double t,d,l1,g1,g2,f1;
  sphaer s;
  d= jd-2451545.0;
  t= d/36525.0+1.0;
  l1= M_2PI*frac(0.700695+0.01136771400*d); /*  9 */
  g1= M_2PI*frac(0.485541+0.01136759566*d); /* 10 */
  f1= M_2PI*frac(0.566441+0.01136762384*d); /* 11 */
  g2= M_2PI*frac(0.140023+0.00445036173*d); /* 13 */
/*                     10=g1  11=f1  13=g2         */
s.lambda=84378.0*    sin(    g1              )+
      10733.0*    sin(2.0*g1              )+
       1892.0*    sin(3.0*g1              )+
       -646.0*    sin(       2.0*f1       )+
        381.0*    sin(4.0*g1              )+
       -306.0*    sin(    g1-2.0*f1       )+
       -274.0*    sin(    g1+2.0*f1       )+
        -92.0*    sin(2.0*g1+2.0*f1       )+
         83.0*    sin(5.0*g1              )+
        -28.0*    sin(3.0*g1+2.0*f1       )+
         25.0*    sin(2.0*g1-2.0*f1       )+
         19.0*    sin(6.0*g1              )+
         -9.0*    sin(4.0*g1+2.0*f1       )+
          8.0*t*  sin(    g1              )+
          7.0*    cos(2.0*g1       -5.0*g2);
s.lambda=mod(l1+s.lambda/648000.0*PI,M_2PI);
s.beta=24134.0*    sin(           f1       )+
       5180.0*    sin(    g1    -f1       )+
       4910.0*    sin(    g1    +f1       )+
       1124.0*    sin(2.0*g1    +f1       )+
        271.0*    sin(3.0*g1    +f1       )+
        132.0*    sin(2.0*g1    -f1       )+
         67.0*    sin(4.0*g1    +f1       )+
         18.0*    sin(3.0*g1    -f1       )+
         17.0*    sin(5.0*g1    +f1       )+
        -10.0*    sin(       3.0*f1       )+
         -9.0*    sin(    g1-3.0*f1       );
s.beta=s.beta/648000.0*PI;
s.r=0.39528                              +
     -0.07834*    cos(    g1              )+
     -0.00795*    cos(2.0*g1              )+
     -0.00121*    cos(3.0*g1              )+
     -0.00022*    cos(4.0*g1              );
  return(s);
}

sphaer ven( double jd )

{ double t,d,gs,l2,g2,f2;
  sphaer s;

  d= jd-2451545.0;
  t= d/36525.0+1.0;
  gs= M_2PI*frac(0.993126+0.00273777850*d); /*  8 */
  l2= M_2PI*frac(0.505498+0.00445046867*d); /* 12 */
  g2= M_2PI*frac(0.140023+0.00445036173*d); /* 13 */
  f2= M_2PI*frac(0.292498+0.00445040017*d); /* 14 */
/*                      8=gs  13=g2  14=f2         */
s.lambda= 2814.0*    sin(           g2       )+
       -181.0*    sin(              2.0*f2)+
        -20.0*t*  sin(           g2       )+
         12.0*    sin(       2.0*g2       )+
        -10.0*    cos(2.0*gs-2.0*g2       )+
          7.0*    cos(3.0*gs-3.0*g2       );
s.lambda=mod(l2+s.lambda/648000.0*PI,M_2PI);
s.beta=12215.0*    sin(                  f2)+
         83.0*    sin(           g2    +f2)+
         83.0*    sin(           g2    -f2);
s.beta=s.beta/648000.0*PI;
s.r=0.72335                              +
     -0.00493*    cos(           g2       );
  return(s);
}
