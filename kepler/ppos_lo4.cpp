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
/* Module: PPOS_LO4.CPP                                                      */
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

sphaer ura( double jd )

{ double t,d,g5,g6,g7,f7,g8,l7;
  sphaer s;

  d= jd-2451545.0;
  t= d/36525.0+1.0;
  g5= M_2PI*frac(0.056531+0.00023080893*d); /* 19 */
  g6= M_2PI*frac(0.882987+0.00009294371*d); /* 22 */
  l7= M_2PI*frac(0.870169+0.00003269438*d); /* 24 */
  g7= M_2PI*frac(0.400589+0.00003269438*d); /* 25 */
  f7= M_2PI*frac(0.664614+0.00003265562*d); /* 26 */
  g8= M_2PI*frac(0.725368+0.00001672092*d); /* 28 */

/*                     19=g5  22=g6  25=g7  26=f7  28=g8 */
s.lambda=19397.0*    sin(                  g7              )+
        570.0*    sin(              2.0*g7              )+
       -536.0*t*  cos(                  g7              )+
        143.0*    sin(           g6-2.0*g7              )+
        110.0*t*  sin(                  g7              )+
        102.0*    sin(           g6-3.0*g7              )+
         76.0*    cos(           g6-3.0*g7              )+
        -49.0*    sin(    g5           -g7              )+
         32.0*t*t                                        +
        -30.0*t*  cos(              2.0*g7              )+
         29.0*    sin(2.0*g5-6.0*g6+3.0*g7              )+
         29.0*    cos(              2.0*g7       -2.0*g8)+
        -28.0*    cos(                  g7           -g8)+
         23.0*    sin(              3.0*g7              )+
        -21.0*    cos(    g5           -g7              )+
         20.0*    sin(                  g7           -g8)+
         20.0*    cos(           g6-2.0*g7              )+
        -19.0*    cos(           g6    -g7              )+
         17.0*    sin(              2.0*g7       -3.0*g8)+
         14.0*    sin(              3.0*g7       -3.0*g8)+
         13.0*    sin(           g6    -g7              )+
        -12.0*t*t*cos(                  g7              )+
        -12.0*    cos(                  g7              );
s.lambda+=  10.0*    sin(              2.0*g7       -2.0*g8)+
         -9.0*    sin(                     2.0*f7       )+
         -9.0*t*t*sin(                  g7              )+
          9.0*    cos(              2.0*g7       -3.0*g8)+
          8.0*t*  cos(           g6-2.0*g7              )+
          7.0*t*  cos(           g6-3.0*g7              )+
         -7.0*t*  sin(           g6-3.0*g7              )+
          7.0*t*  sin(              2.0*g7              )+
          6.0*    sin(2.0*g5-6.0*g6+2.0*g7              )+
          6.0*    cos(2.0*g5-6.0*g6+2.0*g7              )+
          5.0*    sin(           g6-4.0*g7              )+
         -4.0*    sin(              3.0*g7       -4.0*g8)+
          4.0*    cos(              3.0*g7       -3.0*g8)+
         -3.0*    cos(                                g8)+
         -2.0*    sin(                                g8);
s.lambda=mod(l7+s.lambda/648000.0*PI,M_2PI);
s.beta= 2775.0*    sin(                         f7       )+
        131.0*    sin(                  g7    -f7       )+
        130.0*    sin(                  g7    +f7       );
s.beta=s.beta/648000.0*PI;
s.r=19.21216                                           +
     -0.90154*    cos(                  g7              )+
     -0.02488*t*  sin(                  g7              )+
     -0.02121*    cos(              2.0*g7              )+
     -0.00585*    cos(           g6-2.0*g7              )+
     -0.00508*t*  cos(                  g7              )+
     -0.00451*    cos(    g5           -g7              )+
      0.00336*    sin(           g6    -g7              )+
      0.00198*    sin(    g5           -g7              )+
      0.00118*    cos(           g6-3.0*g7              )+
      0.00107*    sin(           g6-2.0*g7              )+
     -0.00103*t*  sin(              2.0*g7              )+
     -0.00081*    cos(              3.0*g7       -3.0*g8);
  return(s);
}

sphaer nep( double jd )

{ double t,d,g5,g6,g7,g8,f8,l7,l8;
  sphaer s;

  d= jd-2451545.0;
  t= d/36525.0+1.0;
  g5= M_2PI*frac(0.056531+0.00023080893*d); /* 19 */
  g6= M_2PI*frac(0.882987+0.00009294371*d); /* 22 */
  l7= M_2PI*frac(0.870169+0.00003269438*d); /* 24 */
  g7= M_2PI*frac(0.400589+0.00003269438*d); /* 25 */
  l8= M_2PI*frac(0.846912+0.00001672092*d); /* 27 */
  g8= M_2PI*frac(0.725368+0.00001672092*d); /* 28 */
  f8= M_2PI*frac(0.480856+0.00001663715*d); /* 29 */
/*                     19=g5  22=g6  24=l7  25=g7  27=l8  28=g8  29=f8 */
s.lambda= 3523.0*    sin(                                       g8       )+
        -50.0*    sin(                                          2.0*f8)+
        -43.0*t*  cos(                                       g8       )+
         29.0*    sin(    g5                                -g8       )+
         19.0*    sin(                                   2.0*g8       )+
        -18.0*    cos(    g5                                -g8       )+
         13.0*    cos(           g6                         -g8       )+
         13.0*    sin(           g6                         -g8       )+
         -9.0*    sin(                     2.0*g7       -3.0*g8       )+
          9.0*    cos(                     2.0*g7       -2.0*g8       )+
         -5.0*    cos(                     2.0*g7       -3.0*g8       )+
         -4.0*t*  sin(                                       g8       )+
          4.0*    cos(                         g7       -2.0*g8       )+
          4.0*t*t*sin(                                       g8       );
s.lambda=mod(l8+s.lambda/648000.0*PI,M_2PI);
s.beta= 6404.0*    sin(                                              f8)+
         55.0*    sin(                                       g8    +f8)+
         55.0*    sin(                                       g8    -f8)+
        -33.0*t*  sin(                                              f8);
s.beta=s.beta/648000.0*PI;
s.r=30.07175                                                         +
     -0.25701*    cos(                                       g8       )+
     -0.00787*    cos(              2.0*l7    -g7-2.0*l8              )+
      0.00409*    cos(    g5                                -g8       )+
     -0.00314*t*  sin(                                       g8       )+
      0.00250*    sin(    g5                                -g8       )+
     -0.00194*    sin(           g6                         -g8       )+
      0.00185*    cos(           g6                         -g8       );
  return(s);
}

sphaer plu( double jd )

{ double t,d,l9,g9,f9;
  sphaer s;

  d= jd-2451545.0;
  t= d/36525.0+1.0;
  l9= M_2PI*frac(0.663854+0.00001115482*d); /* 31 */
  g9= M_2PI*frac(0.041020+0.00001104864*d); /* 32 */
  f9= M_2PI*frac(0.357355+0.00001104864*d); /* 33 */
/*                     32=g9  33=f9  */
s.lambda=101577.0*   sin(    g9       )+
      15517.0*    sin(2.0*g9       )+
      -3593.0*    sin(       2.0*f9)+
       3414.0*    sin(3.0*g9       )+
      -2201.0*    sin(    g9-2.0*f9)+
      -1871.0*    sin(    g9+2.0*f9)+
        839.0*    sin(4.0*g9       )+
       -757.0*    sin(2.0*g9+2.0*f9)+
       -285.0*    sin(3.0*g9+2.0*f9)+
        227.0*t*t*sin(    g9       )+
        218.0*    sin(2.0*g9-2.0*f9)+
        200.0*t*  sin(    g9       );
s.lambda=mod(l9+s.lambda/648000.0*PI,M_2PI);
s.beta=57726.0*    sin(           f9)+
      15257.0*    sin(    g9    -f9)+
      14102.0*    sin(    g9    +f9)+
       3870.0*    sin(2.0*g9    +f9)+
       1138.0*    sin(3.0*g9    +f9)+
        472.0*    sin(2.0*g9    -f9)+
        353.0*    sin(4.0*g9    +f9)+
       -144.0*    sin(    g9-3.0*f9)+
       -119.0*    sin(       3.0*f9)+
       -111.0*    sin(    g9+3.0*f9);
s.beta=s.beta/648000.0*PI;
s.r=40.74638                      +
     -9.58235*    cos(    g9       )+
     -1.16703*    cos(2.0*g9       )+
     -0.22649*    cos(3.0*g9       )+
     -0.04996*    cos(4.0*g9       );
  return(s);
}
