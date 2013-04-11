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
/* Module: PPOS_LO3.CPP                                                      */
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

sphaer jup( double jd )

{ double t,d,g5,g6,g7,l5;
  sphaer s;

  d= jd-2451545.0;
  t= d/36525.0+1.0;
  l5= M_2PI*frac(0.089608+0.00023080893*d); /* 18 */
  g5= M_2PI*frac(0.056531+0.00023080893*d); /* 19 */
  g6= M_2PI*frac(0.882987+0.00009294371*d); /* 22 */
  g7= M_2PI*frac(0.400589+0.00003269438*d); /* 25 */

/*                     18=l5  19=g5  22=g6  25=g7  */
s.lambda=19934.0*    sin(           g5              )+
       5023.0*t                                   +
       2511.0                                     +
       1093.0*    cos(       2.0*g5-5.0*g6       )+
        601.0*    sin(       2.0*g5              )+
       -479.0*    sin(       2.0*g5-5.0*g6       )+
       -185.0*    sin(       2.0*g5-2.0*g6       )+
        137.0*    sin(       3.0*g5-5.0*g6       )+
       -131.0*    sin(           g5-2.0*g6       )+
         79.0*    cos(           g5    -g6       )+
        -76.0*    cos(       2.0*g5-2.0*g6       )+
        -74.0*t*  cos(           g5              )+
         68.0*t*  sin(           g5              )+
         66.0*    cos(       2.0*g5-3.0*g6       )+
         63.0*    cos(       3.0*g5-5.0*g6       )+
         53.0*    cos(           g5-5.0*g6       )+
         49.0*    sin(       2.0*g5-3.0*g6       )+
        -43.0*t*  sin(       2.0*g5-5.0*g6       )+
        -37.0*    cos(           g5              )+
         25.0*    sin(2.0*l5                     )+
         25.0*    sin(       3.0*g5              )+
        -23.0*    sin(           g5-5.0*g6       )+
        -19.0*t*  cos(       2.0*g5-5.0*g6       );
s.lambda+=  17.0*    cos(       2.0*g5-4.0*g6       )+
         17.0*    cos(       3.0*g5-3.0*g6       )+
        -14.0*    sin(           g5    -g6       )+
        -13.0*    sin(       3.0*g5-4.0*g6       )+
         -9.0*    cos(2.0*l5                     )+
          9.0*    cos(                  g6       )+
         -9.0*    sin(                  g6       )+
         -9.0*    sin(       3.0*g5-2.0*g6       )+
          9.0*    sin(       4.0*g5-5.0*g6       )+
          9.0*    sin(       2.0*g5-6.0*g6+3.0*g7)+
         -8.0*    cos(      4.0*g5-10.0*g6       )+
          7.0*    cos(       3.0*g5-4.0*g6       )+
         -7.0*    cos(           g5-3.0*g6       )+
         -7.0*    sin(      4.0*g5-10.0*g6       )+
         -7.0*    sin(           g5-3.0*g6       )+
          6.0*    cos(       4.0*g5-5.0*g6       )+
         -6.0*    sin(       3.0*g5-3.0*g6       )+
          5.0*    cos(              2.0*g6       )+
         -4.0*    sin(       4.0*g5-4.0*g6       )+
         -4.0*    cos(              3.0*g6       )+
          4.0*    cos(       2.0*g5    -g6       )+
         -4.0*    cos(       3.0*g5-2.0*g6       )+
         -4.0*t*  cos(       2.0*g5              );
s.lambda+=   3.0*t*  sin(       2.0*g5              )+
          3.0*    cos(              5.0*g6       )+
          3.0*    cos(      5.0*g5-10.0*g6       )+
          3.0*    sin(              2.0*g6       )+
         -2.0*    sin(2.0*l5    -g5              )+
          2.0*    sin(2.0*l5    +g5              )+
         -2.0*t*  sin(       3.0*g5-5.0*g6       )+
         -2.0*t*  sin(           g5-5.0*g6       );
s.lambda=mod(l5+s.lambda/648000.0*PI,M_2PI);
s.beta=-4692.0*    cos(           g5              )+
        259.0*    sin(           g5              )+
        227.0                                     +
       -227.0*    cos(       2.0*g5              )+
         30.0*t*  sin(           g5              )+
         21.0*t*  cos(           g5              )+
         16.0*    sin(       3.0*g5-5.0*g6       )+
        -13.0*    sin(           g5-5.0*g6       )+
        -12.0*    cos(       3.0*g5              )+
         12.0*    sin(       2.0*g5              )+
          7.0*    cos(       3.0*g5-5.0*g6       )+
         -5.0*    cos(           g5-5.0*g6       );
s.beta=s.beta/648000.0*PI;
s.r=5.20883                                     +
     -0.25122*    cos(           g5              )+
     -0.00604*    cos(       2.0*g5              )+
      0.00260*    cos(       2.0*g5-2.0*g6       )+
     -0.00170*    cos(       3.0*g5-5.0*g6       )+
     -0.00106*    sin(       2.0*g5-2.0*g6       )+
     -0.00091*t*  sin(           g5              )+
     -0.00084*t*  cos(           g5              )+
      0.00069*    sin(       2.0*g5-3.0*g6       )+
     -0.00067*    sin(           g5-5.0*g6       )+
      0.00066*    sin(       3.0*g5-5.0*g6       )+
      0.00063*    sin(           g5    -g6       )+
     -0.00051*    cos(       2.0*g5-3.0*g6       )+
     -0.00046*    sin(           g5              )+
     -0.00029*    cos(           g5-5.0*g6       )+
      0.00027*    cos(           g5-2.0*g6       )+
     -0.00022*    cos(       3.0*g5              )+
     -0.00021*    sin(       2.0*g5-5.0*g6       );
  return(s);
}

sphaer sat( double jd )

{ double t,d,g5,g6,l6,g7;
  sphaer s;

  d= jd-2451545.0;
  t= d/36525.0+1.0;
  g5= M_2PI*frac(0.056531+0.00023080893*d); /* 19 */
  l6= M_2PI*frac(0.133295+0.00009294371*d); /* 21 */
  g6= M_2PI*frac(0.882987+0.00009294371*d); /* 22 */
  g7= M_2PI*frac(0.400589+0.00003269438*d); /* 25 */
/*                     19=g5  21=l6  22=g6  25=g7  */
s.lambda=23045.0*    sin(                  g6       )+
       5014.0*t                                   +
      -2689.0*    cos(2.0*g5       -5.0*g6       )+
       2507.0                                     +
       1177.0*    sin(2.0*g5       -5.0*g6       )+
       -826.0*    cos(2.0*g5       -4.0*g6       )+
        802.0*    sin(              2.0*g6       )+
        425.0*    sin(    g5       -2.0*g6       )+
       -229.0*t*  cos(                  g6       )+
       -153.0*    cos(2.0*g5       -6.0*g6       )+
       -142.0*t*  sin(                  g6       )+
       -114.0*    cos(                  g6       )+
        101.0*t*  sin(2.0*g5       -5.0*g6       )+
        -70.0*    cos(       2.0*l6              )+
         67.0*    sin(       2.0*l6              )+
         66.0*    sin(2.0*g5       -6.0*g6       )+
         60.0*t*  cos(2.0*g5       -5.0*g6       )+
         41.0*    sin(    g5       -3.0*g6       )+
         39.0*    sin(              3.0*g6       )+
         31.0*    sin(    g5           -g6       )+
         31.0*    sin(2.0*g5       -2.0*g6       )+
        -29.0*    cos(2.0*g5       -3.0*g6       )+
        -28.0*    sin(2.0*g5       -6.0*g6+3.0*g7);
s.lambda+=  28.0*    cos(    g5       -3.0*g6       )+
         22.0*t*  sin(2.0*g5       -4.0*g6       )+
        -22.0*    sin(                  g6-3.0*g7)+
         20.0*    sin(2.0*g5       -3.0*g6       )+
         20.0*    cos(4.0*g5      -10.0*g6       )+
         19.0*    cos(              2.0*g6-3.0*g7)+
         19.0*    sin(4.0*g5      -10.0*g6       )+
        -17.0*t*  cos(              2.0*g6       )+
        -16.0*    cos(                  g6-3.0*g7)+
        -12.0*    sin(2.0*g5       -4.0*g6       )+
         12.0*    cos(    g5                     )+
        -12.0*    sin(              2.0*g6-2.0*g7)+
        -11.0*t*  sin(              2.0*g6       )+
        -11.0*    cos(2.0*g5       -7.0*g6       )+
         10.0*    sin(              2.0*g6-3.0*g7)+
         10.0*    cos(2.0*g5       -2.0*g6       )+
          9.0*    sin(4.0*g5       -9.0*g6       )+
         -8.0*    sin(                  g6-2.0*g7)+
         -8.0*    cos(       2.0*l6    +g6       )+
          8.0*    cos(       2.0*l6    -g6       )+
          8.0*    cos(                  g6    -g7)+
         -8.0*    sin(       2.0*l6    -g6       )+
          7.0*    sin(       2.0*l6    +g6       );
s.lambda+=  -7.0*    cos(    g5       -2.0*g6       )+
         -7.0*    cos(              2.0*g6       )+
         -6.0*t*  sin(4.0*g5      -10.0*g6       )+
          6.0*t*  cos(4.0*g5      -10.0*g6       )+
	  6.0*t*  sin(2.0*g5       -6.0*g6       )+
         -5.0*    sin(3.0*g5       -7.0*g6       )+
         -5.0*    cos(3.0*g5       -3.0*g6       )+
         -5.0*    cos(              2.0*g6-2.0*g7)+
          5.0*    sin(3.0*g5       -4.0*g6       )+
          5.0*    sin(2.0*g5       -7.0*g6       )+
          4.0*    sin(3.0*g5       -3.0*g6       )+
          4.0*    sin(3.0*g5       -5.0*g6       )+
          4.0*t*  cos(    g5       -2.0*g6       )+
          3.0*t*  cos(2.0*g5       -4.0*g6       )+
          3.0*    cos(2.0*g5       -6.0*g6+3.0*g7)+
         -3.0*t*  sin(       2.0*l6              )+
          3.0*t*  cos(2.0*g5       -6.0*g6       )+
         -3.0*t*  cos(       2.0*l6              )+
          3.0*    cos(3.0*g5       -7.0*g6       )+
          3.0*    cos(4.0*g5       -9.0*g6       )+
          3.0*    sin(3.0*g5       -6.0*g6       )+
          3.0*    sin(2.0*g5           -g6       )+
          3.0*    sin(    g5       -4.0*g6       );
s.lambda+=   2.0*    cos(              3.0*g6-3.0*g7)+
          2.0*t*  sin(    g5       -2.0*g6       )+
          2.0*    sin(              4.0*g6       )+
         -2.0*    cos(3.0*g5       -4.0*g6       )+
         -2.0*    cos(2.0*g5           -g6       )+
         -2.0*    sin(2.0*g5       -7.0*g6+3.0*g7)+
          2.0*    cos(    g5       -4.0*g6       )+
          2.0*    cos(4.0*g5      -11.0*g6       )+
         -2.0*    sin(                  g6    -g7);
s.lambda=mod(l6+s.lambda/648000.0*PI,M_2PI);
s.beta= 8297.0*    sin(                  g6       )+
      -3346.0*    cos(                  g6       )+
        462.0*    sin(              2.0*g6       )+
       -189.0*    cos(              2.0*g6       )+
        185.0                                     +
         79.0*t*  cos(                  g6       )+
        -71.0*    cos(2.0*g5       -4.0*g6       )+
         46.0*    sin(2.0*g5       -6.0*g6       )+
        -45.0*    cos(2.0*g5       -6.0*g6       )+
         29.0*    sin(              3.0*g6       )+
        -20.0*    cos(2.0*g5       -3.0*g6       )+
         18.0*t*  sin(                  g6       )+
        -14.0*    cos(2.0*g5       -5.0*g6       )+
        -11.0*    cos(              3.0*g6       )+
        -10.0*t                                   +
          9.0*    sin(    g5       -3.0*g6       )+
          8.0*    sin(    g5           -g6       )+
         -6.0*    sin(2.0*g5       -3.0*g6       )+
          5.0*    sin(2.0*g5       -7.0*g6       )+
         -5.0*    cos(2.0*g5       -7.0*g6       )+
          4.0*    sin(2.0*g5       -5.0*g6       )+
         -4.0*t*  sin(              2.0*g6       )+
         -3.0*    cos(    g5           -g6       )+
          3.0*    cos(    g5       -3.0*g6       )+
          3.0*t*  sin(2.0*g5       -4.0*g6       )+
          3.0*    sin(    g5       -2.0*g6       )+
          2.0*    sin(              4.0*g6       )+
         -2.0*    cos(2.0*g5       -2.0*g6       );
s.beta=s.beta/648000.0*PI;
s.r=9.55774                                     +
     -0.53252*    cos(                  g6       )+
     -0.01878*    sin(2.0*g5       -4.0*g6       )+
     -0.01482*    cos(              2.0*g6       )+
      0.00817*    sin(    g5           -g6       )+
     -0.00539*    cos(    g5       -2.0*g6       )+
     -0.00524*t*  sin(                  g6       )+
      0.00349*    sin(2.0*g5       -5.0*g6       )+
      0.00347*    sin(2.0*g5       -6.0*g6       )+
      0.00328*t*  cos(                  g6       )+
     -0.00225*    sin(                  g6       )+
      0.00149*    cos(2.0*g5       -6.0*g6       )+
     -0.00126*    cos(2.0*g5       -2.0*g6       )+
      0.00104*    cos(    g5           -g6       )+
      0.00101*    cos(2.0*g5       -5.0*g6       )+
      0.00098*    cos(    g5       -3.0*g6       )+
     -0.00073*    cos(2.0*g5       -3.0*g6       )+
     -0.00062*    cos(              3.0*g6       )+
      0.00042*    sin(              2.0*g6-3.0*g7)+
      0.00041*    sin(2.0*g5       -2.0*g6       )+
     -0.00040*    sin(    g5       -3.0*g6       )+
      0.00040*    cos(2.0*g5       -4.0*g6       )+
     -0.00028*t                                   +
     -0.00023*    sin(    g5                     )+
      0.00020*    sin(2.0*g5       -7.0*g6       );
  return(s);
}
