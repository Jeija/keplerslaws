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
/* Module: PPOS_LO2.CPP                                                      */
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

sphaer mon( double jd )

{ double t,d,lm,gm,fm,dm,km,ls,gs,l2;
  sphaer s;

  d= jd-2451545.0;
  t= d/36525.0+1.0;
  lm= M_2PI*frac(0.606434+0.03660110129*d); /*  1 */
  gm= M_2PI*frac(0.374897+0.03629164709*d); /*  2 */
  fm= M_2PI*frac(0.259091+0.03674819520*d); /*  3 */
  dm= M_2PI*frac(0.827362+0.03386319198*d); /*  4 */
  km= M_2PI*frac(0.347343-0.00014709391*d); /*  5 */
  ls= M_2PI*frac(0.779072+0.00273790931*d); /*  7 */
  gs= M_2PI*frac(0.993126+0.00273777850*d); /*  8 */
  l2= M_2PI*frac(0.505498+0.00445046867*d); /* 12 */
/*                      2=gm   3=fm   4=dm   5=km   7=ls   8=gs   12=l2 */
s.lambda=22640.0*    sin(    gm                                           )+
      -4586.0*    sin(    gm       -2.0*dm                             )+
       2370.0*    sin(              2.0*dm                             )+
        769.0*    sin(2.0*gm                                           )+
       -668.0*    sin(                                       gs        )+
       -412.0*    sin(       2.0*fm                                    )+
       -212.0*    sin(2.0*gm       -2.0*dm                             )+
       -206.0*    sin(    gm       -2.0*dm                  +gs        )+
        192.0*    sin(    gm       +2.0*dm                             )+
        165.0*    sin(              2.0*dm                  -gs        )+
        148.0*    sin(    gm                                -gs        )+
       -125.0*    sin(                  dm                             )+
       -110.0*    sin(    gm                                +gs        )+
        -55.0*    sin(       2.0*fm-2.0*dm                             )+
        -45.0*    sin(    gm+2.0*fm                                    )+
         40.0*    sin(    gm-2.0*fm                                    )+
        -38.0*    sin(    gm       -4.0*dm                             )+
         36.0*    sin(3.0*gm                                           )+
        -31.0*    sin(2.0*gm       -4.0*dm                             )+
         28.0*    sin(    gm       -2.0*dm                  -gs        )+
        -24.0*    sin(              2.0*dm                  +gs        )+
         19.0*    sin(    gm           -dm                             )+
         18.0*    sin(                  dm                  +gs        );
s.lambda+=  15.0*    sin(    gm       +2.0*dm                  -gs        )+
         14.0*    sin(2.0*gm       +2.0*dm                             )+
         14.0*    sin(              4.0*dm                             )+
        -13.0*    sin(3.0*gm       -2.0*dm                             )+
        -11.0*    sin(    gm                    +16.0*ls       -18.0*l2)+
         10.0*    sin(2.0*gm                                -gs        )+
          9.0*    sin(    gm-2.0*fm-2.0*dm                             )+
          9.0*    cos(    gm                    +16.0*ls       -18.0*l2)+
         -9.0*    sin(2.0*gm       -2.0*dm                  +gs        )+
         -8.0*    sin(    gm           +dm                             )+
          8.0*    sin(              2.0*dm              -2.0*gs        )+
         -8.0*    sin(2.0*gm                                +gs        )+
         -7.0*    sin(                                   2.0*gs        )+
         -7.0*    sin(    gm       -2.0*dm              +2.0*gs        )+
          7.0*    sin(                        km                       )+
         -6.0*    sin(    gm-2.0*fm+2.0*dm                             )+
         -6.0*    sin(       2.0*fm+2.0*dm                             )+
         -4.0*    sin(    gm       -4.0*dm                  +gs        )+
          4.0*t*  cos(    gm                   +16.0*ls        -18.0*l2)+
         -4.0*    sin(2.0*gm+2.0*fm                                    )+
          4.0*t*  sin(    gm                   +16.0*ls        -18.0*l2)+
          3.0*    sin(    gm       -3.0*dm                             )+
         -3.0*    sin(    gm       +2.0*dm                  +gs        );
s.lambda+=  -3.0*    sin(2.0*gm       -4.0*dm                  +gs        )+
          3.0*    sin(    gm                            -2.0*gs        )+
          3.0*    sin(    gm       -2.0*dm              -2.0*gs        )+
         -2.0*    sin(2.0*gm       -2.0*dm                  -gs        )+
         -2.0*    sin(       2.0*fm-2.0*dm                  +gs        )+
          2.0*    sin(    gm       +4.0*dm                             )+
          2.0*    sin(4.0*gm                                           )+
          2.0*    sin(              4.0*dm                  -gs        )+
          2.0*    sin(2.0*gm           -dm                             );
s.lambda=mod(lm+s.lambda/648000.0*PI,M_2PI);
s.beta=18461.0*    sin(           fm                                    )+
       1010.0*    sin(    gm    +fm                                    )+
       1000.0*    sin(    gm    -fm                                    )+
       -624.0*    sin(           fm-2.0*dm                             )+
       -199.0*    sin(    gm    -fm-2.0*dm                             )+
       -167.0*    sin(    gm    +fm-2.0*dm                             )+
        117.0*    sin(           fm+2.0*dm                             )+
         62.0*    sin(2.0*gm    +fm                                    )+
         33.0*    sin(    gm    -fm+2.0*dm                             )+
         32.0*    sin(2.0*gm    -fm                                    )+
        -30.0*    sin(           fm-2.0*dm                  +gs        )+
        -16.0*    sin(2.0*gm    +fm-2.0*dm                             )+
         15.0*    sin(    gm    +fm+2.0*dm                             )+
         12.0*    sin(           fm-2.0*dm                  -gs        )+
         -9.0*    sin(    gm    -fm-2.0*dm                  +gs        )+
         -8.0*    sin(           fm           +km                      )+
          8.0*    sin(           fm+2.0*dm                  -gs        )+
         -7.0*    sin(    gm    +fm-2.0*dm                  +gs        )+
          7.0*    sin(    gm    +fm                         -gs        )+
         -7.0*    sin(    gm    +fm-4.0*dm                             )+
         -6.0*    sin(          +fm                         +gs        )+
         -6.0*    sin(       3.0*fm                                    )+
          6.0*    sin(    gm    -fm                         -gs        );
s.beta+=  -5.0*    sin(           fm    +dm                             )+
         -5.0*    sin(    gm    +fm                         +gs        )+
         -5.0*    sin(    gm    -fm                         +gs        )+
          5.0*    sin(           fm                         -gs        )+
          5.0*    sin(           fm    -dm                             )+
          4.0*    sin(3.0*gm    +fm                                    )+
         -4.0*    sin(           fm-4.0*dm                             )+
         -3.0*    sin(    gm    -fm-4.0*dm                             )+
          3.0*    sin(    gm-3.0*fm                                    )+
         -2.0*    sin(2.0*gm    -fm-4.0*dm                             )+
         -2.0*    sin(       3.0*fm-2.0*dm                             )+
          2.0*    sin(2.0*gm    -fm+2.0*dm                             )+
          2.0*    sin(    gm    -fm+2.0*dm                  -gs        )+
          2.0*    sin(2.0*gm    -fm-2.0*dm                             )+
          2.0*    sin(3.0*gm    -fm                                    );
s.beta=s.beta/648000.0*PI;
s.r=60.36298                                                          +
     -3.27746*    cos(    gm                                           )+
     -0.57994*    cos(    gm       -2.0*dm                             )+
     -0.46357*    cos(              2.0*dm                             )+
     -0.08904*    cos(2.0*gm                                           )+
      0.03865*    cos(2.0*gm       -2.0*dm                             )+
     -0.03237*    cos(              2.0*dm                  -gs        )+
     -0.02688*    cos(    gm       +2.0*dm                             )+
     -0.02358*    cos(    gm       -2.0*dm                  +gs        )+
     -0.02030*    cos(    gm                                -gs        )+
      0.01719*    cos(                  dm                             )+
      0.01671*    cos(    gm                                +gs        )+
      0.01247*    cos(    gm-2.0*fm                                    )+
      0.00704*    cos(                                      +gs        )+
      0.00529*    cos(              2.0*dm                  +gs        )+
     -0.00524*    cos(    gm       -4.0*dm                             )+
      0.00398*    cos(    gm       -2.0*dm                  -gs        )+
     -0.00366*    cos(3.0*gm                                           )+
     -0.00295*    cos(2.0*gm       -4.0*dm                             )+
     -0.00263*    cos(                  dm                  +gs        )+
      0.00249*    cos(3.0*gm       -2.0*dm                             )+
     -0.00221*    cos(    gm       +2.0*dm                  -gs        )+
      0.00185*    cos(       2.0*fm-2.0*dm                             );
s.r+=
     -0.00161*    cos(              2.0*dm              -2.0*gs        )+
      0.00147*    cos(    gm+2.0*fm-2.0*dm                             )+
     -0.00142*    cos(              4.0*dm                             )+
      0.00139*    cos(2.0*gm       -2.0*dm                  +gs        )+
     -0.00118*    cos(    gm       -4.0*dm                  +gs        )+
     -0.00116*    cos(2.0*gm       +2.0*dm                             )+
     -0.00110*    cos(2.0*gm                                -gs        );
  return(s);
}

sphaer mar( double jd )

{ double t,d,gs,g2,l4,g4,f4,g5;
  sphaer s;

  d= jd-2451545.0;
  t= d/36525.0+1.0;
  gs= M_2PI*frac(0.993126+0.00273777850*d); /*  8 */
  g2= M_2PI*frac(0.140023+0.00445036173*d); /* 13 */
  l4= M_2PI*frac(0.987353+0.00145575328*d); /* 15 */
  g4= M_2PI*frac(0.053856+0.00145561327*d); /* 16 */
  f4= M_2PI*frac(0.849694+0.00145569465*d); /* 17 */
  g5= M_2PI*frac(0.056531+0.00023080893*d); /* 19 */
/*                      8=gs  13=g2  16=g4  17=f4  19=g5  */
s.lambda=38451.0*    sin(                  g4              )+
       2238.0*    sin(              2.0*g4              )+
        181.0*    sin(              3.0*g4              )+
        -52.0*    sin(                     2.0*f4       )+
         37.0*t*  sin(                  g4              )+
        -22.0*    cos(                  g4       -2.0*g5)+
        -19.0*    sin(                  g4           -g5)+
         17.0*    cos(                  g4           -g5)+
         17.0*    sin(              4.0*g4              )+
        -16.0*    cos(              2.0*g4       -2.0*g5)+
         13.0*    cos(    gs       -2.0*g4              )+
        -10.0*    sin(                  g4-2.0*f4       )+
        -10.0*    sin(                  g4+2.0*f4       )+
          7.0*    cos(    gs           -g4              )+
         -7.0*    cos(2.0*gs       -3.0*g4              )+
         -5.0*    sin(           g2-3.0*g4              )+
         -5.0*    sin(    gs           -g4              )+
         -5.0*    sin(    gs       -2.0*g4              )+
         -4.0*    cos(2.0*gs       -4.0*g4              )+
          4.0*t*  sin(              2.0*g4              )+
          4.0*    cos(                                g5)+
          3.0*    cos(           g2-3.0*g4              )+
          3.0*    sin(              2.0*g4       -2.0*g5);
s.lambda=mod(l4+s.lambda/648000.0*PI,M_2PI);
s.beta= 6603.0*    sin(                         f4       )+
        622.0*    sin(                  g4    -f4       )+
        615.0*    sin(                  g4    +f4       )+
         64.0*    sin(              2.0*g4    +f4       );
s.beta=s.beta/648000.0*PI;
s.r=1.53031                                            +
     -0.14170*    cos(                  g4              )+
     -0.00660*    cos(              2.0*g4              )+
     -0.00047*    cos(              3.0*g4              );
  return(s);
}
