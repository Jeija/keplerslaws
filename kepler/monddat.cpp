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
/* Module: MONDDAT.CPP                                                       */
/* Version 2.0                                                               */
/* Last modified: Jan 16, 1994                                               */
/*                                                                           */
/* Version 2.0: Stoerungen durch Planeten werden beruecksichtigt             */
/*****************************************************************************/
/* Abkuerzungen (Literatur):                                                 */
/*                                                                           */
/* GdE   : Montenbruck: Grundlagen der Ephemeridenberechnung                 */
/* AmdPC : Montenbruck,Pfleger: Astronomie mit dem Personal Computer         */
/* PAwyPC: Duffet-Smith: Practical Astronomy with your Personal Computer     */
/* AAyyyy: Astronomical Almanac of the Year yyyy                             */
/* AsAl  : Meeus: Astronomical Algorithms                                    */
/*****************************************************************************/

#include <math.h>
#include "astromat.h"
#include "monddat.h"

/**************************/
/* Funktionen             */
/**************************/
/* mond_per_apo()         */
/* mond_perigaeum()       */
/* mond_apogaeum()        */
/* mond_mittlere_phase()  */
/* mond_neu               */
/* mond_erstes_viertel()  */
/* mond_voll()            */
/* mond_letztes_viertel() */
/* mond_alter()           */
/* mond_apsiden_phasen()  */
/**************************/

double mond_per_apo( double k, double T )
{  return 2451534.6698+27.55454988*k
		      - 0.0006886*T*T
		      - 0.000001098*T*T*T
		      + 0.0000000052*T*T*T*T;
}

double mond_perigaeum( double jd_start )
{  double k,T,jd,D,F,M,jdk;

   k=ceil((jd_start-PER_MON)/ANO_MON);
   T=k/1325.55;
   jd=mond_per_apo(k,T);
   if (jd<jd_start)
   {  k=k+1.0;
      T=k/1325.55;
      jd=mond_per_apo(k,T);
   }

   D=171.9179*DEG2RAD+335.9106046*DEG2RAD*k
		     -  0.0100250*DEG2RAD*T*T
		     -  0.00001156*DEG2RAD*T*T*T
		     +  0.000000055*DEG2RAD*T*T*T*T;
   F=316.6109*DEG2RAD+364.5287911*DEG2RAD*k
		     -  0.0125131*DEG2RAD*T*T
		     -  0.0000148*DEG2RAD*T*T*T;
   M=347.3477*DEG2RAD+ 27.1577721*DEG2RAD*k
		     -  0.0008323*DEG2RAD*T*T
		     -  0.0000010*DEG2RAD*T*T*T;

   jdk=	/*   1 */ -1.6769*    sin(  +2*D                  )
	/*   2 */ +0.4589*    sin(  +4*D                  )
	/*   3 */ -0.1856*    sin(  +6*D                  )
	/*   4 */ +0.0883*    sin(  +8*D                  )
	/*   5 */ +(-0.0773+0.00019*T)*    sin(  +2*D              -M  )
	/*   6 */ +(+0.0502-0.00013*T)*    sin(                    +M  )
	/*   7 */ -0.0460*    sin( +10*D                  )
	/*   8 */ +(+0.0422-0.00011*T)*    sin(  +4*D              -M  )
	/*   9 */ -0.0256*    sin(  +6*D              -M  )
	/*  10 */ +0.0253*    sin( +12*D                  )
	/*  11 */ +0.0237*    sin(    +D                  )
	/*  12 */ +0.0162*    sin(  +8*D              -M  )
	/*  13 */ -0.0145*    sin( +14*D                  )
	/*  14 */ +0.0129*    sin(          +2*F          )
	/*  15 */ -0.0112*    sin(  +3*D                  )
	/*  16 */ -0.0104*    sin( +10*D              -M  )
	/*  17 */ +0.0086*    sin( +16*D                  )
	/*  18 */ +0.0069*    sin( +12*D              -M  )
	/*  19 */ +0.0066*    sin(  +5*D                  )
	/*  20 */ -0.0053*    sin(  +2*D    +2*F          );

   jdk+=/* 21 */  -0.0052*    sin( +18*D                  )
	/*  22 */ -0.0046*    sin( +14*D              -M  )
	/*  23 */ -0.0041*    sin(  +7*D                  )
	/*  24 */ +0.0040*    sin(  +2*D              +M  )
	/*  25 */ +0.0032*    sin( +20*D                  )
	/*  26 */ -0.0032*    sin(    +D              +M  )
	/*  27 */ +0.0031*    sin( +16*D              -M  )
	/*  28 */ -0.0029*    sin(  +4*D              +M  )
	/*  29 */ +0.0027*    sin(  +9*D                  )
	/*  30 */ +0.0027*    sin(  +4*D    +2*F          )
	/*  31 */ -0.0027*    sin(  +2*D            -2*M  )
	/*  32 */ +0.0024*    sin(  +4*D            -2*M  )
	/*  33 */ -0.0021*    sin(  +6*D            -2*M  )
	/*  34 */ -0.0021*    sin( +22*D                  )
	/*  35 */ -0.0021*    sin( +18*D              -M  )
	/*  36 */ +0.0019*    sin(  +6*D              +M  )
	/*  37 */ -0.0018*    sin( +11*D                  )
	/*  38 */ -0.0014*    sin(  +8*D              +M  )
	/*  39 */ -0.0014*    sin(  +4*D    -2*F          )
	/*  40 */ -0.0014*    sin(  +6*D    +2*F          );

   jdk+=/*  41 */ +0.0014*    sin(  +3*D              +M  )
	/*  42 */ -0.0014*    sin(  +5*D              +M  )
	/*  43 */ +0.0013*    sin( +13*D                  )
	/*  44 */ +0.0013*    sin( +20*D              -M  )
	/*  45 */ +0.0011*    sin(  +3*D            +2*M  )
	/*  46 */ -0.0011*    sin(  +4*D    +2*F    -2*M  )
	/*  47 */ -0.0010*    sin(    +D            +2*M  )
	/*  48 */ -0.0009*    sin( +22*D              -M  )
	/*  49 */ -0.0008*    sin(          +4*F          )
	/*  50 */ +0.0008*    sin(  +6*D    -2*F          )
	/*  51 */ +0.0008*    sin(  +2*D    -2*F      +M  )
	/*  52 */ +0.0007*    sin(                  +2*M  )
	/*  53 */ +0.0007*    sin(          +2*F      -M  )
	/*  54 */ +0.0007*    sin(  +2*D    +4*F          )
	/*  55 */ -0.0006*    sin(          +2*F    -2*M  )
	/*  56 */ -0.0006*    sin(  +2*D    -2*F    +2*M  )
	/*  57 */ +0.0006*    sin( +24*D                  )
	/*  58 */ +0.0005*    sin(  +4*D    -4*F          )
	/*  59 */ +0.0005*    sin(  +2*D            +2*M  )
	/*  60 */ -0.0004*    sin(    +D              -M  );

   return jd+jdk;
}

double mond_apogaeum( double jd_start )
{  double k,T,jd,D,F,M,jdk;

   k=ceil((jd_start-PER_MON)/ANO_MON-0.5)+0.5;
   T=k/1325.55;
   jd=mond_per_apo(k,T);
   if (jd<jd_start)
   {  k=k+1.0;
      T=k/1325.55;
      jd=mond_per_apo(k,T);
   }
   D=171.9179*DEG2RAD+335.9106046*DEG2RAD*k
		     -  0.0100250*DEG2RAD*T*T
		     -  0.00001156*DEG2RAD*T*T*T
		     +  0.000000055*DEG2RAD*T*T*T*T;
   F=316.6109*DEG2RAD+364.5287911*DEG2RAD*k
		     -  0.0125131*DEG2RAD*T*T
		     -  0.0000148*DEG2RAD*T*T*T;
   M=347.3477*DEG2RAD+ 27.1577721*DEG2RAD*k
		     -  0.0008323*DEG2RAD*T*T
		     -  0.0000010*DEG2RAD*T*T*T;

   jdk= /*   1 */ +0.4392*    sin(  +2*D                  )
	/*   2 */ +0.0684*    sin(  +4*D                  )
	/*   3 */ +(+0.0456-0.00011*T)*  sin(                    +M  )
	/*   4 */ +(+0.0426-0.00011*T)*  sin(  +2*D              -M  )
	/*   5 */ +0.0212*    sin(          +2*F          )
	/*   6 */ -0.0189*    sin(    +D                  )
	/*   7 */ +0.0144*    sin(  +6*D                  )
	/*   8 */ +0.0113*    sin(  +4*D              -M  )
	/*   9 */ +0.0047*    sin(  +2*D    +2*F          )
	/*  10 */ +0.0036*    sin(    +D              +M  )
	/*  11 */ +0.0035*    sin(  +8*D                  )
	/*  12 */ +0.0034*    sin(  +6*D              -M  )
	/*  13 */ -0.0034*    sin(  +2*D    -2*F          )
	/*  14 */ +0.0022*    sin(  +2*D            -2*M  )
	/*  15 */ -0.0017*    sin(  +3*D                  )
	/*  16 */ +0.0013*    sin(  +4*D    +2*F          )
	/*  17 */ +0.0011*    sin(  +8*D              -M  )
	/*  18 */ +0.0010*    sin(  +4*D            -2*M  )
	/*  19 */ +0.0009*    sin( +10*D                  )
	/*  20 */ +0.0007*    sin(  +3*D              +M  );

   jdk+=/*  21 */ +0.0006*    sin(                  +2*M  )
	/*  22 */ +0.0005*    sin(  +2*D              +M  )
	/*  23 */ +0.0005*    sin(  +2*D            +2*M  )
	/*  24 */ +0.0004*    sin(  +6*D    +2*F          )
	/*  25 */ +0.0004*    sin(  +6*D            -2*M  )
	/*  26 */ +0.0004*    sin( +10*D              -M  )
	/*  27 */ -0.0004*    sin(  +5*D                  )
	/*  28 */ -0.0004*    sin(  +4*D    -2*F          )
	/*  29 */ +0.0003*    sin(          +2*F      +M  )
	/*  30 */ +0.0003*    sin( +12*D                  )
	/*  31 */ +0.0003*    sin(  +2*D    +2*F      -M  )
	/*  32 */ -0.0003*    sin(    +D              -M  );

   return jd+jdk;
}

double mond_mittlere_phase( double k, double T )
{  return 2451550.09765+29.530588853*k
		       + 0.0001337*T*T
		       - 0.000000150*T*T*T
		       + 0.00000000073*T*T*T*T;
}

double mond_neu( double jd_start )
/* Neumond berechnen.
   Eingabe:
   jd_start: JD NACH dem der Neumond folgen soll
   Ausgabe:
   JDE des Neumondzeitpunktes
*/
{  double k,T,mm,fm,ms,km,jd,jdk,E;

   k=ceil((jd_start-NEU_MON)/SYN_MON);
   T=k/1236.85;
   jd=mond_mittlere_phase(k,T);

   ms=  2.5534*DEG2RAD+ 29.10535669*DEG2RAD*k
		      -  0.0000218 *DEG2RAD*T*T
		      -  0.00000011*DEG2RAD*T*T*T;
   mm=201.5643*DEG2RAD+385.81693528*DEG2RAD*k
		      +  0.0107438 *DEG2RAD*T*T
		      +  0.00001239*DEG2RAD*T*T*T
		      -  0.000000058*DEG2RAD*T*T*T*T;
   fm=160.7108*DEG2RAD+390.67050274*DEG2RAD*k
		      -  0.0016341 *DEG2RAD*T*T
		      -  0.00000227*DEG2RAD*T*T*T
		      +  0.000000011*DEG2RAD*T*T*T*T;
   km=124.7746*DEG2RAD-  1.56375580*DEG2RAD*k
		      +  0.0020691 *DEG2RAD*T*T
		      +  0.00000215*DEG2RAD*T*T*T;
   E=1.0-0.002516*T-0.0000074*T*T;

   /*                        mm     ms     fm km   */
   jdk=-40720.0    *sin(     mm                 )
       +17241.0*E  *sin(            ms          )
       + 1608.0    *sin(+2.0*mm                 )
       + 1039.0    *sin(              +2.0*fm   )
       +  739.0*E  *sin(     mm    -ms          )
       -  514.0*E  *sin(     mm    +ms          )
       +  208.0*E*E*sin(       +2.0*ms          )
       -  111.0    *sin(     mm       -2.0*fm   )
       -   57.0    *sin(     mm       +2.0*fm   )
       +   56.0*E  *sin(+2.0*mm    +ms          )
       -   42.0    *sin(+3.0*mm                 )
       +   42.0*E  *sin(            ms+2.0*fm   )
       +   38.0*E  *sin(            ms-2.0*fm   )
       -   24.0*E  *sin(+2.0*mm    -ms          )
       -   17.0    *sin(                     +km)
       -    7.0    *sin(     mm+2.0*ms          )
       +    4.0    *sin(+2.0*mm       -2.0*fm   )
       +    4.0    *sin(       +3.0*ms          )
       +    3.0    *sin(     mm    +ms-2.0*fm   )
       +    3.0    *sin(+2.0*mm       +2.0*fm   )
       -    3.0    *sin(     mm    +ms+2.0*fm   )
       +    3.0    *sin(     mm    -ms+2.0*fm   )
       -    2.0    *sin(     mm    -ms-2.0*fm   )
       -    2.0    *sin(+3.0*mm    +ms          )
       +    2.0    *sin(+4.0*mm                 );

   jdk+=  +32.5    *sin( 299.77*DEG2RAD +  0.107408*DEG2RAD*k
                                        -  0.009173*DEG2RAD*T*T )
          +16.5    *sin( 251.88*DEG2RAD +  0.016321*DEG2RAD*k )
          +16.4    *sin( 251.83*DEG2RAD + 26.651886*DEG2RAD*k )
          +12.6    *sin( 349.42*DEG2RAD + 36.412478*DEG2RAD*k )
          +11.0    *sin(  84.66*DEG2RAD + 18.206239*DEG2RAD*k )
           +6.2    *sin( 141.74*DEG2RAD + 53.303771*DEG2RAD*k )
           +6.0    *sin( 207.14*DEG2RAD +  2.453732*DEG2RAD*k )
           +5.6    *sin( 154.84*DEG2RAD +  7.306860*DEG2RAD*k )
           +4.7    *sin(  34.52*DEG2RAD + 27.261239*DEG2RAD*k )
           +4.2    *sin( 207.19*DEG2RAD +  0.121824*DEG2RAD*k )
           +4.0    *sin( 291.34*DEG2RAD +  1.844379*DEG2RAD*k )
           +3.7    *sin( 161.72*DEG2RAD + 24.198154*DEG2RAD*k )
           +3.5    *sin( 239.56*DEG2RAD + 25.513099*DEG2RAD*k )
           +2.3    *sin( 331.55*DEG2RAD +  3.592518*DEG2RAD*k );

   return jd+jdk/100000.0;
}

double mond_erstes_viertel( double jd_start )
/* Erstes Viertel berechnen.
   Eingabe:
   jd_start: JD NACH dem das 1. 1/4 folgen soll
   Ausgabe:
   JDE des 1. 1/4
*/
{  double k,T,mm,fm,ms,km,jd,w,jdk,E;

   k=ceil((jd_start-NEU_MON)/SYN_MON-0.25)+0.25;
   T=k/1236.85;
   jd=mond_mittlere_phase(k,T);

   ms=  2.5534*DEG2RAD+ 29.10535669*DEG2RAD*k
		      -  0.0000218 *DEG2RAD*T*T
		      -  0.00000011*DEG2RAD*T*T*T;
   mm=201.5643*DEG2RAD+385.81693528*DEG2RAD*k
		      +  0.0107438 *DEG2RAD*T*T
		      +  0.00001239*DEG2RAD*T*T*T
		      -  0.000000058*DEG2RAD*T*T*T*T;
   fm=160.7108*DEG2RAD+390.67050274*DEG2RAD*k
		      -  0.0016341 *DEG2RAD*T*T
		      -  0.00000227*DEG2RAD*T*T*T
		      +  0.000000011*DEG2RAD*T*T*T*T;
   km=124.7746*DEG2RAD-  1.56375580*DEG2RAD*k
		      +  0.0020691 *DEG2RAD*T*T
		      +  0.00000215*DEG2RAD*T*T*T;
   E=1.0-0.002516*T-0.0000074*T*T;

   /*                        mm     ms     fm km   */
   jdk=-62801.0    *sin(     mm                 )
       +17172.0*E  *sin(            ms          )
       - 1183.0*E  *sin(     mm    +ms          )
       +  862.0    *sin(+2.0*mm                 )
       +  804.0    *sin(              +2.0*fm   )
       +  454.0*E  *sin(     mm    -ms          )
       +  204.0*E*E*sin(       +2.0*ms          )
       -  180.0    *sin(     mm       -2.0*fm   )
       -   70.0    *sin(     mm       +2.0*fm   )
       -   40.0    *sin(+3.0*mm                 )
       -   34.0*E  *sin(+2.0*mm    -ms          )
       +   32.0*E  *sin(            ms+2.0*fm   )
       +   32.0*E  *sin(            ms-2.0*fm   )
       -   28.0*E*E*sin(     mm+2.0*ms          )
       +   27.0*E  *sin(+2.0*mm    +ms          )
       -   17.0    *sin(                     +km)
       -    5.0    *sin(     mm    -ms-2.0*fm   )
       +    4.0    *sin(+2.0*mm       +2.0*fm   )
       -    4.0    *sin(     mm    +ms+2.0*fm   )
       +    4.0    *sin(     mm-2.0*ms          )
       +    3.0    *sin(     mm    +ms-2.0*fm   )
       +    3.0    *sin(       +3.0*ms          )
       +    2.0    *sin(+2.0*mm       -2.0*fm   )
       +    2.0    *sin(     mm    -ms+2.0*fm   )
       -    2.0    *sin(+3.0*mm    +ms          );

   jdk+=  +32.5    *sin( 299.77*DEG2RAD +  0.107408*DEG2RAD*k
                                        -  0.009173*DEG2RAD*T*T )
          +16.5    *sin( 251.88*DEG2RAD +  0.016321*DEG2RAD*k )
          +16.4    *sin( 251.83*DEG2RAD + 26.651886*DEG2RAD*k )
          +12.6    *sin( 349.42*DEG2RAD + 36.412478*DEG2RAD*k )
          +11.0    *sin(  84.66*DEG2RAD + 18.206239*DEG2RAD*k )
           +6.2    *sin( 141.74*DEG2RAD + 53.303771*DEG2RAD*k )
           +6.0    *sin( 207.14*DEG2RAD +  2.453732*DEG2RAD*k )
           +5.6    *sin( 154.84*DEG2RAD +  7.306860*DEG2RAD*k )
           +4.7    *sin(  34.52*DEG2RAD + 27.261239*DEG2RAD*k )
           +4.2    *sin( 207.19*DEG2RAD +  0.121824*DEG2RAD*k )
           +4.0    *sin( 291.34*DEG2RAD +  1.844379*DEG2RAD*k )
           +3.7    *sin( 161.72*DEG2RAD + 24.198154*DEG2RAD*k )
           +3.5    *sin( 239.56*DEG2RAD + 25.513099*DEG2RAD*k )
           +2.3    *sin( 331.55*DEG2RAD +  3.592518*DEG2RAD*k );


   w=0.00306-0.00038*E*cos(ms)+0.00026*cos(mm)
	-0.00002*cos(mm-ms)+0.00002*cos(mm+ms)+0.00002*cos(2.0*fm);

   return jd+jdk/100000.0+w;
}

double mond_voll( double jd_start )
/* Vollmond berechnen.
   Eingabe:
   jd_start: JD NACH dem der Vollmond folgen soll
   Ausgabe:
   JDE des Vollmondzeitpunktes
*/
{  double k,T,mm,fm,ms,km,jd,jdk,E;

   k=ceil((jd_start-NEU_MON)/SYN_MON-0.5)+0.5;
   T=k/1236.85;
   jd=mond_mittlere_phase(k,T);

   ms=  2.5534*DEG2RAD+ 29.10535669*DEG2RAD*k
		      -  0.0000218 *DEG2RAD*T*T
		      -  0.00000011*DEG2RAD*T*T*T;
   mm=201.5643*DEG2RAD+385.81693528*DEG2RAD*k
		      +  0.0107438 *DEG2RAD*T*T
		      +  0.00001239*DEG2RAD*T*T*T
		      -  0.000000058*DEG2RAD*T*T*T*T;
   fm=160.7108*DEG2RAD+390.67050274*DEG2RAD*k
		      -  0.0016341 *DEG2RAD*T*T
		      -  0.00000227*DEG2RAD*T*T*T
		      +  0.000000011*DEG2RAD*T*T*T*T;
   km=124.7746*DEG2RAD-  1.56375580*DEG2RAD*k
		      +  0.0020691 *DEG2RAD*T*T
		      +  0.00000215*DEG2RAD*T*T*T;
   E=1.0-0.002516*T-0.0000074*T*T;

   /*                        mm     ms     fm km   */
   jdk=-40614.0    *sin(     mm                 )
       +17302.0*E  *sin(            ms          )
       + 1614.0    *sin(+2.0*mm                 )
       + 1043.0    *sin(              +2.0*fm   )
       +  734.0*E  *sin(     mm    -ms          )
       -  515.0*E  *sin(     mm    +ms          )
       +  209.0*E*E*sin(       +2.0*ms          )
       -  111.0    *sin(     mm       -2.0*fm   )
       -   57.0    *sin(     mm       +2.0*fm   )
       +   56.0*E  *sin(+2.0*mm    +ms          )
       -   42.0    *sin(+3.0*mm                 )
       +   42.0*E  *sin(            ms+2.0*fm   )
       +   38.0*E  *sin(            ms-2.0*fm   )
       -   24.0*E  *sin(+2.0*mm    -ms          )
       -   17.0    *sin(                     +km)
       -    7.0    *sin(     mm+2.0*ms          )
       +    4.0    *sin(+2.0*mm       -2.0*fm   )
       +    4.0    *sin(       +3.0*ms          )
       +    3.0    *sin(     mm    +ms-2.0*fm   )
       +    3.0    *sin(+2.0*mm       +2.0*fm   )
       -    3.0    *sin(     mm    +ms+2.0*fm   )
       +    3.0    *sin(     mm    -ms+2.0*fm   )
       -    2.0    *sin(     mm    -ms-2.0*fm   )
       -    2.0    *sin(+3.0*mm    +ms          )
       +    2.0    *sin(+4.0*mm                 );

   jdk+=  +32.5    *sin( 299.77*DEG2RAD +  0.107408*DEG2RAD*k
                                        -  0.009173*DEG2RAD*T*T )
          +16.5    *sin( 251.88*DEG2RAD +  0.016321*DEG2RAD*k )
          +16.4    *sin( 251.83*DEG2RAD + 26.651886*DEG2RAD*k )
          +12.6    *sin( 349.42*DEG2RAD + 36.412478*DEG2RAD*k )
          +11.0    *sin(  84.66*DEG2RAD + 18.206239*DEG2RAD*k )
           +6.2    *sin( 141.74*DEG2RAD + 53.303771*DEG2RAD*k )
           +6.0    *sin( 207.14*DEG2RAD +  2.453732*DEG2RAD*k )
           +5.6    *sin( 154.84*DEG2RAD +  7.306860*DEG2RAD*k )
           +4.7    *sin(  34.52*DEG2RAD + 27.261239*DEG2RAD*k )
           +4.2    *sin( 207.19*DEG2RAD +  0.121824*DEG2RAD*k )
           +4.0    *sin( 291.34*DEG2RAD +  1.844379*DEG2RAD*k )
           +3.7    *sin( 161.72*DEG2RAD + 24.198154*DEG2RAD*k )
           +3.5    *sin( 239.56*DEG2RAD + 25.513099*DEG2RAD*k )
           +2.3    *sin( 331.55*DEG2RAD +  3.592518*DEG2RAD*k );


   return jd+jdk/100000.0;
}

double mond_letztes_viertel( double jd_start )
/* Letztes Viertel berechnen.
   Eingabe:
   jd_start: JD NACH dem das 3. 1/4 folgen soll
   Ausgabe:
   JDE des 3. 1/4
*/
{  double k,T,mm,fm,ms,km,jd,w,jdk,E;

   k=ceil((jd_start-NEU_MON)/SYN_MON-0.75)+0.75;
   T=k/1236.85;
   jd=mond_mittlere_phase(k,T);

   ms=  2.5534*DEG2RAD+ 29.10535669*DEG2RAD*k
		      -  0.0000218 *DEG2RAD*T*T
		      -  0.00000011*DEG2RAD*T*T*T;
   mm=201.5643*DEG2RAD+385.81693528*DEG2RAD*k
		      +  0.0107438 *DEG2RAD*T*T
		      +  0.00001239*DEG2RAD*T*T*T
		      -  0.000000058*DEG2RAD*T*T*T*T;
   fm=160.7108*DEG2RAD+390.67050274*DEG2RAD*k
		      -  0.0016341 *DEG2RAD*T*T
		      -  0.00000227*DEG2RAD*T*T*T
		      +  0.000000011*DEG2RAD*T*T*T*T;
   km=124.7746*DEG2RAD-  1.56375580*DEG2RAD*k
		      +  0.0020691 *DEG2RAD*T*T
		      +  0.00000215*DEG2RAD*T*T*T;
   E=1.0-0.002516*T-0.0000074*T*T;

   /*                        mm     ms     fm km   */
   jdk=-62801.0    *sin(     mm                 )
       +17172.0*E  *sin(            ms          )
       - 1183.0*E  *sin(     mm    +ms          )
       +  862.0    *sin(+2.0*mm                 )
       +  804.0    *sin(              +2.0*fm   )
       +  454.0*E  *sin(     mm    -ms          )
       +  204.0*E*E*sin(       +2.0*ms          )
       -  180.0    *sin(     mm       -2.0*fm   )
       -   70.0    *sin(     mm       +2.0*fm   )
       -   40.0    *sin(+3.0*mm                 )
       -   34.0*E  *sin(+2.0*mm    -ms          )
       +   32.0*E  *sin(            ms+2.0*fm   )
       +   32.0*E  *sin(            ms-2.0*fm   )
       -   28.0*E*E*sin(     mm+2.0*ms          )
       +   27.0*E  *sin(+2.0*mm    +ms          )
       -   17.0    *sin(                     +km)
       -    5.0    *sin(     mm    -ms-2.0*fm   )
       +    4.0    *sin(+2.0*mm       +2.0*fm   )
       -    4.0    *sin(     mm    +ms+2.0*fm   )
       +    4.0    *sin(     mm-2.0*ms          )
       +    3.0    *sin(     mm    +ms-2.0*fm   )
       +    3.0    *sin(       +3.0*ms          )
       +    2.0    *sin(+2.0*mm       -2.0*fm   )
       +    2.0    *sin(     mm    -ms+2.0*fm   )
       -    2.0    *sin(+3.0*mm    +ms          );

   jdk+=  +32.5    *sin( 299.77*DEG2RAD +  0.107408*DEG2RAD*k
                                        -  0.009173*DEG2RAD*T*T )
          +16.5    *sin( 251.88*DEG2RAD +  0.016321*DEG2RAD*k )
          +16.4    *sin( 251.83*DEG2RAD + 26.651886*DEG2RAD*k )
          +12.6    *sin( 349.42*DEG2RAD + 36.412478*DEG2RAD*k )
          +11.0    *sin(  84.66*DEG2RAD + 18.206239*DEG2RAD*k )
           +6.2    *sin( 141.74*DEG2RAD + 53.303771*DEG2RAD*k )
           +6.0    *sin( 207.14*DEG2RAD +  2.453732*DEG2RAD*k )
           +5.6    *sin( 154.84*DEG2RAD +  7.306860*DEG2RAD*k )
           +4.7    *sin(  34.52*DEG2RAD + 27.261239*DEG2RAD*k )
           +4.2    *sin( 207.19*DEG2RAD +  0.121824*DEG2RAD*k )
           +4.0    *sin( 291.34*DEG2RAD +  1.844379*DEG2RAD*k )
           +3.7    *sin( 161.72*DEG2RAD + 24.198154*DEG2RAD*k )
           +3.5    *sin( 239.56*DEG2RAD + 25.513099*DEG2RAD*k )
           +2.3    *sin( 331.55*DEG2RAD +  3.592518*DEG2RAD*k );

   w=0.00306-0.00038*E*cos(ms)+0.00026*cos(mm)
	-0.00002*cos(mm-ms)+0.00002*cos(mm+ms)+0.00002*cos(2.0*fm);

   return jd=jd+jdk/100000.0-w;
}

double mond_alter( double jd )
/* Mondalter (Tage nach Neumond) berechnen
   Eingabe:
   jd : JD der Beobachtung
   Return:
   Mondalter in Tagen
*/
{  double jdt;

   jdt=mond_neu(jd-SYN_MON+1.5);
   if ((jd-jdt)<0) jdt=mond_neu(jd-SYN_MON-1.5);

   return jd-jdt;
}

void mond_apsiden_phasen( double jd_start, double jde[6], short art[6] )
/* Berechnung von den naechsten vier Mondphasen nach jd_start
   Eingabe:
   jd_start : Datum, nach dem die Phasen liegen sollen
   Ausgabe: (nach aufsteigendem jde sortiert)
   jde      : JED der Phasen
   art      : Art des Ereignisses (0=0/4=Neumond,
				   1=1/4=Erstes Viertel,
				   2=2/4=Vollmond,
				   3=3/4=Letztes Viertel
				   4    =Perigaeum
				   5    =Apogaeum)
*/
{  short a,b;

   art[0]=0;            jde[0]=mond_neu(jd_start-1.0);
   if (jde[0]<jd_start) jde[0]=mond_neu(jd_start+1.0);

   art[1]=1;            jde[1]=mond_erstes_viertel(jd_start-1.0);
   if (jde[1]<jd_start) jde[1]=mond_erstes_viertel(jd_start+1.0);

   art[2]=2;            jde[2]=mond_voll(jd_start-1.0);
   if (jde[2]<jd_start) jde[2]=mond_voll(jd_start+1.0);

   art[3]=3;            jde[3]=mond_letztes_viertel(jd_start-1.0);
   if (jde[3]<jd_start) jde[3]=mond_letztes_viertel(jd_start+1.0);

   art[4]=4;            jde[4]=mond_perigaeum(jd_start-2.0);
   if (jde[4]<jd_start) jde[4]=mond_perigaeum(jd_start+2.0);

   art[5]=5;            jde[5]=mond_apogaeum(jd_start-2.0);
   if (jde[5]<jd_start) jde[5]=mond_apogaeum(jd_start+2.0);

   /* Insertion Sort */
   for(a=1;a<6;++a)
   {  double jde_old=jde[a];
      short    art_old=art[a];
      for(b=a-1;b>=0 && jde_old<jde[b];b--)
      {  jde[b+1]=jde[b];
	 art[b+1]=art[b];
      }
      jde[b+1]=jde_old;
      art[b+1]=art_old;
   }
}
