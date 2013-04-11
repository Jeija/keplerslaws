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

/* BUG behoben: bei l,D,M,m,F war fuer T^3 und T^4 DEG2RAD falsch
   geklammert TJK 4.12.1994 */

#include <math.h>
#include "astromat.h"

sphaer mond_genau( double jde )
{
   double T,l,D,M,m,F,E,E2,A1,A2,A3,delta_l,delta_b,delta_r;

   sphaer s;

   T=(jde-J2000)/36525.0;

   l=mod(+218.3164591*DEG2RAD
	 +481267.88134236*DEG2RAD*T
	 -0.0013268*DEG2RAD*T*T
	 +T*T*T/538841.0*DEG2RAD
	 -T*T*T*T/65194000.0*DEG2RAD,M_2PI);

   D=mod(+297.8502042*DEG2RAD
	 +445267.1115168*DEG2RAD*T
	 -0.0016300*DEG2RAD*T*T
	 +T*T*T/545868.0*DEG2RAD
	 -T*T*T*T/113065000.0*DEG2RAD,M_2PI);

   M=mod(+357.5291092*DEG2RAD
	 +35999.0502909 *DEG2RAD*T
	 -0.0001536*DEG2RAD*T*T
	 +T*T*T/24490000.0*DEG2RAD,M_2PI);

   m=mod(+134.9634114*DEG2RAD
	 +477198.8676313*DEG2RAD*T
	 +0.0089970*DEG2RAD*T*T
	 +T*T*T/69699.0*DEG2RAD
	 -T*T*T*T/14712000.0*DEG2RAD,M_2PI);

   F=mod(+93.2720993*DEG2RAD
	 +483202.0175273*DEG2RAD*T
	 -0.0034029*DEG2RAD*T*T
	 -T*T*T/3526000.0*DEG2RAD
	 +T*T*T*T/863310000.0*DEG2RAD,M_2PI);

   E=1.0-0.002516*T-0.0000074*T*T;
   E2=E*E;

   A1=mod(119.75*DEG2RAD +    131.849*DEG2RAD*T,M_2PI);
   A2=mod( 53.09*DEG2RAD + 479264.290*DEG2RAD*T,M_2PI);
   A3=mod(313.45*DEG2RAD + 481266.484*DEG2RAD*T,M_2PI);

//   T = -0.077221081451;
//   l =134.290186*DEG2RAD;
//   D =113.842309*DEG2RAD;
//   M = 97.643514*DEG2RAD;
//   m =  5.150839*DEG2RAD;
//   F =219.889726*DEG2RAD;
//   A1=109.57*DEG2RAD;
//   A2=123.78*DEG2RAD;
//   A3=229.53*DEG2RAD;
//   E =  1.000194;
//   E2=E*E;

   /* Laenge */
   delta_l  =
	     +sin(                +m      ) * ( +6288774.0)
	     +sin(+2.0*D          -m      ) * ( +1274027.0)
	     +sin(+2.0*D                  ) * (  +658314.0)
             +sin(            +2.0*m      ) * (  +213618.0)
             +sin(          +M            ) * (  -185116.0) * E
             +sin(                  +2.0*F) * (  -114332.0)
             +sin(+2.0*D      -2.0*m      ) * (   +58793.0)
             +sin(+2.0*D    -M    -m      ) * (   +57066.0) * E
	     +sin(+2.0*D          +m      ) * (   +53322.0)
	     +sin(+2.0*D    -M            ) * (   +45758.0) * E
             +sin(          +M    -m      ) * (   -40923.0) * E
             +sin(    +D                  ) * (   -34720.0)
             +sin(          +M    +m      ) * (   -30383.0) * E
             +sin(+2.0*D            -2.0*F) * (   +15327.0)
             +sin(                +m+2.0*F) * (   -12528.0)
             +sin(                +m-2.0*F) * (   +10980.0)
             +sin(+4.0*D          -m      ) * (   +10675.0)
             +sin(            +3.0*m      ) * (   +10034.0)
	     +sin(+4.0*D      -2.0*m      ) * (    +8548.0)
             +sin(+2.0*D    +M    -m      ) * (    -7888.0) * E;
   delta_l +=
	     +sin(+2.0*D    +M            ) * (    -6766.0) * E
	     +sin(    +D          -m      ) * (    -5163.0)
             +sin(    +D    +M            ) * (    +4987.0) * E
             +sin(+2.0*D    -M    +m      ) * (    +4036.0) * E
             +sin(+2.0*D      +2.0*m      ) * (    +3994.0)
             +sin(+4.0*D                  ) * (    +3861.0)
             +sin(+2.0*D      -3.0*m      ) * (    +3665.0)
	     +sin(          +M-2.0*m      ) * (    -2689.0) * E
	     +sin(+2.0*D          -m+2.0*F) * (    -2602.0)
             +sin(+2.0*D    -M-2.0*m      ) * (    +2390.0) * E
             +sin(    +D          +m      ) * (    -2348.0)
             +sin(+2.0*D-2.0*M            ) * (    +2236.0) * E2
             +sin(          +M+2.0*m      ) * (    -2120.0) * E
             +sin(      +2.0*M            ) * (    -2069.0) * E2
             +sin(+2.0*D-2.0*M    -m      ) * (    +2048.0) * E2
             +sin(+2.0*D          +m-2.0*F) * (    -1773.0)
             +sin(+2.0*D            +2.0*F) * (    -1595.0)
	     +sin(+4.0*D    -M    -m      ) * (    +1215.0) * E
             +sin(            +2.0*m+2.0*F) * (    -1110.0)
	     +sin(+3.0*D          -m      ) * (     -892.0);
   delta_l +=
	     +sin(+2.0*D    +M    +m      ) * (     -810.0) * E
             +sin(+4.0*D    -M-2.0*m      ) * (     +759.0) * E
             +sin(      +2.0*M    -m      ) * (     -713.0) * E2
             +sin(+2.0*D+2.0*M    -m      ) * (     -700.0) * E2
             +sin(+2.0*D    +M-2.0*m      ) * (     +691.0) * E
             +sin(+2.0*D    -M      -2.0*F) * (     +596.0) * E
	     +sin(+4.0*D          +m      ) * (     +549.0)
	     +sin(            +4.0*m      ) * (     +537.0)
             +sin(+4.0*D    -M            ) * (     +520.0) * E
             +sin(    +D      -2.0*m      ) * (     -487.0)
             +sin(+2.0*D    +M      -2.0*F) * (     -399.0) * E
             +sin(            +2.0*m-2.0*F) * (     -381.0)
             +sin(    +D    +M    +m      ) * (     +351.0) * E
             +sin(+3.0*D      -2.0*m      ) * (     -340.0)
             +sin(+4.0*D      -3.0*m      ) * (     +330.0)
             +sin(+2.0*D    -M+2.0*m      ) * (     +327.0) * E
	     +sin(      +2.0*M    +m      ) * (     -323.0) * E2
             +sin(    +D    +M    -m      ) * (     +299.0) * E
	     +sin(+2.0*D      +3.0*m      ) * (     +294.0);

   /* Breite */
   delta_b  =
	     +sin(                      +F) * ( +5128122.0)
	     +sin(                +m    +F) * (  +280602.0)
	     +sin(                +m    -F) * (  +277693.0)
	     +sin(+2.0*D                -F) * (  +173237.0)
	     +sin(+2.0*D          -m    +F) * (   +55413.0)
	     +sin(+2.0*D          -m    -F) * (   +46271.0)
             +sin(+2.0*D                +F) * (   +32573.0)
             +sin(            +2.0*m    +F) * (   +17198.0)
             +sin(+2.0*D          +m    -F) * (    +9266.0)
             +sin(            +2.0*m    -F) * (    +8822.0)
             +sin(+2.0*D    -M          -F) * (    +8216.0) * E
             +sin(+2.0*D      -2.0*m    -F) * (    +4324.0)
             +sin(+2.0*D          +m    +F) * (    +4200.0)
             +sin(+2.0*D    +M          -F) * (    -3359.0) * E
	     +sin(+2.0*D    -M    -m    +F) * (    +2463.0) * E
	     +sin(+2.0*D    -M          +F) * (    +2211.0) * E
	     +sin(+2.0*D    -M    -m    -F) * (    +2065.0) * E
	     +sin(          +M    -m    -F) * (    -1870.0) * E
	     +sin(+4.0*D          -m    -F) * (    +1828.0)
	     +sin(          +M          +F) * (    -1794.0) * E;
   delta_b +=
	     +sin(                  +3.0*F) * (    -1749.0)
	     +sin(          +M    -m    +F) * (    -1565.0) * E
	     +sin(    +D                +F) * (    -1491.0)
	     +sin(          +M    +m    +F) * (    -1475.0) * E
	     +sin(          +M    +m    -F) * (    -1410.0) * E
	     +sin(          +M          -F) * (    -1344.0) * E
	     +sin(    +D                -F) * (    -1335.0)
	     +sin(            +3.0*m    +F) * (    +1107.0)
	     +sin(+4.0*D                -F) * (    +1021.0)
	     +sin(+4.0*D          -m    +F) * (     +833.0)
	     +sin(                +m-3.0*F) * (     +777.0)
	     +sin(+4.0*D      -2.0*m    +F) * (     +671.0)
	     +sin(+2.0*D            -3.0*F) * (     +607.0)
	     +sin(+2.0*D      +2.0*m    -F) * (     +596.0)
	     +sin(+2.0*D    -M    +m    -F) * (     +491.0) * E
	     +sin(+2.0*D      -2.0*m    +F) * (     -451.0)
	     +sin(            +3.0*m    -F) * (     +439.0)
	     +sin(+2.0*D      +2.0*m    +F) * (     +422.0)
	     +sin(+2.0*D      -3.0*m    -F) * (     +421.0)
	     +sin(+2.0*D    +M    -m    +F) * (     -366.0) * E;
   delta_b +=
	     +sin(+2.0*D    +M          +F) * (     -351.0) * E
	     +sin(+4.0*D                +F) * (     +331.0)
	     +sin(+2.0*D    -M    +m    +F) * (     +315.0) * E
	     +sin(+2.0*D-2.0*M          -F) * (     +302.0) * E2
	     +sin(                +m+3.0*F) * (     -283.0)
	     +sin(+2.0*D    +M    +m    -F) * (     -229.0) * E
	     +sin(    +D    +M          -F) * (     +223.0) * E
	     +sin(    +D    +M          +F) * (     +223.0) * E
	     +sin(          +M-2.0*m    -F) * (     -220.0) * E
	     +sin(+2.0*D    +M    -m    -F) * (     -220.0) * E
	     +sin(    +D          +m    +F) * (     -185.0)
	     +sin(+2.0*D    -M-2.0*m    -F) * (     +181.0) * E
	     +sin(          +M+2.0*m    +F) * (     -177.0) * E
	     +sin(+4.0*D      -2.0*m    -F) * (     +176.0)
	     +sin(+4.0*D    -M    -m    -F) * (     +166.0) * E
	     +sin(    +D          +m    -F) * (     -164.0)
	     +sin(+4.0*D          +m    -F) * (     +132.0)
	     +sin(    +D          -m    -F) * (     -119.0)
             +sin(+4.0*D    -M          -F) * (     +115.0) * E
	     +sin(+2.0*D-2.0*M          +F) * (     +107.0) * E2;

   /* Radiusvektor */
   delta_r  =
	     +cos(                +m      ) * (-20905355.0)
	     +cos(+2.0*D          -m      ) * ( -3699111.0)
             +cos(+2.0*D                  ) * ( -2955968.0)
             +cos(            +2.0*m      ) * (  -569925.0)
             +cos(          +M            ) * (   +48888.0) * E
             +cos(                  +2.0*F) * (    -3149.0)
             +cos(+2.0*D      -2.0*m      ) * (  +246158.0)
             +cos(+2.0*D    -M    -m      ) * (  -152138.0) * E
             +cos(+2.0*D          +m      ) * (  -170733.0)
             +cos(+2.0*D    -M            ) * (  -204586.0) * E
             +cos(          +M    -m      ) * (  -129620.0) * E
             +cos(    +D                  ) * (  +108743.0)
             +cos(          +M    +m      ) * (  +104755.0) * E
             +cos(+2.0*D            -2.0*F) * (   +10321.0)
             +cos(                +m-2.0*F) * (   +79661.0)
             +cos(+4.0*D          -m      ) * (   -34782.0)
             +cos(            +3.0*m      ) * (   -23210.0)
             +cos(+4.0*D      -2.0*m      ) * (   -21636.0)
             +cos(+2.0*D    +M    -m      ) * (   +24208.0) * E;
   delta_r +=
	     +cos(+2.0*D    +M            ) * (   +30824.0) * E
	     +cos(    +D          -m      ) * (    -8379.0)
             +cos(    +D    +M            ) * (   -16675.0) * E
             +cos(+2.0*D    -M    +m      ) * (   -12831.0) * E
             +cos(+2.0*D      +2.0*m      ) * (   -10445.0)
             +cos(+4.0*D                  ) * (   -11650.0)
             +cos(+2.0*D      -3.0*m      ) * (   +14403.0)
             +cos(          +M-2.0*m      ) * (    -7003.0) * E
             +cos(+2.0*D    -M-2.0*m      ) * (   +10056.0) * E
             +cos(    +D          +m      ) * (    +6322.0)
             +cos(+2.0*D-2.0*M            ) * (    -9884.0) * E2
             +cos(          +M+2.0*m      ) * (    +5751.0) * E
             +cos(+2.0*D-2.0*M    -m      ) * (    -4950.0) * E2
             +cos(+2.0*D          +m-2.0*F) * (    +4130.0)
             +cos(+4.0*D    -M    -m      ) * (    -3958.0) * E
             +cos(+3.0*D          -m      ) * (    +3258.0);
   delta_r +=
	     +cos(+2.0*D    +M    +m      ) * (    +2616.0) * E
             +cos(+4.0*D    -M-2.0*m      ) * (    -1897.0) * E
             +cos(      +2.0*M    -m      ) * (    -2117.0) * E2
             +cos(+2.0*D+2.0*M    -m      ) * (    +2354.0) * E2
	     +cos(+4.0*D          +m      ) * (    -1423.0)
             +cos(            +4.0*m      ) * (    -1117.0)
             +cos(+4.0*D    -M            ) * (    -1571.0) * E
             +cos(    +D      -2.0*m      ) * (    -1739.0)
             +cos(            +2.0*m-2.0*F) * (    -4421.0)
             +cos(      +2.0*M    +m      ) * (    +1165.0) * E2
	     +cos(+2.0*D          -m-2.0*F) * (    +8752.0);

   delta_l += + 3958.0 * sin (A1)
	      + 1962.0 * sin (l - F)
	      +  318.0 * sin (A2);

   delta_b += - 2235.0 * sin (l)
	      +  382.0 * sin (A3)
	      +  175.0 * sin (A1 - F)
	      +  175.0 * sin (A1 + F)
	      +  127.0 * sin (l - m)
	      -  115.0 * sin (l + m);

//   printf("\n%.1f %.1f %.1f",delta_l,delta_b,delta_r/1000.0);

   s.lambda=mod(l    +delta_l*0.000001*DEG2RAD,M_2PI);
   s.beta  =          delta_b*0.000001*DEG2RAD;
   s.r     =385000.56+delta_r*0.001;

   return s;
}

/*

void test( void )
{
   sphaer s;
   DMS l,b;

   s=mond_hi(juldat(1992,4,12,0.0));

   l=s.lambda*RAD2DEG; l.format(SS,1);
   b=s.beta  *RAD2DEG; b.format(SS,1);

   printf("\n%12.6f %c %3d %2d %4.1f",s.lambda*RAD2DEG,(l.v<0)?'-':'+',l.i_d,l.i_m,l.d_s);
   printf("\n%12.6f %c %3d %2d %4.1f",s.beta  *RAD2DEG,(b.v<0)?'-':'+',b.i_d,b.i_m,b.d_s);
   printf("\n%10.2f",s.r);

   return 0;
}

*/
