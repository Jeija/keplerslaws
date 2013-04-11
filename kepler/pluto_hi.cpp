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

#include <stdio.h>
#include <math.h>
#include "astromat.h"
#include "cartes.h"
#include "einaus.h"

sphaer pluto_hi( double jd )
/* Pluto-Reihenentwicklung nach

   An accurat representation of the motion of Pluto
   E. Goffin, J.Meeus und C.Steyaert
   Astronomy and Astrophysics 155,323-325 (1986)

   Genauigkeit (1885-2099):
   Laenge und Breite besser als 1"
   Entfernung besser als 0.0001 AE

   Eingabe:
   jd : Julianisches (Ephemeriden-) Datum

   Return:
   Heliozentrische ekliptikale Koordinaten
   fuer das Aequinoktium des Datums (nach AmdPC S.232)
*/
{  double t,j,s,p,l,b,r,T,Pia,pia,pa,c1,c2,c3,s1,s2,s3,x,y,z;
   sphaer plu;

   t=(jd-2415020.0)/36525.0; // Julianische Jahrhunderte seit 1900

   j=238.74*DEG2RAD+3034.9057*DEG2RAD*t;
   s=267.26*DEG2RAD+1222.1138*DEG2RAD*t;
   p= 93.48*DEG2RAD+ 144.9600*DEG2RAD*t;

/*                        j     s     p  */
l    =-19977972.0*sin(                p)
      +19667536.0*cos(                p)
	+987114.0*sin(            2.0*p)
       -4939350.0*cos(            2.0*p)
	+577978.0*sin(            3.0*p)
       +1226898.0*cos(            3.0*p)
	-334695.0*sin(            4.0*p)
	-201966.0*cos(            4.0*p)
	+130519.0*sin(            5.0*p)
	 -29025.0*cos(            5.0*p)
	 -39851.0*sin(            6.0*p)
	 +28968.0*cos(            6.0*p)
	 +20387.0*sin(          s    -p)
	  -9832.0*cos(          s    -p)
	  -3986.0*sin(          s      )
	  -4954.0*cos(          s      )
	  -5817.0*sin(          s    +p)
	  -3365.0*cos(          s    +p)
	  -3903.0*sin(          s+2.0*p)
	  +2895.0*cos(          s+2.0*p);
l+=	   -738.0*sin(          s+3.0*p)
	  +3443.0*cos(          s+3.0*p)
	  +1234.0*sin(      2.0*s-2.0*p)
	   +472.0*cos(      2.0*s-2.0*p)
	  +1101.0*sin(      2.0*s    -p)
	   -894.0*cos(      2.0*s    -p)
	   +625.0*sin(      2.0*s      )
	  -1214.0*cos(      2.0*s      )
	  +2485.0*sin(    j    -s      )
	   -486.0*cos(    j    -s      )
	   +852.0*sin(    j    -s    +p)
	  -1407.0*cos(    j    -s    +p)
	   -948.0*sin(    j      -3.0*p)
	  +1073.0*cos(    j      -3.0*p)
	  -2309.0*sin(    j      -2.0*p)
	  -1024.0*cos(    j      -2.0*p)
	  +7047.0*sin(    j 	     -p)
	   +770.0*cos(    j 	     -p)
	  +1184.0*sin(    j 	       )
	   -344.0*cos(    j 	       );
l+=	   +394.0*sin(    j 	     +p)
	    -55.0*cos(    j 	     +p)
	   +119.0*sin(    j      +2.0*p)
	   -264.0*cos(    j      +2.0*p)
	    -46.0*sin(    j      +3.0*p)
	   -156.0*cos(    j      +3.0*p)
	    -77.0*sin(    j      +4.0*p)
	    -33.0*cos(    j      +4.0*p)
	    -34.0*sin(    j    +s-3.0*p)
	    -26.0*cos(    j    +s-3.0*p)
	    -43.0*sin(    j    +s-2.0*p)
	    -15.0*sin(    j    +s    -p)
	    +21.0*cos(    j    +s    -p)
	    +10.0*sin(2.0*j      -3.0*p)
	    +22.0*cos(2.0*j      -3.0*p)
	    -57.0*sin(2.0*j      -2.0*p)
	    -32.0*cos(2.0*j      -2.0*p)
	   +158.0*sin(2.0*j 	     -p)
	    -43.0*cos(2.0*j 	     -p);

/*                        j     s     p  */
b    = -5323113.0*sin(                p)
      -15024245.0*cos(                p)
       +3497557.0*sin(            2.0*p)
       +1735457.0*cos(            2.0*p)
       -1059559.0*sin(            3.0*p)
	+299464.0*cos(            3.0*p)
	+189102.0*sin(            4.0*p)
	-285383.0*cos(            4.0*p)
	 +14231.0*sin(            5.0*p)
	+101218.0*cos(            5.0*p)
	 -29164.0*sin(            6.0*p)
	 -27461.0*cos(            6.0*p)
	  +4935.0*sin(          s    -p)
	 +11282.0*cos(          s    -p)
	   +312.0*sin(          s      )
	   -128.0*cos(          s      )
	  +2057.0*sin(          s    +p)
	   -904.0*cos(          s    +p)
	    +19.0*sin(          s+2.0*p)
	   -674.0*cos(          s+2.0*p);
b+=	   -307.0*sin(          s+3.0*p)
	   -576.0*cos(          s+3.0*p)
	    -65.0*sin(      2.0*s-2.0*p)
	    +39.0*cos(      2.0*s-2.0*p)
	    -97.0*sin(      2.0*s    -p)
	   +208.0*cos(      2.0*s    -p)
	   -160.0*cos(      2.0*s      )
	   -177.0*sin(    j    -s      )
	   +259.0*cos(    j    -s      )
	    +15.0*sin(    j    -s    +p)
	   +235.0*cos(    j    -s    +p)
	   +578.0*sin(    j      -3.0*p)
	   -293.0*cos(    j      -3.0*p)
	   -294.0*sin(    j      -2.0*p)
	   +694.0*cos(    j      -2.0*p)
	   +156.0*sin(    j 	     -p)
	   +201.0*cos(    j 	     -p)
	   +294.0*sin(    j 	       )
	   +829.0*cos(    j 	       );
b+=	   -123.0*sin(    j 	     +p)
	    -31.0*cos(    j 	     +p);

/*                        j     s     p  */
r    =  6623876.0*sin(                p)
       +6955990.0*cos(                p)
       -1181808.0*sin(            2.0*p)
	 -54836.0*cos(            2.0*p)
	+163227.0*sin(            3.0*p)
	-139603.0*cos(            3.0*p)
	  -3644.0*sin(            4.0*p)
	 +48144.0*cos(            4.0*p)
	  -6268.0*sin(            5.0*p)
	  -8851.0*cos(            5.0*p)
	  +3111.0*sin(            6.0*p)
	   -408.0*cos(            6.0*p)
	   -621.0*sin(          s    -p)
	  +2223.0*cos(          s    -p)
	   +438.0*sin(          s      )
	   +450.0*cos(          s      )
	   -153.0*sin(          s    +p)
	    +61.0*cos(          s    +p)
	     -3.0*sin(          s+2.0*p)
	    +79.0*cos(          s+2.0*p);
r+=	    +50.0*sin(          s+3.0*p)
	    +54.0*cos(          s+3.0*p)
		 -sin(      2.0*s-2.0*p)
	    -22.0*cos(      2.0*s-2.0*p)
	    +84.0*sin(      2.0*s    -p)
	    -48.0*cos(      2.0*s    -p)
	    -30.0*sin(      2.0*s      )
	    +61.0*cos(      2.0*s      )
	    +26.0*sin(    j    -s      )
	    -39.0*cos(    j    -s      )
	    -19.0*sin(    j    -s    +p)
	    -40.0*cos(    j    -s    +p)
	   -321.0*sin(    j      -3.0*p)
	    +42.0*cos(    j      -3.0*p)
	   +797.0*sin(    j      -2.0*p)
	   -792.0*cos(    j      -2.0*p)
	     -4.0*sin(    j 	     -p)
	  +4564.0*cos(    j 	     -p)
	   +852.0*sin(    j 	       )
	   +855.0*cos(    j 	       );
r+=	    -88.0*sin(    j 	     +p)
	    -82.0*cos(    j 	     +p)
	    +21.0*sin(    j      +2.0*p)
	    -12.0*cos(    j      +2.0*p)
	    -14.0*sin(    j      +3.0*p)
	     +6.0*cos(    j      +3.0*p)
	     -6.0*sin(2.0*j      -3.0*p)
		 +cos(2.0*j      -3.0*p)
	    +13.0*sin(2.0*j      -2.0*p)
	    -23.0*cos(2.0*j      -2.0*p)
	    +25.0*sin(2.0*j 	     -p)
	   +107.0*cos(2.0*j 	     -p)
	    +25.0*sin(2.0*j 	       )
	    +16.0*cos(2.0*j 	       );

   plu.lambda=93.297471*DEG2RAD+144.9600*DEG2RAD*t+0.000001*DEG2RAD*l;
   plu.beta  =-3.909434*DEG2RAD                   +0.000001*DEG2RAD*b;
   plu.r     =40.724725                           +0.000001*r;

   T=t-0.5;  // Julianische Jahrhunderte seit B1950

   Pia= 3.044;
   pia= 0.000228*T;
   pa =(0.0243764+5.39E-6*T)*T;

   c1=cos(pia); c2=cos(plu.beta); c3=cos(Pia-plu.lambda);
   s1=sin(pia); s2=sin(plu.beta); s3=sin(Pia-plu.lambda);

   x=c2*c3;
   y=c1*c2*s3-s1*s2;
   z=s1*c2*s3+c1*s2;

   plu.beta=atan(z/sqrt((1.0-z)*(1.0+z)));
   if (x>0) plu.lambda=M_2PI*frac((Pia+pa-atan(y/x))/M_2PI);
   else     plu.lambda=M_2PI*frac((Pia+pa-atan(y/x))/M_2PI+0.5);

   return plu;
}
