/*
 * vsop87.cc, originally distributed along with xephem 3.2.2,
 * MODIFIED BY Tobias Kramer 1999, NO LONGER COMPATIBLE WITH XEPHEM
 */

/* VSOP87 planetary theory
 *
 * currently uses version VSOP87D:
 * heliocentric spherical, mean ecliptic of date.
 *
 * calculation of rates (daily changes) is optional;
 * see header file for the necessary #define's
 *
 * rough orientation on calculation time, miliseconds
 * on an HP 715/75, all planets Mercury to Neptune, prec=0.0:
 *
 *      terms	with rates	without rates
 * 	3598	11		7.1
 *      31577	51		44
 *
 * with secular terms for JD 2232395.0  19/12/1399 0h TDB:
 *
 *	FULL PRECISION code (31577 terms), milliseconds
 *	prec	terms	rates	no rates
 *	1e-8	15086	62	36
 *	1e-7	10105	44	25
 *	1e-6	3725	20	13
 *	1e-5	1324	11	7.8
 *	1e-4	443	7.0	6.0
 *	1e-3	139	6.0	5.0
 *
 *	REDUCED PRECISION code (3598 terms), milliseconds
 *	prec	terms	rates	no rates
 *	1e-7	2463	9.9	5.5
 *	1e-6	1939	8.0	4.5
 *	1e-5	1131	4.9	2.9
 *	1e-4	443	2.2	1.5
 *	1e-3	139	1.0	0.9
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "astromat.h"
#include "vsop87.h"

#define MJD_J2000 (2451545.0 - 2415020.0)


#define VSOP_A1000	365250.0	/* days per millenium */
#define VSOP_MAXALPHA	5		/* max degree of time */

/******************************************************************
 * adapted from BdL FORTRAN Code; stern, modified: Tobias Kramer
 *
 *    Reference : Bureau des Longitudes - PBGF9502
 *
 *    Object :  calculate a VSOP87 position for a given time.
 *
 *    Input :
 *
 *    mjd      modified julian date, counted from J1900.0
 *             time scale : dynamical time TDB.
 *
 *    obj	object number, NB: not for pluto
 *                Mercury    1
 *                Venus      2
 *                Earth      3
 *                Mars       4
 *                Jupiter    5
 *                Saturn     6
 *                Uranus     7
 *                Neptune    8
 *
 *
 *    prec     relative precision
 *
 *             if prec is equal to 0 then the precision is the precision
 *                p0 of the complete solution VSOP87.
 *                Mercury    p0 =  0.6 10**-8
 *                Venus      p0 =  2.5 10**-8
 *                Earth      p0 =  2.5 10**-8
 *                Mars       p0 = 10.0 10**-8
 *                Jupiter    p0 = 35.0 10**-8
 *                Saturn     p0 = 70.0 10**-8
 *                Uranus     p0 =  8.0 10**-8
 *                Neptune    p0 = 42.0 10**-8
 *
 *             if prec is not equal to 0, let us say in between p0 and
 *             10**-3, the precision is :
 *                for the positions :
 *                - prec*a0 au for the distances.
 *                - prec rad for the other variables.
 *                for the velocities :
 *                - prec*a0 au/day for the distances.
 *                - prec rad/day for the other variables.
 *                  a0 is the semi-major axis of the body.
 *
 *    Output :
 *
 *    ret[6]     array of the results (double).
 *
 *             for spherical coordinates :
 *                 1: longitude (rd)
 *                 2: latitude (rd)
 *                 3: radius (au)
 *		#if VSOP_GETRATE:
 *                 4: longitude velocity (rad/day)
 *                 5: latitude velocity (rad/day)
 *                 6: radius velocity (au/day)
 *
 *    return:     error index (int)
 *                 0: no error.
 *		   2: object out of range [MERCURY .. NEPTUNE, SUN]
 *		   3: precision out of range [0.0 .. 1e-3]
 ******************************************************************/
static int vn_mercury[][3] = {
	/* addresses for mercury l, b, r  */
	/* T^0 */ { 0, 144, 217, },
	/* T^1 */ { 75, 177, 246, },
	/* T^2 */ { 112, 192, 258, },
	/* T^3 */ { 126, 203, 266, },
	/* T^4 */ { 135, 211, 273, },
	/* T^5 */ { 143, 216, 277, },
	/* end */ { 144, 217, 0, },
	/* termination */ { 0, }
};
static int vn_venus[][3] = {
	/* addresses for venus l, b, r  */
	/* T^0 */ { 0, 91, 134, },
	/* T^1 */ { 47, 109, 156, },
	/* T^2 */ { 70, 121, 163, },
	/* T^3 */ { 81, 125, 166, },
	/* T^4 */ { 85, 129, 168, },
	/* T^5 */ { 88, 133, 169, },
	/* end */ { 91, 134, 0, },
	/* termination */ { 0, }
};
static int vn_earth[][3] = {
	/* addresses for earth l, b, r  */
	/* T^0 */ { 0, 247, 274, },
	/* T^1 */ { 117, 267, 345, },
	/* T^2 */ { 188, 272, 381, },
	/* T^3 */ { 232, 274, 398, },
	/* T^4 */ { 240, 0, 402, },
	/* T^5 */ { 244, 0, 405, },
	/* end */ { 247, 0, 406, },
	/* termination */ { 0, }
};
static int vn_mars[][3] = {
	/* addresses for mars l, b, r  */
	/* T^0 */ { 0, 343, 433, },
	/* T^1 */ { 124, 383, 528, },
	/* T^2 */ { 222, 409, 612, },
	/* T^3 */ { 289, 419, 665, },
	/* T^4 */ { 321, 427, 679, },
	/* T^5 */ { 336, 432, 688, },
	/* end */ { 343, 433, 691, },
	/* termination */ { 0, }
};
static int vn_jupiter[][3] = {
	/* addresses for jupiter l, b, r  */
	/* T^0 */ { 0, 392, 542, },
	/* T^1 */ { 109, 428, 672, },
	/* T^2 */ { 194, 468, 785, },
	/* T^3 */ { 271, 504, 903, },
	/* T^4 */ { 341, 526, 1001, },
	/* T^5 */ { 383, 538, 1047, },
	/* end */ { 392, 542, 1056, },
	/* termination */ { 0, }
};
static int vn_saturn[][3] = {
	/* addresses for saturn l, b, r  */
	/* T^0 */ { 0, 592, 811, },
	/* T^1 */ { 165, 645, 1099, },
	/* T^2 */ { 303, 698, 1345, },
	/* T^3 */ { 418, 744, 1584, },
	/* T^4 */ { 500, 777, 1741, },
	/* T^5 */ { 566, 803, 1805, },
	/* end */ { 592, 811, 1833, },
	/* termination */ { 0, }
};
static int vn_uranus[][3] = {
	/* addresses for uranus l, b, r  */
	/* T^0 */ { 0, 355, 462, },
	/* T^1 */ { 152, 398, 832, },
	/* T^2 */ { 253, 432, 1122, },
	/* T^3 */ { 315, 449, 1314, },
	/* T^4 */ { 347, 460, 1369, },
	/* T^5 */ { 354, 462, 1380, },
	/* end */ { 355, 0, 0, },
	/* termination */ { 0, }
};
static int vn_neptune[][3] = {
	/* addresses for neptune l, b, r  */
	/* T^0 */ { 0, 125, 182, },
	/* T^1 */ { 58, 148, 356, },
	/* T^2 */ { 90, 165, 462, },
	/* T^3 */ { 112, 175, 534, },
	/* T^4 */ { 122, 180, 556, },
	/* T^5 */ { 124, 181, 563, },
	/* end */ { 125, 182, 0, },
	/* termination */ { 0, }
};

int vsop87( double mjd, int obj, double prec, double *ret )
{
   static int (*vn_map[])[3] = {		/* indexes */
		vn_mercury, vn_venus, vn_earth, vn_mars, vn_jupiter,
		vn_saturn, vn_uranus, vn_neptune,
   };
   static double a0[] = {	/* semimajor axes; for precision ctrl only */
	    0.39, 0.72, 1.0, 1.5, 5.2, 9.6, 19.2, 30.1, 39.5,
   };
#include"vsop87_size.h"
#ifdef USE_PRECOMPILED_SIZE
   static int vn_size[] = {
      6648,4056,9744,16584,25344,43992,33120,13512,
   };
#endif
   double t[VSOP_MAXALPHA+1];			/* powers of time */
   double t_abs[VSOP_MAXALPHA+1];		/* powers of abs(time) */
   double q;					/* aux for precision control */
   int i, cooidx, alpha;			/* misc indexes */

   int (*vn_obj)[3] = vn_map[obj-1];
      
   /* load VSOP87 data into memory */
   char pfad[255];
   FILE *PDataFile;
   double *PData[8];
   int planet;
   
   if ((obj<=0)||(obj>=9))
   {
      fprintf(stderr,"%s\n%s\n%s\n","VSOP87 Planet ungueltig","","Programmende ...");
      exit(3);
   }
   
   strcpy(pfad,"VSOP87.DAT");
   PDataFile = fopen (pfad, "rb");
   if (PDataFile == NULL)
   {
      fprintf(stderr,"%s\n%s\n%s\n","Datei nicht gefunden:",pfad,"Programmende ...");
      exit(2);
   }
   for(planet=0;planet<8;planet++)
   {
      PData[planet]=(double*)malloc(vn_size[planet]);
      if (PData[planet] == NULL)
      {
	 fprintf(stderr,"%s\n%s\n%s\n","Kein Speicher frei fuer Planetenterme!","","Programmende ...");
	 exit(3);
      }
      if(fread(PData[planet],vn_size[planet],1,PDataFile)!=1)
      {
	 fprintf(stderr,"%s\n%s\n%s\n","Lesefehler:",pfad,"Programmende ...");
	 exit(2);
      }
   }
   fclose(PDataFile); 

   if (prec < 0.0 || prec > 1e-3)
       return(3);

   /* zero result array */
   for (i = 0; i < 6; ++i) ret[i] = 0.0;

   /* time and its powers */
   t[0] = 1.0;
   t[1] = (mjd - MJD_J2000)/VSOP_A1000;
   for (i = 2; i <= VSOP_MAXALPHA; ++i) t[i] = t[i-1] * t[1];
   t_abs[0] = 1.0;
   for (i = 1; i <= VSOP_MAXALPHA; ++i) t_abs[i] = fabs(t[i]);

   /* precision control */
   q = -log10(prec + 1e-35) - 2;	/* decades below 1e-2 */
   q = VSOP_ASCALE * prec / 10.0 / q;	/* reduce threshold progressively
					* for higher precision */

   /* do the term summation; first the spatial dimensions */
   for (cooidx = 0; cooidx < 3; ++cooidx)
   {
       /* then the powers of time */
       for (alpha = 0; vn_obj[alpha+1][cooidx] ; ++alpha)
       {
	   double p, term, termdot;

	   /* precision threshold */
	   p = q/(t_abs[alpha] + alpha * t_abs[alpha-1] * 1e-4 + 1e-35);
	   if (cooidx == 2)	/* scale by semimajor axis for radius */
	       p *= a0[obj];

	   term = termdot = 0.0;
	   for (i = vn_obj[alpha][cooidx]; i < vn_obj[alpha+1][cooidx]; ++i)
	   {
	       double a, b, c, arg;

	       a = PData[obj-1][3*i+0];
	       if (a < p) continue;	/* ignore small terms */

	       b = PData[obj-1][3*i+1];
	       c = PData[obj-1][3*i+2];

	       arg = b + c * t[1];
	       term += a * cos(arg);
#if VSOP_GETRATE
	       termdot += -c * a * sin(arg);
#endif
	   }

	   ret[cooidx] += t[alpha] * term;
#if VSOP_GETRATE
	   ret[cooidx + 3] += t[alpha] * termdot +
		   ((alpha > 0) ? alpha * t[alpha - 1] * term : 0.0);
#endif
       } /* alpha */
   } /* cooidx */

   for (i = 0; i < 6; ++i) ret[i] /= VSOP_ASCALE;

   /* reduce longitude to 0..2pi */
   ret[0] -= floor(ret[0]/(2.*PI)) * (2.*PI);

#if VSOP_GETRATE
   /* convert millenium rate to day rate */
   for (i = 3; i < 6; ++i) ret[i] /= VSOP_A1000;
#endif

   /* reduction from dynamical equinox of VSOP87 to FK5;
    */
   if (prec < 5e-7) {		/* 5e-7 rad = 0.1 arc seconds */
       double L1, c1, s1;
       L1 = ret[0] - DEG2RAD*(13.97 * t[1] - 0.031 * t[2]);
       c1 = cos(L1); s1 = sin(L1);
       ret[0] += DEG2RAD*(-0.09033 + 0.03916 * (c1 + s1) * tan(ret[1]))/3600.0;
       ret[1] += DEG2RAD*(0.03916 * (c1 - s1))/3600.0;
   }

   for(planet=0;planet<8;planet++)
   {
      free(PData[planet]);
   }

   return (0);
}
