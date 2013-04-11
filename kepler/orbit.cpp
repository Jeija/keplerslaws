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
/* Module: ORBIT.CPP                                                         */
/* Version 3.0                                                               */
/* Last modified: 10. 3.1996                                                 */
// Probleme bei Komet Hyakutake (e=0.999846,a=1494.30519,m=0)
// Die Zahl der erlaubten Iterationen wurde auf 25 erhoeht
/* BIG BUG: (behoben in Version 2.0)
   in kepler muss fuer die mittlere Anomalie M gelten:
   0 <= M < M_2PI
   (sonst arbeitet die Iteration in pos_ell() nicht korrekt!)
   Bei dieser Gelegenheit wurden Fehlermeldungen eingebunden, fuer den Fall
   dass die Zahl der max. erlaubten Iterationen erreicht ist.
*/
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
#include <stdio.h>
#include "astromat.h"
// #include "extvar.h"
// #include "gscreen.h"
// #include "screen.h"
#include "orbit.h"

static const double c_eps=1e-11; /* Genauigkeit bei Iteration */

/**************/
/* Funktionen */
/**************/
/* pos_ell()  */
/* pos_par()  */
/* pos_hyp()  */
/* gaussvek() */
/* bahn2ekl() */
/* kepler()   */
/**************/

void pos_ell( double  a,  double  exz, double m,
	      double* x,  double* y,
	      double* vx, double* vy )
/* Orts- und Geschwindigkeitsvektor (ELLIPSE)
   Eingabe:
   a    : Grosse Bahnhalbachse in AE
   exz  : Bahnexzentrizitaet
   m    : Mittlere Anomalie in RAD
   Ausgabe:
   x,y  : Ortsvektor in AE
   vx,vy: Geschwindigkeit in AE/d
*/
{  double e,f,r,k,w,cos_e,sin_e;
   short i;

   if (exz<0.8) e=m; else e=PI;

   f=e-exz*sin(e)-m;
   for(i=0;(fabs(f)>c_eps)&&(i<25);i++)
   {  e=e-f/(1.0-exz*cos(e));
      f=e-exz*sin(e)-m;
   }

   if (i==25)
   {
      char s[40];
      sprintf(s,"a=%.5f e=%.5f m=%.5f",a,exz,m);

      fprintf(stderr,"%s\n%s\n%s\n","Konvergenz-Problem in pos_ell()! Werte",
      "notieren und Autor benachrichtigen:",
      s);
   }

   cos_e=cos(e); sin_e=sin(e);
   k=KGAUSS/sqrt(a);
   w=sqrt((1.0-exz)*(1.0+exz));
   r=1.0-exz*cos_e;
   (*x)=a*(cos_e-exz);
   (*y)=a*w*sin_e;
   (*vx)=-k*sin_e/r;
   (*vy)=+k*cos_e/r*w;
}

void pos_par( double  q,  double  d, double d0,
	      double* x,  double *y,
	      double* vx, double *vy )
/* Orts- und Geschwindigkeitsvektor (PARABEL)
   Eingabe:
   q    : Periheldistanz in AE
   d    : Beobachtungszeitpunkt JD
   d0   : Perhiheldurchgangszeitpunkt JD
   Ausgabe:
   x,y  : Ortsvektor in AE
   vx,vy: Geschwindigkeit in AE/d
*/
{  double A,B,r,u,u2,sqrt_2q=sqrt(2.0*q),f=KGAUSS/sqrt_2q;

   A=1.5*KGAUSS/(q*sqrt_2q)*(d-d0);
   B=root3(A+sqrt(A*A+1.0));
   u=B-1.0/B;     /* tangens( wahre Anomalie/2 ) */
   u2=u*u;        /* u hoch 2                    */
   r=q*(1.0+u2);  /* Entfernung von der Sonne    */
   (*x)=q*(1.0-u2);
   (*y)=q*2.0*u;
   (*vx)=-(*y)/r*f;
   (*vy)=((*x)/r+1.0)*f;
}

void pos_hyp( double  a,  double  exz, double d, double d0,
	      double* x,  double* y,
	      double* vx, double* vy )
/* Orts- und Geschwindigkeitsvektor (HYPERBEL)
   Eingabe:
   a    : Grosse Bahnhalbachse in AE
   exz  : Bahnexzentrizitaet
   d    : Beobachtungszeitpunkt JD
   d0   : Perhiheldurchgangszeitpunkt JD
   Ausgabe:
   x,y  : Ortsvektor in AE
   vx,vy: Geschwindigkeit in AE/d
*/
{  double k,mh,h,f,cosh_h,sinh_h,w,r;
   short i;

   a=fabs(a);
   k=KGAUSS/sqrt(a);
   mh=k*(d-d0)/a;

   h=root3(6.0*mh); /* Startwert nach GdE S.70 */

/* h=log(2.0*fabs(mh)/exz+1.8); */ /* Startwert nach AmdPC S.61 */
/* if (mh<0) h=-h;              */

   f=exz*sinh(h)-h-mh;

   for(i=0;(fabs(f)>c_eps /* (1.0+fabs(h+mh)) */ )&&(i<25);i++)
   {  h=h-f/(exz*cosh(h)-1.0);
      f=exz*sinh(h)-h-mh;
   }

   if (i==25)
   {
      char s[40];
      sprintf(s,"a=%.5f e=%.5f d-d0=%.3f",a,exz,d-d0);

      fprintf(stderr,"%s\n%s\n%s\n","Konvergenz-Problem in pos_hyp()! Werte",
      "notieren und Autor benachrichtigen:",
      s);
   }

   cosh_h=cosh(h); sinh_h=sinh(h);
   w=sqrt((exz+1.0)*(exz-1.0));
   r=exz*cosh_h-1.0;
   (*x)=a*(exz-cosh_h);
   (*y)=a*w*sinh_h;
   (*vx)=-k*sinh_h/r;
   (*vy)=+k*cosh_h/r*w;
}

void gaussvek( double k, double i, double w, vektor& P, vektor& Q )
/* Berechnung der Gaussschen Vektoren aus ekliptikalen Bahnelementen
   (Aequinoktium wie k,i und w)
   Eingabe:
   k: Knotenlaenge des aufsteigenden Knotens in RAD
   i: Bahnneigung in RAD
   w: Argument des Perihels in RAD
   Ausgabe:
   P: Gaussscher Vektor
   Q: Gaussscher Vektor
   Anwendung:
   p_helio_ekl= x_bahn*P +y_bahn*Q; heliozentrische ekliptikale Koordinaten
   v_helio_ekl=vx_bahn*P+vy_bahn*Q; Geschwindigkeitsvektor
*/
{  double c1=cos(w), c2=cos(i), c3=cos(k),
	  s1=sin(w), s2=sin(i), s3=sin(k);

   P[0] = +c1*c3 - s1*c2*s3;
   P[1] = +c1*s3 + s1*c2*c3;
   P[2] = +s1*s2;

   Q[0] = -s1*c3 - c1*c2*s3;
   Q[1] = -s1*s3 + c1*c2*c3;
   Q[2] = +c1*s2;
}

vektor bahn2ekl( double x, double y, const vektor& P, const vektor& Q )
/* Transformation von Koordinaten im System der Bahnebene
   in ekliptikale Koordinaten (Aequinoktium wie das der Gaussschen Vektoren)
   Eingabe:
   x,y: Koordinaten (Geschwindigkeit) in der Bahnebene
   P,Q: Gausssche Vektoren
   Return:
   heliozentrische ekliptikale Koordinaten
*/
{
   return x * P + y * Q;
}

vektor kepler( double e, double q, double d, double d0,
               const vektor& P, const vektor& Q,
               vektor& v_helio_ekl )
/* Berechnet die Heliozentrische ekliptikale Position eines Gestirns aus
   Bahnelementen (Aequinoktium wie das der Gaussschen Vektoren)
   Eingabe:
   e  : Exzentrizitaet
   q  : Periheldistanz
   d  : JD Beobachtungszeitpunkt
   d0 : JD Periheldurchgang
   P,Q: Gausssche Vektoren
   Ausgabe:
   v_helio_ekl: Geschwindigkeitsvektor
   Return:
   p_helio_ekl: Heliozentrische ekliptikale Position
*/
{  double x,y,vx,vy;
   vektor p_helio_ekl;

   if (e<1.0)
   {  double a,M;
      a=q/(1.0-e);
      M=KGAUSS*(d-d0)*1.0/sqrt(a*a*a);
// BIG BUG 19. 4.1994: mod eingefuehrt! (Iteration in pos_ell versagte)
      M=mod(M,M_2PI);
      pos_ell(a,e,M,&x,&y,&vx,&vy);
   }
   else if (e==1.0)
   {  pos_par(q,d,d0,&x,&y,&vx,&vy);
   }
   else
   {  double a;
      a=q/(1.0-e);
      pos_hyp(a,e,d,d0,&x,&y,&vx,&vy);
   }
   p_helio_ekl= x*P+ y*Q;
   v_helio_ekl=vx*P+vy*Q;

   return p_helio_ekl;
}
