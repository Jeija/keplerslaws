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
/* Module: ASTROMAT.CPP                                                      */
/* Version 1.0                                                               */
/* Last modified: March 15, 1993                                             */
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

/********************************/
/* Zusaetzliche Mathefunktionen */
/********************************/
/* frac()                       */
/* mod()                        */
/* root3()                      */
/********************************/

double frac( double x )
/* Dezimalteil von x mit Vorzeichen
   Eingabe:
   x
   Return:
   Dezimalteil von x mit Vorzeichen
*/
{  double dummy;

   return(modf(x,&dummy));
}

double mod( double x, double y )
/* x auf 0..y reduzieren
   Eingabe:
   x: Zahl die reduziert werden soll
   y: Zahl soll auf Bereich 0..y reduziert werden
   Return:
   x auf 0..y reduziert
*/
{  double t=fmod(x,y);
   if (t>=0) return t; else return t+y;
}

double root3( double x )
/* Dritte Wurzel von x
   Eingabe:
   x
   Return:
   Dritte Wurzel von x
*/
{  double e;

   if (x==0.0) return(0.0);
   else
   {  e=pow(fabs(x),1.0/3.0);
      if (x<0.0) e=-e;
      return e;
   }
}

/*******************************************/
/* Transformationen der Koordinatensysteme */
/*******************************************/
/* helio2geo()                             */
/* geo2helio()                             */
/* m_ekl2equ()                             */
/* m_equ2ekl()                             */
/* m_equ2hor()                             */
/* m_hor2equ()                             */
/* geo2topo()                              */
/* topo2geo()                              */
/*******************************************/

vektor helio2geo( const vektor& helio, const vektor& sonne )
/* heliozentrisch -> geozentrisch
   (gleiches Aequinoktium und Bezugssystem)
   Eingabe:
   helio: heliozentrische Planetenkoordinaten
   sonne: geozentrische Sonnenkoordinaten
   Return:
   geozentrische Planetenkoordinaten
*/
{
   return helio + sonne;
}

vektor geo2helio( const vektor& geo, const vektor& sonne )
/* geozentrisch -> heliozentrisch
   (gleiches Aequinoktium und Bezugssystem)
   Eingabe:
   geo:   geozentrische Planetenkoordinaten
   sonne: geozentrische Sonnenkoordinaten
   Return:
   heliozentrische Planetenkoordinatern
*/
{
   return geo - sonne;
}

matrix m_ekl2equ( double ekls )
/* Matrix geozentrisch ekliptikal -> geozentrisch aequatorial
   Eingabe:
   ekls: Winkel unter dem sich Aequator und Ekliptik schneiden
   Return:
   Matrix nach PAwyC p.50
*/
{  double s_eps=sin(ekls),c_eps=cos(ekls);

   return matrix( +1.0,   +0.0,   +0.0,
		  +0.0, +c_eps, -s_eps,
		  +0.0, +s_eps, +c_eps );
}

matrix m_equ2ekl( double ekls )
/* Matrix geozentrisch aequatorial -> geozentrisch ekliptikal
   Eingabe:
   ekls: Winkel unter dem sich Aequator und Ekliptik schneiden
   Return:
   Matrix nach PAwyC p.50
*/
{  double c_eps=cos(ekls),s_eps=sin(ekls);

   return matrix( +1.0,   +0.0,   +0.0,
		  +0.0, +c_eps, +s_eps,
		  +0.0, -s_eps, +c_eps );
}

matrix m_equ2hor( double theta, double phi )
/* Matrix aequatorial -> horizontal
   Der Suedpunkt hat das Azimut 0
   Eingabe:
   theta: Ortsssternzeit in RAD
   phi  : geographische Breite (Noerdlich > 0) in RAD
   Return:
   Matrix nach PAwyC p.50
*/
{  double cos_th=cos(theta),cos_phi=cos(phi),
	  sin_th=sin(theta),sin_phi=sin(phi);

   /* tau,delta -> a,h */
   matrix A( +sin_phi, +0.0, -cos_phi,
	     	 +0.0, +1.0,     +0.0,
	     +cos_phi, +0.0, +sin_phi );

   /* alpha,delta -> tau,delta */
   matrix B( +cos_th, +sin_th, +0.0,
	     +sin_th, -cos_th, +0.0,
		+0.0,    +0.0, +1.0 );

   /* alpha,delta -> a,h */
   return A*B;
}

matrix m_hor2equ( double theta, double phi )
/* Matrix horizontal -> aequatorial
   Der Suedpunkt hat das Azimut 0
   Eingabe:
   theta: Ortsssternzeit in RAD
   phi  : geographische Breite (Noerdlich > 0) in RAD
   Return:
   Matrix nach PAwyC p.50
*/
{  double cos_th=cos(theta),cos_phi=cos(phi),
	  sin_th=sin(theta),sin_phi=sin(phi);

   matrix A( +sin_phi, +0.0, +cos_phi,
		 +0.0, +1.0,     +0.0,
	     -cos_phi, +0.0, +sin_phi );

   matrix B( +cos_th, +sin_th, +0.0,
	     +sin_th, -cos_th, +0.0,
		+0.0,    +0.0, +1.0);

   return B*A;
}

vektor geo2topo( const vektor& geo, const vektor& ort )
/* geozentrisch aequatorial -> topozentrische aequatorial
   Eingabe:
   geo: geozentrische aequatoriale Planetenkoordinaten
   ort: geozentrische aequatoriale Koordinaten des Beobachters
   Return:
   topozentrische aequatoriale Planetenkoordinaten
*/
{
   return geo - ort;
}

vektor topo2geo( const vektor& topo, const vektor& ort )
/* topozentrisch aequatorial -> geozentrische aequatorial
   Eingabe:
   topo: topozentrische aequatoriale Planetenkoordinatengeozentrische
   ort : geozentrische aequatoriale Koordinaten des Beobachters
   Return:
   geozentrische aequatoriale Planetenkoordinaten
*/
{
   return topo + ort;
}

/************************************/
/* Koordinatensysteme - Korrekturen */
/************************************/
/* ekls_m()                         */
/* geozentrischer_ort()             */
/* annual_aberration_ekl()          */
/* nutation_ekl()                   */
/* m_nutat_equ                      */
/* m_praez_ekl                      */
/* m_praez_equ                      */
/************************************/

double ekls_m( double jde )
/* Berechnet die mittlere Ekliptikschiefe in RAD
   Eingabe:
   jde: Julianisches Ephemeriden Datum
   Return:
   mittlere Ekliptikschiefe nach AmdPC S.15
*/
{  double d=(jde-J2000);

   return 23.43929111*DEG2RAD-(46.8150*B2RAD/36525.0
			     +(0.00059*B2RAD/36525.0/36525.0
			     -0.001813*B2RAD/36525.0/36525.0/36525.0*d)*d)*d;
}

vektor geozentrischer_ort( double theta, double phi_graph, double h )
/* Geozentrische aequatoriale KARTESICHE Koordinaten des Beobachters
   auf dem Meeresspiegelniveau
   Eingabe:
   phi_graph: geographische Breite
   theta    : Ortssternzeit
   h        : Hoehe in Metern ueber NN
   Return:
   ort      : geozentrische aequatoriale Koordinaten des Beobachters
*/
{  double rcpz,rspz,u,h_,f=(1.0/298.257);

   h_=h/6378140.0;
   u=atan((1.0-f)*tan(phi_graph));

   rspz= (1 - f) * sin (u) + h_ * sin (phi_graph);
   rcpz= cos (u) + h_ * cos (phi_graph);

   return vektor( rcpz*RHO_AE * cos(theta),
		  rcpz*RHO_AE * sin(theta),
		  rspz*RHO_AE );
}

void annual_aberration_ekl( double jde, double lam, double bet,
                            double lam_sonne,
                            double* d_lam, double* d_bet )
/* the component of stellar aberration resulting from the
   motion of the Earth about the Sun

   Eingabe:
   jde      : Julianisches Ephemeridendatum
   lam      : ekliptikale Laenge
   bet      : ekliptikale Breite
   lam_sonne: geozentrisch ekliptikale Laenge der Sonne (zur Zeit jde)
              [geometrisch]
   Ausgabe:
   d_lam    : Korrektur zur ekliptikalen Laenge
   d_bet    : Korrektur zur ekliptikalen Breite
*/
{  double T,e,pi;

   T=(jde-J2000)/36525.0;
   e=0.016708617-0.000042037*T-0.0000001236*T*T;
   pi=(102.93735+0.71953*T+0.00046*T*T)*DEG2RAD;

   (*d_lam)=(-KABERRATION*B2RAD*cos(lam_sonne-lam)
             +e*KABERRATION*B2RAD*cos(pi-lam))/cos(bet);
   (*d_bet)=-KABERRATION*B2RAD*sin(bet)*(sin(lam_sonne-lam)-e*sin(pi-lam));
}

void nutation_ekl( double jde, double* d_lam, double* d_eps )
/* the short-period oscillations in the motion of the pole of rotation
   of the Earth that is undergoing torque from external gravitational forces

   Eingabe:
   jde  : Julianisches Ephemeridendatum
   Ausgabe:
   d_lam: Korrektur zur ekliptikalen Laenge
   d_eps: Korrektur zur mittleren Ekliptikschiefe
*/
{  double d=(jde-J2000),l,m,k,L,M,k2,L2,l2;

   l=218.316*DEG2RAD+481267.881*DEG2RAD/36525.0*d;   l2=2.0*l;
   m=134.963*DEG2RAD+477198.867*DEG2RAD/36525.0*d;
   k=125.045*DEG2RAD-  1934.136*DEG2RAD/36525.0*d;   k2=2.0*k;
   L=280.466*DEG2RAD+ 36000.770*DEG2RAD/36525.0*d;   L2=2.0*L;
   M=357.528*DEG2RAD+ 35999.050*DEG2RAD/36525.0*d;

   (*d_lam)=-17.200*B2RAD*sin(k) +0.206*B2RAD*sin(k2)
	    - 1.319*B2RAD*sin(L2)+0.143*B2RAD*sin(M)
	    - 0.227*B2RAD*sin(l2)+0.071*B2RAD*sin(m);

   (*d_eps)=  9.203*B2RAD*cos(k) -0.090*B2RAD*cos(k2)
	    + 0.574*B2RAD*cos(L2)+0.022*B2RAD*cos(L2+M)
	    + 0.098*B2RAD*cos(l2)+0.020*B2RAD*cos(l2-k);
}

matrix m_nutat_equ( double jde, double ekls )
/* Nutationsmatrix fuer aequatorialen Koordinaten
   Eingabe:
   ekls: mittlere Ekliptikschiefe
   jde : JDE fuer die wahren aequatorialen Koordinaten
   Return:
   Nutationsmatrix nach AA1991 B20

   Kombination mit Praezessionsmatrix:
   (matrix)PN=m_nut_equ(jd)*m_praez_equ(jd0,jd);
*/
{  double c_eps,s_eps,d_lam,d_eps;

   nutation_ekl(jde,&d_lam,&d_eps);
   c_eps=d_lam*cos(ekls);
   s_eps=d_lam*sin(ekls);

   return matrix( +1.0,   -c_eps, -s_eps,
		  +c_eps,   +1.0, -d_eps,
		  +s_eps, +d_eps,   +1.0 );
}

matrix m_praez_ekl( double jd0, double jd )
/* Praezessionsmatrix zur Transformation ekliptikaler Koordinaten (KARTESISCH)
   es wird das FK5 System verwendet
   Eingabe:
   jd0: Ausgangsaequinoktium
   jd : Neues Aequinoktium
   Return:
   Praezessionsmatrix nach AmdPC S.20
*/
{  double T =(jd -jd0  )/36525.0,
	  T0=(jd0-J2000)/36525.0,
	  Pi,pi,p,c1,c2,c3,s1,s2,s3;

   Pi=629554.98*B2RAD+3289.48*B2RAD*T0+0.61*B2RAD*T0*T0
     +( -869.81*B2RAD-0.5  *B2RAD*T0)*T+0.04 *B2RAD*T*T;
   pi=(  47.003*B2RAD-0.067*B2RAD*T0)*T-0.033*B2RAD*T*T;
   p =(5029.097*B2RAD+2.222*B2RAD*T0)*T+1.111*B2RAD*T*T;

   c1=cos(Pi+p); c2=cos(pi); c3=cos(Pi);
   s1=sin(Pi+p); s2=sin(pi); s3=sin(Pi);

   return matrix( +c1*c3+s1*c2*s3, +c1*s3-s1*c2*c3, -s1*s2,
		  +s1*c3-c1*c2*s3, +s1*s3+c1*c2*c3, +c1*s2,
			   +s2*s3,          -s2*c3,    +c2 );
}

matrix m_praez_equ( double jd0, double jd )
/* Berechnet die Praezessionsmatrix fuer aequatoriale Koordinaten (KARTESISCH)
   es wird das FK5 System verwendet
   Eingabe:
   jd0: Ausgangsaequinoktium
   jd : Neues Aequinoktium
   Return:
   Praezessionsmatrix nach AmdPC S.20
   Kombination mit Nuatationsmatrix:
   (matrix)PN=m_nut_equ(jd)*m_praez_equ(jd0,jd);
*/
{  double T =(jd -jd0  )/36525.0,
	  T0=(jd0-J2000)/36525.0,
	  zeta,z,theta,c1,c2,c3,s1,s2,s3;

   zeta =(2306.218*B2RAD+1.397*B2RAD*T0)*T+0.302*B2RAD*T*T+0.018*B2RAD*T*T*T;
   theta=(2004.311*B2RAD-0.853*B2RAD*T0)*T-0.427*B2RAD*T*T-0.042*B2RAD*T*T*T;
   z    =zeta+0.793*B2RAD*T*T;

   c1=cos(z); c2=cos(theta); c3=cos(zeta);
   s1=sin(z); s2=sin(theta); s3=sin(zeta);

   return matrix( -s1*s3+c1*c2*c3, -s1*c3-c1*c2*s3, -c1*s2,
		  +c1*s3+s1*c2*c3, +c1*c3-s1*c2*s3, -s1*s2,
		           +s2*c3,          -s2*s3,    +c2 );

}

/****************/
/* Zeitrechnung */
/****************/
/* juldat()     */
/* kaldat()     */
/* lmst()       */
/* lmst2ut()    */
/****************/

double juldat( short y, short m, short d, double t )
/* Berechnet das Julianische Datum (jd)
   Eingabe:
   y,m,d: Jahr, Monat, Tag
   t    : Weltzeit in Dezimalstunden (double)
   Return:
   JD nach Sky & Telescope (august, 1991, p.183)
   ACHTUNG:
   Es gilt die ASTRONOMISCHE Jahresrechnung (mit Jahr 0)
   Fehler:
   - wenn (d.m.y < 1.1.-4712)
   - wenn (04.10.1582 < d.m.y < 15.10.1582) falsches JD (Julianisch)
*/
{  signed long JD,Y=y,M=m,D=d;

   if (y>1582) /* Sicher gregorianischer Kalender */
   {  JD=367*Y-7*(Y+(M+9)/12)/4-3*((Y+(M-9)/7)/100+1)/4+275*M/9+D+1721029L;
      return (double)JD+t/24.0-0.5;
   }
   else if (y<1582) /* Sicher Julianischer Kalender */
   {  JD=367*Y-7*(Y+5001+(M-9)/7)/4+275*M/9+D+1729777L;
      return (double)JD+t/24.0-0.5;
   }
   else if ((m*100+d)>=1015) /* d.m.y>=15.10.1582 Gregorianischer Kalender) */
   {  JD=367*Y-7*(Y+(M+9)/12)/4-3*((Y+(M-9)/7)/100+1)/4+275*M/9+D+1721029L;
      return (double)JD+t/24.0-0.5;
   }
   else /* Julianischer Kalender (Unsinn vom 4.10.1582 bis zum 15.10.1582) */
   {  JD=367*Y-7*(Y+5001+(M-9)/7)/4+275*M/9+D+1729777L;
      return (double)JD+t/24.0-0.5;
   }
}

void kaldat( double jd, short* Y, short* M, short* D, double* T )
/* Berechnet aus dem Julianischem Datum (jd) das Buergerliche Datum
   Eingabe:
   jd     : Julianisches Datum
   Ausgabe:
   Y,M,D,T: Jahr, Monat, Tag, Weltzeit (in Dezimalstunden) nach GdE S.51
*/
{  unsigned long a,b,c,d,e,f;
   a=(unsigned long)floor(jd+0.5);

   if (a>=2299161l)
   {  b=(unsigned long)floor((a-1867216.25)/36524.25);
      c=a+b-(unsigned long)floor(b/4.0)+1525;
   }
   else
      c=a+1524l;
   d=(unsigned long)floor((c-122.1)/365.25);
   e=(unsigned long)floor(365.25*d);
   f=(unsigned long)floor((c-e)/30.6001);
   (*T)=frac(jd+0.5)*24.0;
   (*D)=c-e-(unsigned long)floor(30.6001*f);
   (*M)=f-1l-12l*(unsigned long)floor(f/14.0);
   (*Y)=d-4715l-(unsigned long)floor((7.0+(*M))/10.0);
}

double lmst( double jd, double laenge )
/* Berechnet die mittlere Ortssternzeit
   Eingabe:
   jd    : Julianisches Datum
   laenge: Geographische Laenge des Ortes in RAD
   (laenge < 0 wenn oestlich von Greenwich)
   Return:
   mittlere Ortssternzeit in RAD (local mean siderial time) nach GdE S.55
*/
{  double UT,jd0,gmst;

   UT  =frac (jd+0.5)*24.0;  // Weltzeit in Stunden
   jd0 =floor(jd+0.5)-0.5;   // JD um 0 Uhr Weltzeit
   gmst=15.0*DEG2RAD*6.656306+15.0*DEG2RAD*0.0657098242*(jd0-2445700.5)
			     +15.0*DEG2RAD*1.0027379093*UT;

   return mod(gmst-laenge,M_2PI);
}

double lmst2ut( double jd, double laenge, double sz )
/* LMST -> Weltzeit (an einem gegebenen Datum)
   Eingabe:
   jd    : (Beliebiges) Julianisches Datum am Tage der Umwandlung
   laenge: geographische Laenge in RAD (laenge < 0 wenn oestlich von Greenwich)
   sz    : mittlere Ortssternzeit in RAD
   Return:
   Weltzeit, die der LMST entspricht in RAD
   ACHTUNG:
   Die Umwandlung ist nicht eindeutig, da ein Sterntag kuerzer als ein
   Sonnentag ist. Nur ein Ergebnis wird dann zurueckgeliefert.
*/
{  double jd0,gmst0;

   jd0 =floor(jd+0.5)-0.5;                          // JD um 0 Uhr Weltzeit
   gmst0=15.0*DEG2RAD*6.656306
	+15.0*DEG2RAD*0.0657098242*(jd0-2445700.5); // GMST(0h UT)

   return mod((sz+laenge-gmst0),M_2PI)*0.9972695663;
}
