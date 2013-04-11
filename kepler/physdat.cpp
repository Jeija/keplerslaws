/***************************************************************************
                          physdat.cpp  -  description
                             -------------------
    begin                : Sat Apr 15 2000
    copyright            : (C) 2000 by tkramer
    email                : tkramer@ph.tum.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
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
/* Module: PHYSDAT.CPP                                                       */
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

static char rcsid[] = "$Id: physdat.cpp,v 1.1 2000/04/16 08:38:03 tkramer Exp $";

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "physdat.h"

/*****************************************************************************/
/* ACHTUNG: Alle Koordinaten sind auf den MOMENTANEN FRUEHLINGSPUNKT bezogen */
/*          es muss ausserdem das JDE anstelle des jd verwendet werden       */
/*	    zudem muss die LICHTLAUFZEIT beruecksichtigt werden              */
/*****************************************************************************/

/**************************/
/* Unabhaengig von planet */
/**************************/
/* elong_helio()          */
/* elongation()           */
/* phasenwin()            */
/* poswin_sonne()         */
/* phase()                */
/* beldef()               */
/* poswin_achse()         */
/* zm_zb()                */
/* hor_par()              */
/**************************/

double elong_helio( const vektor& p_h, const vektor& s_g )
/* Elongation Planet-Sonne im gleichen Koordinatensystem mit unterschiedlichem
   Ursprung
   Eingabe:
   p_h: Planet HELIOzentrisch
   s_g: Sonne GEOzentrisch
   Return:
   Elongation
*/
{  double S2 = (s_g|s_g),
	  P2 = (p_h|p_h),
	  PS = (p_h|s_g);

   return acos( (S2+PS) / sqrt( S2 * (S2 + P2 + 2.0 * PS ) ) );
}

double elongation( const vektor& p, const vektor& s )
/* Elongation Planet-Sonne (beide im GLEICHEM Koordinatensystem)
   Eingabe:
   p: Planet (beide im GLEICHEM Koordinatensystem)
   s: Sonne  (beide im GLEICHEM Koordinatensystem)
   Return:
   Elongation
*/
{
   return acos( (p|s) / sqrt( (p|p) * (s|s) ) );
}

double phasenwin( const vektor& p_h, const vektor& p_g )
/* Phasenwinkel
   Eingabe:
   p_h: Planet (heliozentrisch)
   p_g: Planet (geozentrisch)
   Return:
   Phasenwinkel
*/
{
   return acos( (p_h|p_g) / sqrt( (p_h|p_h) * (p_g|p_g) ) );
}

double poswin_sonne( const vektor& p_g, const vektor& s_g )
/* Positionswinkel der Sonne
   Eingabe:
   p_g: Planet geozentrisch aequatorial
   s_g: Sonne geozentrisch aequatorial
   Return:
   Positionswinkel der Sonne
*/
{  double P=1.0/sqrt(p_g|p_g);

   return atan2( P * (s_g[1] * p_g[0] - s_g[0] * p_g[1]),
                 s_g[2] - P * P * p_g[2] * (s_g|p_g) );
}

double phase( double i )
/* Phase (Anteil der beleuchteten Flaeche an der Gesamtflaeche)
   Eingabe:
   i: Phasenwinkel
   Return:
   Phase
*/
{
   return 0.5 * ( 1.0 + cos(i) );
}

double beldef( double i, double d )
/* Beleuchtungsdefekt
   Eingabe:
   i: Phasenwinkel
   d: (Scheinbarer) Durchmesser des Planeten
   Return:
   Beleuchtungsdefekt
*/
{
   return d * 0.5 * ( 1.0 - cos(i) );
}

double poswin_achse( const vektor& p_g, const vektor& a )
/* Positionswinkel der Rotationsachse
   Eingabe:
   p: Planet geozentrisch aequatorial
   a: Nordpol der Rotationsachse aequatorial
   Return:
   Positionswinkel der Rotationsachse
*/
{  double P=1.0/sqrt(p_g|p_g);

   return atan2(P * (a[1] * p_g[0] - a[0] * p_g[1]),
                a[2] - P * P * p_g[2] * (a|p_g) );
}

double zm_zb( const vektor& p, const vektor& a, double w, short s, double* zb )
/* Zentralmeridian und planetozentrische Erdbreite
   Eingabe:
   p: Planet geozentrisch aequatorial
   a: Nordpol der Rotationsachse aequatorial
   w: Laenge des Nullmeridians
   s: Rotationssinn des Planeten (W'(t)<0 oder W'(t)>0)
   Return:
   Zentralmeridian
   Ausgabe:
   zb:planetoZENTRISCHE Erdbreite
*/
{  double ap=(a|p),a2=(a|a),p2=(p|p),K;

   (*zb)=-asin( ap / sqrt(a2*p2) );
   K=atan2( ap * a[2] - a2 * p[2],
            (p[0] * a[1] - p[1] * a[0]) );
   if (s>0) return mod( w - K,M_2PI);
   else     return mod( K - w,M_2PI);
}

double hor_par( double r )
/* Eingabe:
   r: geozentrischer Abstand eines Planeten von der Erde (in AE)
   Return:
   Horizontparallaxe in RAD
*/
{
   double x=RHO_AE/r;

   return (x<1.0) ? asin(RHO_AE/r) : -99.9 ;
}

/**************************/
/* Abhaengig von planet   */
/**************************/
/* rot_sinn()             */
/* libration()            */
/* durchmesser()          */
/* helligkeit()           */
/* iau1982_j2000()        */
/* iau1988_j2000()        */
/* iau_datum()            */
/* rot_achse_null_mer()   */
/**************************/

short rot_sinn( PLANET planet )
/* Rotationssinn des Planeten (W(t)'>0) oder (W(t)'<0)
   Eingabe:
   planet: Planet
   Return:
   Rotationssinn
*/
{  switch(planet)
   {  
      case so :	return -1;
      case me : return +1;
      case ve : return -1;
      case ma : return +1;
      case ju : return +1;
      case sa : return +1;
      case ur : return -1;
      case ne : return +1;
      case pl : return -1;
      case mo : return -1;
      default : printf("\nNicht implementiert!!!"); exit(1);
   }
   return (1);
}

void libration( double jde, vektor m_ekl, double* l, double* b )
/* Libration des Mondes nach AsAl
   Eingabe:
   jde  : Julianisches Ephemeriden Datum
   m_ekl: Mond geozentrisch ekliptikal
   Ausgabe:
   l    : Selenographische Laenge des sub-earth point
   b    : Selenographische Breite des sub-earth point
*/
{  double T,k,f,sin_k,cos_k,sin_i,cos_i,r;
   const double I=1.54242*DEG2RAD; // IAU-value mean lunar equator-ecliptic

   T=(jde-J2000)/36525.0;
   k=125.04452*DEG2RAD-  1934.136261*DEG2RAD*T
    +0.0020708*DEG2RAD*T*T+T*T*T/450000.0*DEG2RAD;
   f= 93.27191*DEG2RAD+483202.017538*DEG2RAD*T
    -0.0036825*DEG2RAD*T*T+T*T*T/327270.0*DEG2RAD;
   sin_k=sin(k);  sin_i=sin(I);
   cos_k=cos(k);  cos_i=cos(I);

   r=sqrt(m_ekl|m_ekl);
   (*l)=atan2( m_ekl[1]*cos_k*cos_i-m_ekl[0]*sin_k*cos_i-m_ekl[2]*sin_i,
               m_ekl[0]*cos_k+m_ekl[1]*sin_k)-f;
   (*b)=asin((-m_ekl[1]*cos_k*sin_i+m_ekl[0]*sin_k*sin_i-m_ekl[2]*cos_i)/r);
}

double durchmesser( PLANET planet, double ep )
/* Scheinbarer Durchmesser eines Planeten
   Quelle: AA 1992 E43
   Eingabe:
   planet: Nummer des Planeten
   ep    : geozentrischen Entfernung (in AE)
*/
{  double s;
   switch(planet)
   {  case so   : s=959.63*B2RAD;   break;
      case me   : s=  3.36*B2RAD;   break;
      case ve   : s=  8.34*B2RAD;   break;
      case ma   : s=  4.68*B2RAD;   break;
      case ju   : /* Aequator */
      case ju_a : s= 98.44*B2RAD;   break;
      case ju_p : s= 92.06*B2RAD;   break;
      case sa   : /* Aequator */
      case sa_a : s= 82.73*B2RAD;   break;
      case sa_p : s= 73.82*B2RAD;   break;
      case sa_ri: s=124.8 *B2RAD;   break;
      case sa_ra: s=187.7 *B2RAD;   break;
      case ur   : s= 35.02*B2RAD;   break;
      case ne   : s= 33.50*B2RAD;   break;
      case pl   : s=  2.07*B2RAD;   break;
      case mo   : s=1738000.0/AE_M; break; // Radius (m) / AE (m)
      default   : printf("\nNicht implementiert!!!"); exit(1);
    }
   return(2.0*s/ep);
}

double helligkeit( PLANET planet, double sp, double ep, double i, double b )
/* Scheinbare Helligkeit eines Planeten
   (source: AsAl p.270 values used in the AA since 1984)
   Eingabe:
   planet: Planet
   sp    : Entfernung Sonne-Planet in AE
   ep    : Entfernung Erde-Planet in AE
   i     : Phasenwinkel in RAD
   b     : planetographische Erdbreite (nur bei Saturn)
   Return:
   scheinbare Helligkeit in Groessenklassen
*/
{  double m0,s_b,ig;

   ig=i*RAD2DEG; /* Phasenwinkel in Grad */

   switch (planet)
    { case me : m0=-0.42+0.0380*ig-0.000273*ig*ig+0.000002  *ig*ig*ig; break;
      case ve : m0=-4.40+0.0009*ig+0.000239*ig*ig-0.00000065*ig*ig*ig; break;
      case ma : m0=-1.52+0.016*ig; break;
      case ju : m0=-9.40+0.005*ig; break;
      case sa : s_b=sin(fabs(b));
		m0=-8.88+0.044*ig-2.60*s_b+1.25*s_b*s_b; break;
      case ur : m0=-7.19; break;
      case ne : m0=-6.87; break;
      case pl : m0=-1.00; break;
      default : printf("\nNicht implementiert!!!"); exit(1);
    }
   return(m0+5.0*log10(sp*ep));
}

void iau1982_j2000( PLANET planet, double jd,
                    double* alpha0, double* delta0, double* W )
/* Alle Daten direkt nach dem IAU 1982 report
   (Saturn System I nach GdE S.114)

   Eingabe:
   planet: Nummer des Planeten
   jd    : Julianisches Datum
   Ausgabe:
   alpha0: Rektaszension der Rotationsachse (J2000)
   delta0: Deklination der Rotationsachse (J2000)
   W     : Winkel des Nullmeridians (J2000)
*/
{  double T,d,e1,e2,e3,e4,e5;

   T=(jd-J2000)/36525.0;
   d=jd-J2000;

   switch (planet)
   {  case so    : (*alpha0)=285.96;
	           (*delta0)= 63.96;
                   (*W)     = 84.11+14.1844000*d;
                   break;
      case me    : (*alpha0)=281.02-0.003*T;
	           (*delta0)= 61.45-0.005*T;
                   (*W)     =329.71+6.1385025*d;
	           break;
      case ve    : (*alpha0)=272.78;
	           (*delta0)= 67.21;
                   (*W)     =159.91-1.4814205*d;
                   break;
      case ma    : (*alpha0)=317.681-0.108*T;
	           (*delta0)= 52.886-0.061*T;
                   (*W)     =176.655+350.8919830*d;
                   break;
      case ju_iii: ;
      case ju    : (*alpha0)=268.05-0.009*T;
	           (*delta0)= 64.49+0.003*T;
  /* System III */ (*W)     =284.95+870.5360000*d;
                   break;
      case ju_i  : (*alpha0)=268.05-0.009*T;
	           (*delta0)= 64.49+0.003*T;
  /* System I */   (*W)     = 67.10+877.900*d;
                   break;
      case ju_ii : (*alpha0)=268.05-0.009*T;
	           (*delta0)= 64.49+0.003*T;
  /* System II */  (*W)     = 43.30+870.270*d;
                   break;
      case sa_iii: ;
      case sa    : (*alpha0)= 40.66-0.036*T;
	           (*delta0)= 83.52-0.004*T;
  /* System III */ (*W)     = 38.90+810.7939024*d;
                   break;
      case sa_i  : (*alpha0)= 40.66-0.036*T;
	           (*delta0)= 83.52-0.004*T;
  /* System I */   (*W)     =227.2037+844.3000000*d;
                   break;
      case ur    : (*alpha0)=257.43;
	           (*delta0)=-15.10;
  /* System III */ (*W)     =261.62-554.9130000*d;
                   break;
      case ne    : (*alpha0)=295.33;
	           (*delta0)= 40.65;
                   (*W)     =107.21+468.7500000*d;
                   break;
      case pl    : (*alpha0)=311.63;
	           (*delta0)=  4.18;
                   (*W)     =252.66-56.3640000*d;
                   break;
      case mo    : e1=(125.045- 0.052992*d)*DEG2RAD;
	           e2=(249.390- 0.105984*d)*DEG2RAD;
	           e3=(196.694-13.012000*d)*DEG2RAD;
	           e4=(176.630+13.340716*d)*DEG2RAD;
                   e5=(358.219- 0.985600*d)*DEG2RAD;

	           (*alpha0)=270.000 -3.878*sin(e1)-0.120*sin(e2)
			             +0.070*sin(e3)-0.017*sin(e4);
	           (*delta0)= 66.534 +1.543*cos(e1)+0.024*cos(e2)
				     -0.028*cos(e3)+0.007*cos(e4);
		   (*W)     = 38.314+13.1763581*d  +3.558*sin(e1)
				     +0.121*sin(e2)-0.064*sin(e3)
				     +0.016*sin(e4)+0.025*sin(e5);
                   break;
      default    : printf("\nNicht implementiert!!!"); exit(1);
   }

   (*alpha0)*=DEG2RAD;
   (*delta0)*=DEG2RAD;
   (*W)     *=DEG2RAD;
}

void iau1988_j2000( PLANET planet, double jd,
                    double* alpha0, double* delta0, double* W )
/* Alle Daten direkt nach dem IAU 1988 report
   (Saturn System I nach GdE S.114)

   Eingabe:
   planet: Nummer des Planeten
   jd    : Julianisches Datum
   Ausgabe:
   alpha0: Rektaszension der Rotationsachse (J2000)
   delta0: Deklination der Rotationsachse (J2000)
   W     : Winkel des Nullmeridians (J2000)
*/
{  double T,d,N,e1,e2,e3,e4,e5;

   T=(jd-J2000)/36525.0;
   d=jd-J2000;

   switch (planet)
   {  case so    : (*alpha0)=286.13;
	           (*delta0)= 63.87;
                   (*W)     = 84.10+14.1844000*d;
                   break;
      case me    : (*alpha0)=281.01-0.003*T;
	           (*delta0)= 61.45-0.005*T;
                   (*W)     =329.71+6.1385025*d;
	           break;
      case ve    : (*alpha0)=272.69;
	           (*delta0)= 67.17;
                   (*W)     =160.39-1.4813291*d;
                   break;
      case ma    : (*alpha0)=317.681-0.108*T;
	           (*delta0)= 52.886-0.061*T;
                   (*W)     =176.868+350.8919830*d;
                   break;
      case ju_iii: ;
      case ju    : (*alpha0)=268.05-0.009*T;
	           (*delta0)= 64.49+0.003*T;
  /* System III */ (*W)     =284.95+870.5360000*d;
                   break;
      case ju_i  : (*alpha0)=268.05-0.009*T;
	           (*delta0)= 64.49+0.003*T;
  /* System I */   (*W)     = 67.10+877.900*d;
                   break;
      case ju_ii : (*alpha0)=268.05-0.009*T;
	           (*delta0)= 64.49+0.003*T;
  /* System II */  (*W)     = 43.30+870.270*d;
                   break;
      case sa_iii: ;
      case sa    : (*alpha0)= 40.58-0.036*T;
	           (*delta0)= 83.54-0.004*T;
  /* System III */ (*W)     = 38.90+810.7939024*d;
                   break;
      case sa_i  : (*alpha0)= 40.58-0.036*T;
	           (*delta0)= 83.54-0.004*T;
  /* System I */   (*W)     =227.2037+844.3000000*d;
                   break;
      case ur    : (*alpha0)=257.43;
	           (*delta0)=-15.10;
  /* System III */ (*W)     =203.81-501.1600928*d;
                   break;
      case ne    : N=(359.28+54.308*T)*DEG2RAD;
                   (*alpha0)=298.72+  2.58*sin(N)-0.04*sin(2.0*N);
	           (*delta0)= 42.63-  1.90*cos(N)+0.01*cos(2.0*N);
                   (*W)     =313.66+483.7625981*d-1.75*sin(N)+0.04*sin(2.0*N);
                   break;
      case pl    : (*alpha0)=311.63;
	           (*delta0)=  4.18;
                   (*W)     =252.66-56.3640000*d;
                   break;
      case mo    : e1=(125.045- 0.052992*d)*DEG2RAD;
	           e2=(249.390- 0.105984*d)*DEG2RAD;
	           e3=(196.694-13.012000*d)*DEG2RAD;
	           e4=(176.630+13.340716*d)*DEG2RAD;
                   e5=(358.219- 0.985600*d)*DEG2RAD;

	           (*alpha0)=270.000 -3.878*sin(e1)-0.120*sin(e2)
			             +0.070*sin(e3)-0.017*sin(e4);
	           (*delta0)= 66.534 +1.543*cos(e1)+0.024*cos(e2)
				     -0.028*cos(e3)+0.007*cos(e4);
		   (*W)     = 38.314+13.1763581*d  +3.558*sin(e1)
				     +0.121*sin(e2)-0.064*sin(e3)
				     +0.016*sin(e4)+0.025*sin(e5);
                   break;
      default    : printf("\nNicht implementiert!!!"); exit(1);
   }

   (*alpha0)*=DEG2RAD;
   (*delta0)*=DEG2RAD;
   (*W)     *=DEG2RAD;
}

void iau1991_j2000( PLANET planet, double jd,
                    double* alpha0, double* delta0, double* W )
/* Alle Daten direkt nach dem IAU 1991 report
   (Saturn System I nach GdE S.114)

   Eingabe:
   planet: Nummer des Planeten
   jd    : Julianisches Datum
   Ausgabe:
   alpha0: Rektaszension der Rotationsachse (J2000)
   delta0: Deklination der Rotationsachse (J2000)
   W     : Winkel des Nullmeridians (J2000)
*/
{  double T,d,N,e1,e2,e3,e4,e5;

   T=(jd-J2000)/36525.0;
   d=jd-J2000;

   switch (planet)
   {  case so    : (*alpha0)=286.13;
		   (*delta0)= 63.87;
		   (*W)     = 84.10+14.1844000*d;
		   break;
      case me    : (*alpha0)=281.01-0.003*T;
		   (*delta0)= 61.45-0.005*T;
		   (*W)     =329.71+6.1385025*d;
		   break;
      case ve    : (*alpha0)=272.76;
		   (*delta0)= 67.16;
		   (*W)     =160.20-1.4813688*d;
		   break;
      case ma    : (*alpha0)=317.681-0.108*T;
		   (*delta0)= 52.886-0.061*T;
		   (*W)     =176.868+350.8919830*d;
		   break;
      case ju_iii: ;
      case ju    : (*alpha0)=268.05-0.009*T;
		   (*delta0)= 64.49+0.003*T;
  /* System III */ (*W)     =284.95+870.5360000*d;
		   break;
      case ju_i  : (*alpha0)=268.05-0.009*T;
		   (*delta0)= 64.49+0.003*T;
  /* System I */   (*W)     = 67.10+877.900*d;
		   break;
      case ju_ii : (*alpha0)=268.05-0.009*T;
		   (*delta0)= 64.49+0.003*T;
  /* System II */  (*W)     = 43.30+870.270*d;
		   break;
      case sa_iii: ;
      case sa    : (*alpha0)= 40.58-0.036*T;
		   (*delta0)= 83.54-0.004*T;
  /* System III */ (*W)     = 38.90+810.7939024*d;
		   break;
      case sa_i  : (*alpha0)= 40.58-0.036*T;
		   (*delta0)= 83.54-0.004*T;
  /* System I */   (*W)     =227.2037+844.3000000*d;
		   break;
      case ur    : (*alpha0)=257.43;
		   (*delta0)=-15.10;
  /* System III */ (*W)     =203.81-501.1600928*d;
		   break;
      case ne    : N=(357.85+52.316*T)*DEG2RAD;
		   (*alpha0)=299.36+  0.70*sin(N);
		   (*delta0)= 43.46-  0.51*cos(N);
		   (*W)     =253.18+536.3128492*d-0.48*sin(N);
		   break;
      case pl    : (*alpha0)=313.02;
		   (*delta0)=  9.09;
		   (*W)     =236.77-56.3623195*d;
		   break;
      case mo    : e1=(125.045- 0.052992*d)*DEG2RAD;
		   e2=(250.090- 0.105984*d)*DEG2RAD;
		   e3=(260.008-13.012001*d)*DEG2RAD;
		   e4=(176.625+13.340716*d)*DEG2RAD;
		   e5=(357.529- 0.985600*d)*DEG2RAD;

		   (*alpha0)=270.000 +0.003*T
				     -3.878*sin(e1)-0.120*sin(e2)
				     +0.070*sin(e3)-0.017*sin(e4);
		   (*delta0)= 66.541 +0.013*T
				     +1.543*cos(e1)+0.024*cos(e2)
				     -0.028*cos(e3)+0.007*cos(e4);
		   (*W)     = 38.317+13.1763582*d  +3.558*sin(e1)
				     +0.121*sin(e2)-0.064*sin(e3)
				     +0.016*sin(e4)+0.025*sin(e5);
		   break;
      default    : printf("\nNicht implementiert!!!"); exit(1);
   }

   (*alpha0)*=DEG2RAD;
   (*delta0)*=DEG2RAD;
   (*W)     *=DEG2RAD;
}

void iau_datum( double jd, double alpha2000, double delta2000, double W2000,
                double* alpha, double* delta, double* W )
/* Daten direkt nach dem IAU report von J2000 auf das Datum transformieren

   Eingabe:
   planet   : Nummer des Planeten
   jd       : Julianisches Datum
   alpha2000: Rektaszension der Rotationsachse (J2000)
   delta2000: Deklination der Rotationsachse (J2000)
   W2000    : Nullmeridians (J2000)
   Ausgabe:
   alpha    : Rektaszension der Rotationsachse (des Datums)
   delta    : Deklination der Rotationsachse (des Datums)
   W        : Nullmeridians (des Datums)
*/
{  double T,zeta,theta,z,delta_W,X,Y,Z,q,r;

   T=(jd-J2000)/36525.0;
   zeta =(2306.2181*T+0.30188*T*T+0.017998*T*T*T)*B2RAD;
   z    =(2306.2181*T+1.09468*T*T+0.018203*T*T*T)*B2RAD;
   theta=(2004.3109*T-0.42665*T*T-0.041833*T*T*T)*B2RAD;

   X=cos(theta)*cos(delta2000)*cos(alpha2000+zeta)-sin(theta)*sin(delta2000);
   Y=cos(delta2000)*sin(alpha2000+zeta);
   Z=sin(theta)*cos(delta2000)*cos(alpha2000+zeta)+cos(theta)*sin(delta2000);

   q=X*X+Y*Y;
   if (q==0.0) (*alpha)=0.0; else (*alpha)=atan2(Y,X);
   r=sqrt(q+Z*Z);
   q=sqrt(q);
   if (r==0.0) (*delta)=0.0; else (*delta)=atan2(Z,q);
   (*alpha)+=z;
   delta_W=asin(-sin(theta)*sin(alpha2000+zeta)/cos((*delta)));
   (*W)=W2000+delta_W;
}

void rot_achse_null_mer( PLANET planet, double jd, vektor& a_equ, double* w )
/* Rotationsachse und Nullmeridian
   Eingabe:
   planet: Planet
   jd    : Julianisches Datum
   Ausgabe:
   a_equ : Rotationsachse (Aequinoktium des Datums)
   w     : Nullmeridian (Aequinoktium des Datums)
*/
{  double alpha2000,delta2000,W2000,alpha,delta;

/*
   iau1982_j2000(planet,jd,&alpha2000,&delta2000,&W2000);
   iau1988_j2000(planet,jd,&alpha2000,&delta2000,&W2000);
*/
   iau1991_j2000(planet,jd,&alpha2000,&delta2000,&W2000);

   iau_datum(jd,alpha2000,delta2000,W2000,&alpha,&delta,w);

   a_equ=(vektor)sphaer(alpha,delta,1.0);
}

double eins_minus_f_quadrat( PLANET planet )
{
   double f;

   switch (planet)
   {
      case so :	f=0.0; break;
      case me : f=0.0; break;
      case ve : f=0.0; break;
      case er : f=0.0; break;
      case ma : f=0.0065;  break;
      case ju : f=0.06487; break;
      case sa : f=0.09796; break;
      case ur : f=0.02293; break;
      case ne : f=0.0171;  break;
      case pl : f=0.0; break;
      case mo : f=0.0; break;
      default : printf("\nNicht implementiert!!!"); exit(1);
   }
   return (1.0-f)*(1.0-f);
}
