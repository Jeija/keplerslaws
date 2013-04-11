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

#if !defined (_ASTROMAT_H)
#define _ASTROMAT_H

#include "cartes.h"

// check if PI is defined ...
#include <math.h>
#if !defined PI
// Mathematische Konstanten
const double PI=3.14159265358979323846264338327950;
#endif

// Astronomische Konstanten
const double J2000  =2451545.000;   // JD des Aequinoktium J2000
const double B1950  =2433282.423;   // JD des Aequinoktium B1950

// IAU(1976) System of Astronomical Constants
const double KGAUSS     =0.01720209895;      // Gausssche Gravitationskonstante
const double KGAUSS2    =0.0002959122083;    // KGAUSS*KGAUSS
const double KABERRATION=20.49552;           // Aberrationskonstante
const double AE_M       =149597870000.0;     // AE in Metern
const double RHO_M      =6378140.0;          // Erd-Aequator-Radius in Metern
const double LZAE       =499.004782/86400.0; // Lichtlaufzeit fuer 1 AE in d

const double RHO_AE =RHO_M/AE_M;    // Erd-Aequator-Radius in AE

// Umrechnungsfaktoren
const double M_2PI  =PI+PI;         // 2*pi
const double PI_2   =0.5*PI;        // 0.5*pi
const double DEG2RAD=PI/180.0;      // grad*DEG2RAD=bogenmass;
const double RAD2DEG=180.0/PI;      // bogenmass*RAD2DEG=grad;
const double B2RAD  =PI/648000.0;   // bogensekunden*B2RAD=bogenmass;
const double RAD2B  =648000.0/PI;   // bogenmass*RAD2B=bogensekunden;
const double RAD2H  =12.0/PI;       // bogenmass*RAD2H=dezimalstunden;
const double H2RAD  =PI/12.0;       // dezimalstunden*H2RAD=bogenmass;
const double ST2SO  =0.99726956633; // sternzeit*ST2SO=sonnenzeit;
const double SO2ST  =1.00273790935; // sonnenzeit*SO2ST=sternzeit;

// Funktions-Prototypen

double frac( double x );
double mod( double x, double y );
double root3( double x );

vektor helio2geo( const vektor& helio, const vektor& sonne );
vektor geo2helio( const vektor& geo, const vektor& sonne );
matrix m_ekl2equ( double ekls );
matrix m_equ2ekl( double ekls );
matrix m_equ2hor( double theta, double phi );
matrix m_hor2equ( double theta, double phi );
vektor geo2topo( const vektor& geo, const vektor& ort );
vektor topo2geo( const vektor& topo, const vektor& ort );

double ekls_m( double jde );
vektor geozentrischer_ort( double theta, double phi_graph, double h );
void annual_aberration_ekl( double jde, double lam, double bet,
                            double lam_sonne,
                            double* d_lam, double* d_bet );
void nutation_ekl( double jde, double* d_lam, double* d_eps );
matrix m_nutat_equ( double jde, double ekls );
matrix m_praez_ekl( double jd0, double jd );
matrix m_praez_equ( double jd0, double jd );

double juldat( short y, short m, short d, double t );
void kaldat( double jd, short* Y, short* M, short* D, double* T );
double lmst( double jd, double laenge );
double lmst2ut( double jd, double laenge, double lmst );

#endif /* !defined (_ASTROMAT_H) */
