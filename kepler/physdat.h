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

#if !defined (_PHYSDAT_H)
#define _PHYSDAT_H

#include "astromat.h"
#include "enums.h"

double elong_helio( const vektor& p_h, const vektor& s_g );
double elongation( const vektor& p, const vektor& s );
double phasenwin( const vektor& p_h, const vektor& p_g );
double poswin_sonne( const vektor& p_g, const vektor& s_g );
double phase( double i );
double beldef( double i, double d );
double poswin_achse( const vektor& p_g, const vektor& a );
double zm_zb( const vektor& p, const vektor& a, double w, short s, double* zb );
double hor_par( double r );

short rot_sinn( PLANET planet );
void libration( double jde, vektor m_ekl, double* l, double* b );
double durchmesser( PLANET planet, double ep );
double helligkeit( PLANET planet, double sp, double ep, double i, double b );
void iau1982_j2000( PLANET planet, double jd,
                    double* alpha0, double* delta0, double* W );
void iau1988_j2000( PLANET planet, double jd,
		    double* alpha0, double* delta0, double* W );
void iau1991_j2000( PLANET planet, double jd,
		    double* alpha0, double* delta0, double* W );
void iau_datum( double jd, double alpha2000, double delta2000, double W2000,
                double* alpha, double* delta, double* W );
void rot_achse_null_mer( PLANET planet, double jd, vektor& a_equ, double* w );
double eins_minus_f_quadrat( PLANET planet );

#endif /* !defined (_PHYSDAT_H) */
