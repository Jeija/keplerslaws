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

#if !defined (_POSITION_H)
#define _POSITION_H

#include "astromat.h"
#include "enums.h"
#include "ppos_lo.h"

vektor element_gekl_geom( const ELEMENT& e, double jde, GENAUIGKEIT g );
vektor element_gequ_app( const ELEMENT& e, double jde, GENAUIGKEIT g  );
void element_pos( const ELEMENT& e,
                  double jde,
                  KOORDINATEN_TYP koor_typ,
                  AEQUINOKTIUM aqu_mo,
                  double aqu_jd,
                  GENAUIGKEIT g,

	          vektor& k_hekl_g, vektor& k_gekl_g, vektor& k_gequ_g,
                  vektor& s_gekl_g, vektor& k_gekl_k, vektor& k_gequ_k,
	          double* delta, double* tau, matrix& matrix_ekl2equ
                );
double element_helligkeit( const ELEMENT& element, double i, double sp, double ep );
sphaer planet_lo( PLANET planet, double jde );
sphaer mond_lo( double jde );
sphaer planet_hi( PLANET planet, double jde );
sphaer mond_hi( double jde );
vektor planet_gekl_geom( PLANET planet, double jde, GENAUIGKEIT g  );
vektor planet_gequ_app( PLANET planet, double jde, GENAUIGKEIT g  );
void planet_pos( PLANET planet,
                 double jde,
                 KOORDINATEN_TYP korr_typ,
                 AEQUINOKTIUM aqu_mo,
                 double aqu_jd,
                 GENAUIGKEIT g,

	         vektor& p_hekl_g, vektor& p_gekl_g, vektor& p_gequ_g,
                 vektor& s_gekl_g, vektor& p_gekl_k, vektor& p_gequ_k,
	         double* delta, double* tau, matrix& matrix_ekl2equ
               );

#endif /* !defined (_POSITION_H) */
