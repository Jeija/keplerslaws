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

#ifndef _GV_COMM_H
#define _GV_COMM_H

#ifdef WIN32
#include "glut.h"
#else
#include <GL/glut.h>
#endif
#include "astromat.h"
#include "definiti.h"

#define MAX_PLANETS  10
#define MAX_ELEMENTS 10

struct PLANETEN {
   PLANET  p[MAX_PLANETS];
   int     on[MAX_PLANETS];
   int n;
};
struct ELEMENTE {
   ELEMENT e[MAX_ELEMENTS];
   int     on[MAX_ELEMENTS];
   int n;
};
struct VEGA_STERN {
   signed long lam,bet,mag;
   int nr;
};
enum KAMERA_MODUS { FIXOBJ, FOLLOWOBJ, KAMERA_MODUS_ENDE };
enum OBJEKT_MODUS { MESH, KUGEL, DISK, OBJEKT_MODUS_ENDE };

extern char      PLANET_NAME[][8];

void fourmatrix( matrix t, vektor s, double m[16] );
void gv_init( const PLANETEN& planet, const ELEMENTE& elem );
void gv_sonne( double jd );
void gv_erde ( const tschebyscheff& tps, double jd );
void gv_zentrallinie( const tschebyscheff& tps, double jd );
void gv_mond ( const tschebyscheff& tps, const tschebyscheff& tpm, double jd );
void gv_SonneErdeMond_init( void );
void gv_planets( const PLANETEN& planet,
                 const ELEMENTE& elem, 
                 double jd, double RADIUS, double AEQUINOKTIUM );
void gv_orbits( AUS_EIN a, const PLANETEN& planet, const ELEMENTE& elem, double RADIUS, double AEQUINOKTIUM );
void gv_apsiden( AUS_EIN a, const PLANETEN& planet, const ELEMENTE& elem, double RADIUS, double AEQUINOKTIUM );
void gv_focus( AUS_EIN a, const PLANETEN& planet, const ELEMENTE& elem, double RADIUS, double AEQUINOKTIUM );
void gv_camera( AUS_EIN reset, KAMERA_MODUS mode, vektor target, vektor position, vektor fixobj, vektor fixdir );
void gv_objekt( const char* objekt_name, double RADIUS, OBJEKT_MODUS mode,
                AUS_EIN label, double size, double lift );
void gv_sterne( AUS_EIN ausein, double JD, double AEQUINOKTIUM );
void gv_polygon( const char* name, int n, vektor v[] );
void gv_strahl( const char* name, vektor start, vektor end );
void gv_cone( const char *geom, const vektor& spitze, const vektor& basis, double basis_radius );

void og_init( const PLANETEN& planet, const ELEMENTE& elem );
void og_planets( const PLANETEN& planet,
                 const ELEMENTE& elem, 
                 double jd, double RADIUS, double AEQUINOKTIUM );
void og_orbit( const ELEMENT& elem, double AEQUINOKTIUM );
void og_apsiden( AUS_EIN a, const PLANETEN& planet, const ELEMENTE& elem, double RADIUS, double AEQUINOKTIUM );
void og_focus( const ELEMENT& elem, double AEQUINOKTIUM );
void og_camera( KAMERA_MODUS mode, vektor target, vektor position, vektor fixobj, vektor fixdir );
void og_sterne( AUS_EIN ausein, double JD, double AEQUINOKTIUM );
void og_polygon( int n, vektor v[] );
void og_strahl( vektor start, vektor end );
void og_erde( const vektor& v, double jd, int list, int n, vektor *ve );
void og_mond( const vektor& v, double jd, int list );
void og_kegel(const vektor& spitze, const vektor& basis, double basis_radius);

#endif
