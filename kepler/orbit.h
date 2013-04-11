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

#if !defined (_ORBIT_H)
#define _ORBIT_H

void pos_ell( double a,   double exz, double m,
	      double* x,  double* y,
	      double* vx, double* vy );
void pos_par( double q, double d, double d0,
	      double* x,  double* y,
	      double* vx, double* vy );
void pos_hyp( double a,   double exz, double d, double d0,
	      double* x,  double* y,
	      double* vx, double* vy );
void gaussvek( double k, double i, double w, vektor& P, vektor& Q );
vektor bahn2ekl( double x, double y, const vektor& P, const vektor& Q );
vektor kepler( double e, double q, double d, double d0,
               const vektor& P, const vektor& Q,
               vektor& v_helio_ekl );

#endif /* !defined (_ORBIT_H) */
