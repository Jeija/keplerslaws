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

#if !defined (_TSCHEBY_H)
#define _TSCHEBY_H

#include <stdio.h>

class tschebyscheff {
private:
   int Grad;                          // Grad des TP
   double *L_coeff,*B_coeff,*R_coeff; // Koeffizienten
   double b_minus_a_2,a_plus_b_2;     // Intervall [a,b] mit
				      // a = (a_plus_b_2-b_minus_a_2)/2
				      // b = (a_plus_b_2+b_minus_a_2)/2
public:
   tschebyscheff( void )
   {
      Grad=-1;
      L_coeff=NULL;
      B_coeff=NULL;
      R_coeff=NULL;
   }
   tschebyscheff( double a, double b, int n, sphaer (*f)(double) );
   tschebyscheff( PLANET planet, double a, double b, int n, sphaer (*f)(PLANET,double) );
   ~tschebyscheff();

   double lambda( double x ) const;
   double beta  ( double x ) const;
   double r     ( double x ) const;
   sphaer pos   ( double t ) const;
   void   tschebyscheffpoly( double a, double b, int n, sphaer (*f)(double) );
   void   tschebyscheffpoly( PLANET planet, double a, double b, int n, sphaer (*f)(PLANET,double) );
};

#endif /* !defined (_TSCHEBY_H) */

