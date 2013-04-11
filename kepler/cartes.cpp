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
/* Module: CARTES.CPP                                                        */
/* Version 1.0                                                               */
/* Last modified: March 15, 1993                                             */
/*****************************************************************************/

#include <math.h>
#include <stdio.h>
#include "cartes.h"

/***********************/
/* siehe auch CARTES.H */
/***********************/

/* inline */ vektor::vektor()
{
   a[0] = 0.0;
   a[1] = 0.0;
   a[2] = 0.0;
}

/* inline */ vektor::vektor( const vektor& u )
{
   a[0] = u.a[0];
   a[1] = u.a[1];
   a[2] = u.a[2];
}

/* inline */ vektor::vektor( double x, double y, double z )
{
   a[0] = x;
   a[1] = y;
   a[2] = z;
}

vektor::vektor( const sphaer& s )
{
   double q = s.r * cos(s.beta);

   a[0] = q * cos(s.lambda);
   a[1] = q * sin(s.lambda);
   a[2] = s.r * sin(s.beta);
}

/* inline */ vektor::~vektor() {}

/* inline */ vektor& vektor::operator =  ( const vektor& u )
{
   a[0] = u.a[0];
   a[1] = u.a[1];
   a[2] = u.a[2];
   return *this;
}

/* inline */ vektor& vektor::operator *= ( double t )
{
   a[0] *= t;
   a[1] *= t;
   a[2] *= t;
   return *this;
}

/* inline */ vektor& vektor::operator /= ( double t )
{
   a[0] /= t;
   a[1] /= t;
   a[2] /= t;
   return *this;
}

/* inline */ vektor& vektor::operator += ( const vektor& u )
{
   a[0] += u.a[0];
   a[1] += u.a[1];
   a[2] += u.a[2];
   return *this;
}

/* inline */ vektor& vektor::operator -= ( const vektor& u )
{
   a[0] -= u.a[0];
   a[1] -= u.a[1];
   a[2] -= u.a[2];
   return *this;
}

/* inline */ vektor operator - ( const vektor& u )
{
   return vektor( -u.a[0],
                  -u.a[1],
                  -u.a[2] );
}

/* inline */ vektor operator * ( const vektor& u, double t )
{
   return vektor( u.a[0] * t,
                  u.a[1] * t,
                  u.a[2] * t );
}

/* inline */ vektor  operator * (double t, const vektor& u )
{
   return vektor( t * u.a[0],
                  t * u.a[1],
                  t * u.a[2] );
}

/* inline */ vektor operator * ( const matrix& A, const vektor& u )
{
   return vektor( A.m[0][0] * u.a[0] + A.m[0][1] * u.a[1] + A.m[0][2] * u.a[2],
                  A.m[1][0] * u.a[0] + A.m[1][1] * u.a[1] + A.m[1][2] * u.a[2],
                  A.m[2][0] * u.a[0] + A.m[2][1] * u.a[1] + A.m[2][2] * u.a[2]
                );
}

/* inline */ vektor operator / ( const vektor& u, double t )
{
   return vektor( u.a[0] / t,
                  u.a[1] / t,
                  u.a[2] / t );
}

/* inline */ vektor operator + ( const vektor& u, const vektor& v )
{
   return vektor( u.a[0] + v.a[0],
                  u.a[1] + v.a[1],
                  u.a[2] + v.a[2] );
}

/* inline */ vektor operator - ( const vektor& u, const vektor& v )
{
   return vektor( u.a[0] - v.a[0],
                  u.a[1] - v.a[1],
                  u.a[2] - v.a[2] );
}

/* inline */ vektor operator ^ ( const vektor &u, const vektor& v )
{
   return vektor( u.a[1]*v.a[2]-u.a[2]*v.a[1],
		  u.a[2]*v.a[0]-u.a[0]*v.a[2],
		  u.a[0]*v.a[1]-u.a[1]*v.a[0] );
}

/* inline */ double operator | ( const vektor& u, const vektor& v )
{
   return ( u.a[0]*v.a[0]+u.a[1]*v.a[1]+u.a[2]*v.a[2] );
}

void vektor::write( void )
{  short i;

   for(i=0;i<=2;i++) printf("%+12.5f",a[i]);
   printf("\n");
}

/* inline */ matrix::matrix()
{
   m[0][0] = 1.0; m[0][1] = 0.0; m[0][2] = 0.0;
   m[1][0] = 0.0; m[1][1] = 1.0; m[1][2] = 0.0;
   m[2][0] = 0.0; m[2][1] = 0.0; m[2][2] = 1.0;
}

/* inline */ matrix::matrix( const matrix& A )
{
    m[0][0] = A.m[0][0]; m[0][1] = A.m[0][1]; m[0][2] = A.m[0][2];
    m[1][0] = A.m[1][0]; m[1][1] = A.m[1][1]; m[1][2] = A.m[1][2];
    m[2][0] = A.m[2][0]; m[2][1] = A.m[2][1]; m[2][2] = A.m[2][2];
}

/* inline */ matrix::matrix( double a00, double a01, double a02,
	               double a10, double a11, double a12,
		       double a20, double a21, double a22 )
{
   m[0][0] = a00; m[0][1] = a01; m[0][2] = a02;
   m[1][0] = a10; m[1][1] = a11; m[1][2] = a12;
   m[2][0] = a20; m[2][1] = a21; m[2][2] = a22;
}

/* inline */ matrix::~matrix() {}

/* inline */ matrix& matrix::operator =  ( const matrix& A )
{
   m[0][0] = A.m[0][0]; m[0][1] = A.m[0][1]; m[0][2] = A.m[0][2];
   m[1][0] = A.m[1][0]; m[1][1] = A.m[1][1]; m[1][2] = A.m[1][2];
   m[2][0] = A.m[2][0]; m[2][1] = A.m[2][1]; m[2][2] = A.m[2][2];
   return *this;
}

/* inline */ matrix& matrix::operator *= ( double t )
{
   m[0][0] *= t; m[0][1] *= t; m[0][2] *= t;
   m[1][0] *= t; m[1][1] *= t; m[1][2] *= t;
   m[2][0] *= t; m[2][1] *= t; m[2][2] *= t;
   return *this;
}

/* inline */ matrix& matrix::operator /= ( double t )
{
   m[0][0] /= t; m[0][1] /= t; m[0][2] /= t;
   m[1][0] /= t; m[1][1] /= t; m[1][2] /= t;
   m[2][0] /= t; m[2][1] /= t; m[2][2] /= t;
   return *this;
}

/* inline */ matrix& matrix::operator += ( const matrix& A )
{
   m[0][0] += A.m[0][0]; m[0][1] += A.m[0][1]; m[0][2] += A.m[0][2];
   m[1][0] += A.m[1][0]; m[1][1] += A.m[1][1]; m[1][2] += A.m[1][2];
   m[2][0] += A.m[2][0]; m[2][1] += A.m[2][1]; m[2][2] += A.m[2][2];
   return *this;
}

/* inline */ matrix& matrix::operator -= ( const matrix& A )
{  m[0][0] -= A.m[0][0]; m[0][1] -= A.m[0][1]; m[0][2] -= A.m[0][2];
   m[1][0] -= A.m[1][0]; m[1][1] -= A.m[1][1]; m[1][2] -= A.m[1][2];
   m[2][0] -= A.m[2][0]; m[2][1] -= A.m[2][1]; m[2][2] -= A.m[2][2];
   return *this;
}

/* inline */ matrix operator - ( const matrix& A )
{
   return matrix( -A.m[0][0], -A.m[0][1], -A.m[0][2],
                  -A.m[1][0], -A.m[1][1], -A.m[1][2],
                  -A.m[2][0], -A.m[2][1], -A.m[2][2] );
}

/* inline */ matrix operator * ( const matrix& A, const matrix& B )
{
   return matrix( A.m[0][0]*B.m[0][0]+A.m[0][1]*B.m[1][0]+A.m[0][2]*B.m[2][0],
                  A.m[0][0]*B.m[0][1]+A.m[0][1]*B.m[1][1]+A.m[0][2]*B.m[2][1],
                  A.m[0][0]*B.m[0][2]+A.m[0][1]*B.m[1][2]+A.m[0][2]*B.m[2][2],

                  A.m[1][0]*B.m[0][0]+A.m[1][1]*B.m[1][0]+A.m[1][2]*B.m[2][0],
                  A.m[1][0]*B.m[0][1]+A.m[1][1]*B.m[1][1]+A.m[1][2]*B.m[2][1],
                  A.m[1][0]*B.m[0][2]+A.m[1][1]*B.m[1][2]+A.m[1][2]*B.m[2][2],

                  A.m[2][0]*B.m[0][0]+A.m[2][1]*B.m[1][0]+A.m[2][2]*B.m[2][0],
                  A.m[2][0]*B.m[0][1]+A.m[2][1]*B.m[1][1]+A.m[2][2]*B.m[2][1],
                  A.m[2][0]*B.m[0][2]+A.m[2][1]*B.m[1][2]+A.m[2][2]*B.m[2][2]
                );
}

/* inline */ matrix operator * ( const matrix& A, double t )
{
   return matrix( A.m[0][0] * t, A.m[0][1] * t, A.m[0][2] * t,
                  A.m[1][0] * t, A.m[1][1] * t, A.m[1][2] * t,
                  A.m[2][0] * t, A.m[2][1] * t, A.m[2][2] * t );
}

/* inline */ matrix operator * ( double t, const matrix& A )
{
   return matrix( t * A.m[0][0], t * A.m[0][1], t * A.m[0][2],
                  t * A.m[1][0], t * A.m[1][1], t * A.m[1][2],
                  t * A.m[2][0], t * A.m[2][1], t * A.m[2][2] );
}

/* inline */ matrix operator / ( const matrix& A, double t )
{
   return matrix( A.m[0][0] / t, A.m[0][1] / t, A.m[0][2] / t,
                  A.m[1][0] / t, A.m[1][1] / t, A.m[1][2] / t,
                  A.m[2][0] / t, A.m[2][1] / t, A.m[2][2] / t );
}

/* inline */ matrix operator + ( const matrix& A, const matrix& B )
{
   return matrix(
      A.m[0][0] + B.m[0][0], A.m[0][1] + B.m[0][1], A.m[0][2] + B.m[0][2],
      A.m[1][0] + B.m[1][0], A.m[1][1] + B.m[1][1], A.m[1][2] + B.m[1][2],
      A.m[2][0] + B.m[2][0], A.m[2][1] + B.m[2][1], A.m[2][2] + B.m[2][2]
                );
}

/* inline */ matrix operator - ( const matrix& A, const matrix& B )
{
   return matrix(
      A.m[0][0] - B.m[0][0], A.m[0][1] - B.m[0][1], A.m[0][2] - B.m[0][2],
      A.m[1][0] - B.m[1][0], A.m[1][1] - B.m[1][1], A.m[1][2] - B.m[1][2],
      A.m[2][0] - B.m[2][0], A.m[2][1] - B.m[2][1], A.m[2][2] - B.m[2][2]
                );
}

/* inline */ double det( const matrix& A )
{
   return ( +A.m[0][0] * A.m[1][1] * A.m[2][2]
            -A.m[0][0] * A.m[1][2] * A.m[2][1]
            +A.m[0][2] * A.m[1][0] * A.m[2][1]
            -A.m[0][1] * A.m[1][0] * A.m[2][2]
            +A.m[0][1] * A.m[1][2] * A.m[2][0]
            -A.m[0][2] * A.m[1][1] * A.m[2][0]
          );
}

/* inline */ matrix transponiert( const matrix& A )
{
   return matrix( A.m[0][0], A.m[1][0], A.m[2][0],
                  A.m[0][1], A.m[1][1], A.m[2][1],
                  A.m[0][2], A.m[1][2], A.m[2][2] );

}

matrix invertiert( const matrix& A )
{
   double d;

   d=1.0/det(A);

   return matrix( +d*(A.m[1][1]*A.m[2][2]-A.m[1][2]*A.m[2][1]),
                  -d*(A.m[0][1]*A.m[2][2]-A.m[0][2]*A.m[2][1]),
                  +d*(A.m[0][1]*A.m[1][2]-A.m[0][2]*A.m[1][1]),

                  -d*(A.m[1][0]*A.m[2][2]-A.m[1][2]*A.m[2][0]),
                  +d*(A.m[0][0]*A.m[2][2]-A.m[0][2]*A.m[2][0]),
                  -d*(A.m[0][0]*A.m[1][2]-A.m[0][2]*A.m[1][0]),

                  +d*(A.m[1][0]*A.m[2][1]-A.m[1][1]*A.m[2][0]),
                  -d*(A.m[0][0]*A.m[2][1]-A.m[0][1]*A.m[2][0]),
                  +d*(A.m[0][0]*A.m[1][1]-A.m[0][1]*A.m[1][0])
                );
}

void matrix::write( void )
{
   short i,k;

   for(i=0;i<=2;i++)
   {
      for(k=0;k<=2;k++)
         printf("%+12.5f",m[i][k]);
      printf("\n");
   }
}

/* inline */ sphaer::sphaer()
{
   lambda=0.0;
   beta=0.0;
   r=1.0;
}

/* inline */ sphaer::sphaer( const sphaer& s )
{
   lambda=s.lambda;
   beta=s.beta;
   r=s.r;
}

/* inline */ sphaer::sphaer( double LAMBDA, double BETA, double R )
{
   lambda=LAMBDA;
   beta  =BETA;
   r     =R;
}

/* inline */ sphaer::~sphaer() {}

/* inline */ sphaer& sphaer::operator =  ( const sphaer& s )
{
   lambda=s.lambda;
   beta  =s.beta;
   r     =s.r;
   return *this;
}

sphaer::sphaer( const vektor& u )
{  double q=u.a[0]*u.a[0]+u.a[1]*u.a[1];

   if (q==0.0) lambda=0.0;
   else        lambda=atan2(u.a[1],u.a[0]);

   r=sqrt(q+u.a[2]*u.a[2]);
   q=sqrt(q);

   if (r==0.0) beta=0.0;
   else        beta=atan2(u.a[2],q);
}

void sphaer::write( void )
{
   printf("%+12.5f%+12.5f%+12.5f\n",lambda,beta,r);
}
