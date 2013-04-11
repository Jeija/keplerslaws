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

#if !defined (_CARTES_H)
#define _CARTES_H

class vektor;
class matrix;
class sphaer;

class vektor {
friend sphaer;
friend matrix;
private:
   double a[3];
public:
   vektor();
   vektor( const vektor& u );
   vektor( double x, double y, double z );
   vektor( const sphaer& s );
   ~vektor();

   vektor& operator  = ( const vektor& u );
   vektor& operator *= ( double t );
   vektor& operator /= ( double t );
   vektor& operator += ( const vektor& u );
   vektor& operator -= ( const vektor& u );
   inline double  operator [] ( short i ) const {  return a[i];  }
   inline double& operator [] ( short i ) {  return a[i];  }

   friend vektor  operator -  ( const vektor& u );
   friend vektor  operator *  ( const vektor& u, double t );
   friend vektor  operator *  ( double t, const vektor& u );
   friend vektor  operator *  ( const matrix& A, const vektor& u );
   friend vektor  operator /  ( const vektor& u, double t );
   friend vektor  operator +  ( const vektor& u, const vektor& v );
   friend vektor  operator -  ( const vektor& u, const vektor& v );
   friend vektor  operator ^  ( const vektor& u, const vektor& v );
   friend double  operator |  ( const vektor& u, const vektor& v );

   void write( void );
};

class matrix {
friend vektor;
private:
   double m[3][3];
public:
   matrix();
   matrix( const matrix& A );
   matrix( double a00, double a01, double a02,
	   double a10, double a11, double a12,
	   double a20, double a21, double a22 );
   ~matrix();

   matrix& operator  = ( const matrix& A );
   matrix& operator *= ( double t );
   matrix& operator /= ( double t );
   matrix& operator += ( const matrix& A );
   matrix& operator -= ( const matrix& A );
   inline double* operator [] ( short i ) {  return m[i];  }

   friend matrix  operator -  ( const matrix& A );
   friend matrix  operator *  ( const matrix& A, const matrix& B );
   friend vektor  operator *  ( const matrix& A, const vektor& u );
   friend matrix  operator *  ( const matrix& A, double t );
   friend matrix  operator *  ( double t, const matrix& A );
   friend matrix  operator /  ( const matrix& A, double t );
   friend matrix  operator +  ( const matrix& A, const matrix& B );
   friend matrix  operator -  ( const matrix& A, const matrix& B );

   friend double det          ( const matrix& A );
   friend matrix transponiert ( const matrix& A );
   friend matrix invertiert   ( const matrix& A );

   void write( void );
};

class sphaer {
public:
   double lambda,beta,r;

   sphaer();
   sphaer( const sphaer& s );
   sphaer( double LAMBDA, double BETA, double R );
   sphaer( const vektor& u );
   ~sphaer();

   sphaer& operator  = ( const sphaer& s );

   void write( void );
};

#endif /* !defined (_CARTES_H) */
