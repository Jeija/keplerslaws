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

#if !defined (_EINAUS_H)
#define _EINAUS_H

#include "enums.h"

class HMS;
class DMS;

class HMS {
public: double h;
	short    i_h;    // i-stunden
	short    i_m;    // i-minuten
	short    i_s;    // i-sekunden
	double d_h;
	double d_m;
	double d_s;
	short    c;      // uebertrag fuer 24h

   HMS();
   HMS(double);
   void format(WINKEL);
   void format(WINKEL,short);
};

class DMS {
public: double d;
	short    i_d;   // i-grad
	short    i_m;   // i-minuten
	short    i_s;   // i-sekunden
	double d_d;
	double d_m;
	double d_s;
	short    v;     // Vorzeichen

   DMS();
   DMS(double);
   void format(WINKEL);
   void format(WINKEL,short);
   void format(WINKEL,short,short);
   void format(WINKEL,short,short,short);
   char* vd2s(char* S);
};

double deg( double rad );
double rad( double deg );

#endif /* !defined (_EINAUS_H) */
