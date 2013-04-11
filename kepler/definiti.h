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

#if !defined (_DEFINITI_H)
#define _DEFINITI_H

#include "enums.h"

#define MAX_ELEMENT_NAME 80

struct ZEIT { short h,m,s; };
struct DATUM { short j,m,t; };
struct LEISTE { char *wt; DATUM date; double et_ut,jd,l,b; short dh,h; ZEIT zone,welt,stern; };
struct TABELLE_ECKDATEN { double jd_start,jd_ende; short intervall; };
struct ELEMENT { AUS_EIN      komet; /* Komet oder Kleinplanet            */
                 char name[MAX_ELEMENT_NAME]; /* Name des Objekts */
		 AEQUINOKTIUM a_mo;  /* Aequinoktium der Bahnelemente     */
		 double       a_jd,  /* Aequinoktium der Bahnelemente     */
			      e,     /* Bahnexzentrizitaet                */
			      q,     /* Periheldistanz                    */
			      a,     /* grosse Bahnhalbachse              */
			      n,     /* mittlere taegliche Bewegung       */
			      t0,    /* Zeitpunkt des Periheldurchgangs   */
			      t,     /* Zeitpunkt fuer den M gilt         */
			      M,     /* Mittlere Anomalie zum Zeitpunkt t */
			      om,    /* Argument des Perihels             */
			      pi,    /* Laenge des Perihels               */
		              kl,    /* Knotenlaenge (aufsteigender)      */
			      i,     /* Bahnneigung                       */
			      mag0,  /* 1. Helligkeitsparameter           */
			      mag1;  /* 2. Helligkeitsparameter           */
	       };

struct EREIGNIS { double jd; char text[80]; };

struct FINSTERNIS { double jde,groesse,gamma,dpp,dup,dut; FINSTER_ART art; PLANET p; };

#endif /* !defined (_DEFINITIO_H) */
