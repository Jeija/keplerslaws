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
/*
   The planetary terms (distributed in PLANETS.NDX and PLANETS.DAT) are
   not copyrighted.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "astromat.h"
#include "definiti.h"
#include "vsop87.h"

/* geometric heliocentric position of planet, mean ecliptic of date
 * (not corrected for light-time)
 */
sphaer planet_pos_xephem( PLANET planet, double jde )
{
   double mjd=jde-2415020.0;
   double prec=0.0;
   double ret[6];
   
   sphaer S;
   
   vsop87(mjd, int(planet), prec, ret);
   S.lambda= ret[0];
   S.beta  = ret[1];
   S.r     = ret[2];
   
   return S;
}
