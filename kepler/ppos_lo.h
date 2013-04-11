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

#if !defined (_PPOS_LO_H)
#define _PPOS_LO_H

#include "astromat.h"

sphaer son( double jd );
sphaer mer( double jd );
sphaer ven( double jd );
sphaer mon( double jd );
sphaer mar( double jd );
sphaer jup( double jd );
sphaer sat( double jd );
sphaer ura( double jd );
sphaer nep( double jd );
sphaer plu( double jd );

#endif /* !defined (_PPOS_LO_H) */
