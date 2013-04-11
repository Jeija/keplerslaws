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

#ifndef _TP_POS_H
#define _TP_POS_H

class tp_position {
private:
   tschebyscheff *TP_PLANET[11]; // TP fuer alle Planeten + Mond
   int    TP_GRAD[11];      // Grad der TP im Intervall TP_INTERVALL
   double TP_INTERVALL[11]; // s.o.
   short  max_index[11];    // Anzahl der Aufteilungen
   PLANET p_erster;         // Ab diesem Planet werden TP berechnet
   PLANET p_letzter;        // Letzter Planet, fuer den TP berechnet werden

   double jde_start;        // Anfangsdatum
   double jde_ende;         // Enddatum
public:
   tp_position( double jde_s, double jde_e, PLANET p1, PLANET p2, sphaer (*f)(PLANET,double) );
   ~tp_position();

   sphaer pos( PLANET planet, double jde ) const;
};

#endif
