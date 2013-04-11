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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "astromat.h"
#include "definiti.h"
#include "position.h"
#include "tscheby.h"
#include "tp_pos.h"

// #define TEST_TSCHEBY
// sphaer planet_pos_hi( short planet, double jde );

static void error( const char *s ) { fprintf(stderr,"%s\n",s); }

tp_position::tp_position( double jde_s, double jde_e, PLANET p1, PLANET p2, sphaer (*f)(PLANET,double) )
{
   short i,p;
			  /*     so    me    ve    er    ma    ju    sa    ur    ne    pl    mo */
   short tp_grad[11]       = {   25,   25,   15,   25,   10,    5,    5,    3,    3,    3,   10};
   double tp_intervall[11] = {200.0, 67.0,200.0,200.0,400.0,400.0,400.0,400.0,400.0,400.0, 20.0};


   for (i=0;i<11;i++)
   {
      TP_GRAD[i]=tp_grad[i];
      TP_INTERVALL[i]=tp_intervall[i];
   }

   jde_start=jde_s;
   jde_ende =jde_e;

   p_erster =p1;
   p_letzter=p2;

   /* TP ausrechnen fuer jde_s bis jde_e */

   for(p=short(p_erster);p<=short(p_letzter);p++)
   {
      short i;

      max_index[p]=short(ceil((jde_ende-jde_start)/TP_INTERVALL[p])); // Anzahl der Teilstuecke

      TP_PLANET[p] = new tschebyscheff[max_index[p]];

      if (!TP_PLANET[p])
      {
	 error("\nKein Speicher frei fuer Tschebyscheff-Polynom-Array!");
	 exit(3);
      }

      for(i=0;i<max_index[p];i++)
      {
	 TP_PLANET[p][i].tschebyscheffpoly(PLANET(p),jde_start+i*TP_INTERVALL[p]-1.0,jde_start+(i+1)*TP_INTERVALL[p]+1.0,TP_GRAD[p],f);
      }
   }
}

tp_position::~tp_position()
{
   short p;

   for(p=short(p_erster);p<=short(p_letzter);p++)
   {
      delete[] TP_PLANET[p];
   }
}

sphaer tp_position::pos( PLANET planet, double jde ) const
{
   sphaer s;
   short index;

   index=short(floor((jde-jde_start)/TP_INTERVALL[short(planet)]));

   /* einige Sicherheitsueberpruefungen */

   if ( (index < 0) || (index >= max_index[short(planet)]) )
   {
      error("Ungueltiger INDEX!\n");
      exit(3);
   }

   s=TP_PLANET[short(planet)][index].pos(jde);

   return s;
}
