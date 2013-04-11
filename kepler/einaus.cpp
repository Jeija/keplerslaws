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
/* Module: EINAUS.CPP                                                        */
/* Version 1.0                                                               */
/* Last modified: February 22, 1993                                          */
/*****************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "astromat.h"
#include "einaus.h"

/************************************/
/* Elementfunktionen der Klasse HMS */
/************************************/
/* HMS()                            */
/* HMS(double)                      */
/* HMS.format(WINKEL)               */
/* HMS.format(WINKEL,short)           */
/************************************/

HMS::HMS()
/* Default-Initialisierung */
{  i_h=i_m=i_s=-9999;
   d_h=d_m=d_s=-9999.0;
}

HMS::HMS(double t)
/* Initialisierung */
{  h=mod(t,24.0);
}

void HMS::format(WINKEL w)
/* Formatierung auf ganze Stunden oder Minuten oder Sekunden */
{  ldiv_t X;
   i_h=i_m=i_s=-9999;
   d_h=d_m=d_s=-9999.0;
   switch(w)
   {  case HH:X=ldiv((long)(h+0.5),24L);
	      i_h=(short)X.rem;
	      c=(short)X.quot;
	      break;                        // auf i-stunden
      case MM:X=ldiv((long)(h*60.0+0.5),60L);
	      i_m=(short)X.rem;
	      X=ldiv(X.quot,24L);
	      i_h=(short)X.rem;
	      c=(short)X.quot;
	      break;                        // auf i-minuten
      case SS:X=ldiv((long)(h*3600.0+0.5),60L);
	      i_s=(short)X.rem;
	      X=ldiv(X.quot,60L);
	      i_m=(short)X.rem;
	      X=ldiv(X.quot,24L);
	      i_h=(short)X.rem;
	      c=(short)X.quot;		    // auf i-sekunden
              break;
      default:assert(0); break;
   }
}

void HMS::format(WINKEL w, short d)
/* Formatierung auf Stunden- oder Minuten- oder Sekundenbruchteile */
{  ldiv_t X;
   double dez=pow(10.0,d),u,v;
   i_h=i_m=i_s=-9999;
   d_h=d_m=d_s=-9999.0;
   switch(w)
   {  case HH:u=floor(h*dez+0.5)/dez;
	      d_h=fmod(u,24.0);
	      c=(short)(u/24.0);
	      break;                        // auf d-stunden
      case MM:u=floor(h*60.0*dez+0.5)/dez;
	      d_m=fmod(u,60.0);
	      v=floor(u/60.0);
	      X=ldiv((long)v,24L);
	      i_h=(short)X.rem;
	      c=(short)X.quot;
	      break;                        // auf d-minuten
      case SS:u=floor(h*3600.0*dez+0.5)/dez;
	      d_s=fmod(u,60.0);
	      v=floor(u/60.0);
	      X=ldiv((long)v,60L);
	      i_m=(short)X.rem;
	      X=ldiv(X.quot,24L);
	      i_h=(short)X.rem;
	      c=(short)X.quot;		    // auf d-sekunden
              break;
      default:assert(0); break;
   }
}

/************************************/
/* Elementfunktionen der Klasse DMS */
/************************************/
/* DMS()                            */
/* DMS(double)                      */
/* DMS.format(WINKEL)               */
/* DMS.format(WINKEL,short)           */
/* DMS.format(WINKEL,short,short)       */
/* DMS.format(WINKEL,short,short,short)   */
/* DMS.vd2s(char)                   */
/************************************/

DMS::DMS()
/* Default-Initialisierung */
{  i_d=i_m=i_s=-9999;
   d_d=d_m=d_s=-9999.0;
}

DMS::DMS(double x)
/* Initialisierung */
{  d=x;
}

void DMS::format(WINKEL w)
/* Formatierung auf ganze Grad oder Minuten oder Sekunden */
{  ldiv_t X;
   long da,u;

   i_d=i_m=i_s=-9999;
   d_d=d_m=d_s=-9999.0;
   switch(w)
   {  case DD:u=(long)floor(d+0.5);
              if (u<0) { da=-u; v=-1; } else { da=u; v=+1; }
//	      v= (u<0) ? (da=-u), -1 : (da=u), +1;
	      i_d=(short)da;
	      break;                        // auf i-grad
      case MM:u=(long)floor(d*60.0+0.5);
              if (u<0) { da=-u; v=-1; } else { da=u; v=+1; }
//	      v= (u<0) ? (da=-u), -1 : (da=u), +1;
	      X=ldiv(da,60L);
	      i_m=(short)X.rem;
	      i_d=(short)X.quot;
	      break;                        // auf i-minuten
      case SS:u=(long)floor(d*3600.0+0.5);
              if (u<0) { da=-u; v=-1; } else { da=u; v=+1; }
//	      v= (u<0) ? (da=-u), -1 : (da=u), +1;
	      X=ldiv(da,60L);
	      i_s=(short)X.rem;
	      X=ldiv(X.quot,60L);
	      i_m=(short)X.rem;
	      i_d=(short)X.quot;
	      break;              	    // auf i-sekunden
      default:assert(0); break;
   }
}

void DMS::format(WINKEL w, short p)
/* Formatierung auf Grad- oder Minuten- oder Sekundenbruchteile */
{  ldiv_t X;
   double da,u;
   double dez=pow(10.0,p);

   i_d=i_m=i_s=-9999;
   d_d=d_m=d_s=-9999.0;
   switch(w)
   {  case DD:u=(long)floor(d*dez+0.5)/dez;
              if (u<0) { da=-u; v=-1; } else { da=u; v=+1; }
//	      v= (u<0) ? (da=-u), -1 : (da=u), +1;
	      d_d=da;
	      break;                        // auf d-grad
      case MM:u=(long)floor(d*dez*60.0+0.5)/dez;
              if (u<0) { da=-u; v=-1; } else { da=u; v=+1; }
//	      v= (u<0) ? (da=-u), -1 : (da=u), +1;
	      d_m=fmod(da,60.0);
	      i_d=(short)floor(da/60.0);
	      break;                        // auf d-minuten
      case SS:u=(long)floor(d*dez*3600.0+0.5)/dez;
              if (u<0) { da=-u; v=-1; } else { da=u; v=+1; }
//	      v= (u<0) ? (da=-u), -1 : (da=u), +1;
	      d_s=fmod(da,60.0);
	      X=ldiv((long)floor(da/60.0),60L);
	      i_m=(short)X.rem;
	      i_d=(short)X.quot;
	      break;              	    // auf i-sekunden
      default:assert(0); break;
   }
}

void DMS::format(WINKEL w, short min, short max)
/* Formatierung auf ganze Grad oder Minuten oder Sekunden
   mit Modulo-Funktion (ACHTUNG: es gilt min <= d <  max)
*/
{  ldiv_t X;
   long da,u;

   max-=min;
   i_d=i_m=i_s=-9999;
   d_d=d_m=d_s=-9999.0;
   switch(w)
   {  case DD:u=(long)(mod(floor(d-min+0.5),max)+(double)min);
              if (u<0) { da=-u; v=-1; } else { da=u; v=+1; }
//	      v= (u<0) ? (da=-u), -1 : (da=u), +1;
	      i_d=(short)da;
	      break;                        // auf i-grad
      case MM:u=(long)(mod(floor((d-min)*60.0+0.5),max*60.0)+min*60.0);
              if (u<0) { da=-u; v=-1; } else { da=u; v=+1; }
//	      v= (u<0) ? (da=-u), -1 : (da=u), +1;
	      X=ldiv(da,60L);
	      i_m=(short)X.rem;
	      i_d=(short)X.quot;
	      break;                        // auf i-minuten
      case SS:u=(long)(mod(floor((d-min)*3600.0+0.5),3600.0*(double)max)
	       +3600.0*double(min));
              if (u<0) { da=-u; v=-1; } else { da=u; v=+1; }
//	      v= (u<0) ? (da=-u), -1 : (da=u), +1;
	      X=ldiv(da,60L);
	      i_s=(short)X.rem;
	      X=ldiv(X.quot,60L);
	      i_m=(short)X.rem;
	      i_d=(short)X.quot;
	      break;              	    // auf i-sekunden
      default:assert(0); break;
   }
}

void DMS::format(WINKEL w, short p, short min, short max)
/* Formatierung auf Grad- oder Minuten- oder Sekundenbruchteile
   mit Modulo-Funktion (ACHTUNG: es gilt fuer d [min,max[)
*/
{  ldiv_t X;
   double da,u;
   double dez=pow(10.0,p);

   i_d=i_m=i_s=-9999;
   d_d=d_m=d_s=-9999.0;
   max-=min;
   switch(w)
   {  case DD:u=floor(d*dez+0.5)/dez;
	      u=mod(u-min,max)+(double)min;
              if (u<0) { da=-u; v=-1; } else { da=u; v=+1; }
//	      v= (u<0) ? (da=-u), -1 : (da=u), +1;
	      d_d=da;
	      break;                        // auf d-grad
      case MM:u=floor(d*dez*60.0+0.5)/dez;
	      u=mod(u-60.0*min,60.0*max)+60.0*min;
              if (u<0) { da=-u; v=-1; } else { da=u; v=+1; }
//	      v= (u<0) ? (da=-u), -1 : (da=u), +1;
	      d_m=fmod(da,60.0);
	      i_d=(short)floor(da/60.0);
	      break;                        // auf d-minuten
      case SS:u=floor(d*dez*3600.0+0.5)/dez;
	      u=mod(u-3600.0*min,3600.0*max)+3600.0*min;
              if (u<0) { da=-u; v=-1; } else { da=u; v=+1; }
//	      v= (u<0) ? (da=-u), -1 : (da=u), +1;
	      d_s=fmod(da,60.0);
	      X=ldiv((long)floor(da/60.0),60L);
	      i_m=(short)X.rem;
	      i_d=(short)X.quot;
	      break;              	    // auf i-sekunden
      default:assert(0); break;
   }
}

char* DMS::vd2s(char* s)
/* Bildet aus i_d und v einen String */
{  if ((i_d==0)&&(v<0))
   {  sprintf(s,"-0");
      return s;
   }
   else
   {  sprintf(s,"%d",v*i_d);
      return s;
   }
}

/**************************/
/* Gradmass <-> Bogenmass */
/**************************/
/* deg()                  */
/* rad()                  */
/**************************/

double deg( double radian )
/* Bogenmass -> Gradmass [0,360[ */
{
   return mod(radian*RAD2DEG,360.0);
}

double rad( double degree )
/* Gradmass -> Bogenmass [0,2*Pi[ */
{
   return mod(degree*DEG2RAD,M_2PI);
}
