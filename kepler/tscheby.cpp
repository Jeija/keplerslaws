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
/* TP-Klasse Version 1.1 vom 25. 3.1994 (C) Tobias Kramer */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "astromat.h"
#include "definiti.h"
#include "tscheby.h"

tschebyscheff::tschebyscheff( double a, double b, int n, sphaer (*f)(double) )
{
   tschebyscheffpoly(a,b,n,f);
}

tschebyscheff::tschebyscheff( PLANET p, double a, double b, int n, sphaer (*f)(PLANET,double) )
{
   tschebyscheffpoly(p,a,b,n,f);
}

tschebyscheff::~tschebyscheff()
{
   delete[] L_coeff; L_coeff = NULL;
   delete[] B_coeff; B_coeff = NULL;
   delete[] R_coeff; R_coeff = NULL;
}

double tschebyscheff::lambda( double x ) const
/*
    TP auswerten fuer lambda.
*/
{
   if (Grad<0)
   {
      fprintf(stderr,"\nFehler in '%s' Zeile %d: Grad < 0!\n",__FILE__,int(__LINE__));
      abort();
   }

   double y,u,v,w;
   int    k;

   u=0;
   v=0;
   y=2*(x-a_plus_b_2)/b_minus_a_2;
   for(k=Grad;k>=1;k--)
   {
      w=u;
      u=y*u-v+(*(L_coeff+k));
      v=w;
   }
   return(y/2.0*u-v+0.5*(*L_coeff));
}

double tschebyscheff::beta( double x ) const
/*
    TP auswerten fuer beta.
*/
{
   if (Grad<0)
   {
      fprintf(stderr,"\nFehler in '%s' Zeile %d: Grad < 0!\n",__FILE__,int(__LINE__));
      abort();
   }

   double y,u,v,w;
   int    k;

   u=0;
   v=0;
   y=2*(x-a_plus_b_2)/b_minus_a_2;
   for(k=Grad;k>=1;k--)
   {
      w=u;
      u=y*u-v+(*(B_coeff+k));
      v=w;
   }
   return(y/2.0*u-v+0.5*(*B_coeff));
}

double tschebyscheff::r( double x ) const
/*
    TP auswerten fuer r.
*/
{
   if (Grad<0)
   {
      fprintf(stderr,"\nFehler in '%s' Zeile %d: Grad < 0!\n",__FILE__,int(__LINE__));
      abort();
   }

   double y,u,v,w;
   int    k;

   u=0;
   v=0;
   y=2*(x-a_plus_b_2)/b_minus_a_2;
   for(k=Grad;k>=1;k--)
   {
      w=u;
      u=y*u-v+(*(R_coeff+k));
      v=w;
   }
   return(y/2.0*u-v+0.5*(*R_coeff));
}

sphaer tschebyscheff::pos( double x ) const
/*
    TP auswerten fuer lambda,beta und r.
*/
{
   if (Grad<0)
   {
      fprintf(stderr,"\nFehler in '%s' Zeile %d: Grad < 0!\n",__FILE__,int(__LINE__));
      abort();
   }

   double y;
   double ul,vl,wl;
   double ub,vb,wb;
   double ur,vr,wr;
   int    k;
   sphaer s;

   ul=ub=ur=0.0;
   vl=vb=vr=0.0;

   y=2*(x-a_plus_b_2)/b_minus_a_2;

   for(k=Grad;k>=1;k--)
   {
      wl=ul;
      wb=ub;
      wr=ur;
      ul=y*ul-vl+(*(L_coeff+k));
      ub=y*ub-vb+(*(B_coeff+k));
      ur=y*ur-vr+(*(R_coeff+k));
      vl=wl;
      vb=wb;
      vr=wr;
   }
   s.lambda=(y/2.0*ul-vl+0.5*(*L_coeff));
   s.beta  =(y/2.0*ub-vb+0.5*(*B_coeff));
   s.r     =(y/2.0*ur-vr+0.5*(*R_coeff));

   return s;
}

void tschebyscheff::tschebyscheffpoly( double a, double b, int n, sphaer (*f)(double) )
/*
    Diese Funktion berechnet die Koeffizienten des TP zur Funktion f
    vom Grad n im Intervall [a,b]. f wird dabei n+1 mal aufgerufen.
*/
{
   int k,j;
   double* h_werte = new double[ 2*n+2 ];
   sphaer* f_werte = new sphaer[ n+1 ];

   if (L_coeff) { delete[] L_coeff; L_coeff = NULL; }
   if (B_coeff) { delete[] B_coeff; B_coeff = NULL; }
   if (R_coeff) { delete[] R_coeff; R_coeff = NULL; }

   L_coeff = new double[ n+1 ];
   B_coeff = new double[ n+1 ];
   R_coeff = new double[ n+1 ];

   if ((h_werte==NULL)||(f_werte==NULL)||
       (L_coeff==NULL)||(B_coeff==NULL)||(R_coeff==NULL))
   {
      fprintf(stderr,"Kein Speicher fuer TP!");
      exit(1);
   }

   Grad=n;
   a_plus_b_2 =(a+b)*0.5;
   b_minus_a_2=(b-a)*0.5;
   *h_werte=1.0;                                            // Hilfstabelle fuer
   *(h_werte+1)=cos(PI/(2*n+2));                      // Stuetzstellen
   for (k=2;k<=2*n+1;k++)                             // x(k)=h_werte(2k+1)
      *(h_werte+k)=2*(*(h_werte+1))*(*(h_werte+k-1))-(*(h_werte+k-2));
   for (k=0;k<=n;k++)                                 // f an Stuetzstellen
   {
      *(f_werte+k)=(*f)(*(h_werte+2*k+1)*b_minus_a_2+a_plus_b_2); //  auswerten
   }
   for (k=1;k<=n;k++) // lambda stetig machen in [-2*Pi,+2*Pi]
   {
      if (((f_werte+k-1)->lambda)<((f_werte+k)->lambda))
	 ((f_werte+k)->lambda)-=M_2PI;
   }
   for (j=0;j<=n;j++)                                 // T-Koeffizienten
   {                                                  // berechenen
      *(h_werte+1)=cos(PI*j/(2*n+2));                     // Hilfstabelle fuer TP
      for (k=2;k<=2*n+1;k++)
	*(h_werte+k)=2*(*(h_werte+1))*(*(h_werte+k-1))-(*(h_werte+k-2));      // T(j,x(k))=h_werte(2k+1)
      *(L_coeff+j)=0;
      *(B_coeff+j)=0;
      *(R_coeff+j)=0;
      for (k=0;k<=n;k++)
      {
	 *(L_coeff+j)=(*(L_coeff+j))+(*(h_werte+2*k+1))*((f_werte+k)->lambda);
	 *(B_coeff+j)=(*(B_coeff+j))+(*(h_werte+2*k+1))*((f_werte+k)->beta);
	 *(R_coeff+j)=(*(R_coeff+j))+(*(h_werte+2*k+1))*((f_werte+k)->r);
      }
      *(L_coeff+j)=(*(L_coeff+j))*2.0/(n+1);
      *(B_coeff+j)=(*(B_coeff+j))*2.0/(n+1);
      *(R_coeff+j)=(*(R_coeff+j))*2.0/(n+1);
   }

   delete[] h_werte;
   delete[] f_werte;
}

void tschebyscheff::tschebyscheffpoly( PLANET planet, double a, double b, int n, sphaer (*f)(PLANET,double) )
/*
    Diese Funktion berechnet die Koeffizienten des TP zur Funktion f
    vom Grad n im Intervall [a,b]. f wird dabei n+1 mal aufgerufen.
*/
{
   int k,j;
   double* h_werte = new double[ 2*n+2 ];
   sphaer* f_werte = new sphaer[ n+1 ];

   if (L_coeff) { delete[] L_coeff; L_coeff = NULL; }
   if (B_coeff) { delete[] B_coeff; B_coeff = NULL; }
   if (R_coeff) { delete[] R_coeff; R_coeff = NULL; }

   L_coeff = new double[ n+1 ];
   B_coeff = new double[ n+1 ];
   R_coeff = new double[ n+1 ];

   if ((h_werte==NULL)||(f_werte==NULL)||
       (L_coeff==NULL)||(B_coeff==NULL)||(R_coeff==NULL))
   {
      fprintf(stderr,"Kein Speicher fuer TP!");
      exit(1);
   }

   Grad=n;
   a_plus_b_2 =(a+b)*0.5;
   b_minus_a_2=(b-a)*0.5;
   *h_werte=1.0;                                            // Hilfstabelle fuer
   *(h_werte+1)=cos(PI/(2*n+2));                      // Stuetzstellen
   for (k=2;k<=2*n+1;k++)                             // x(k)=h_werte(2k+1)
      *(h_werte+k)=2*(*(h_werte+1))*(*(h_werte+k-1))-(*(h_werte+k-2));
   for (k=0;k<=n;k++)                                 // f an Stuetzstellen
   {
      *(f_werte+k)=(*f)(planet,*(h_werte+2*k+1)*b_minus_a_2+a_plus_b_2); //  auswerten
   }
   for (k=1;k<=n;k++) // lambda stetig machen in [-2*Pi,+2*Pi]
   {
      if (((f_werte+k-1)->lambda)<((f_werte+k)->lambda))
	 ((f_werte+k)->lambda)-=M_2PI;
   }
   for (j=0;j<=n;j++)                                 // T-Koeffizienten
   {                                                  // berechenen
      *(h_werte+1)=cos(PI*j/(2*n+2));                     // Hilfstabelle fuer TP
      for (k=2;k<=2*n+1;k++)
	*(h_werte+k)=2*(*(h_werte+1))*(*(h_werte+k-1))-(*(h_werte+k-2));      // T(j,x(k))=h_werte(2k+1)
      *(L_coeff+j)=0;
      *(B_coeff+j)=0;
      *(R_coeff+j)=0;
      for (k=0;k<=n;k++)
      {
	 *(L_coeff+j)=(*(L_coeff+j))+(*(h_werte+2*k+1))*((f_werte+k)->lambda);
	 *(B_coeff+j)=(*(B_coeff+j))+(*(h_werte+2*k+1))*((f_werte+k)->beta);
	 *(R_coeff+j)=(*(R_coeff+j))+(*(h_werte+2*k+1))*((f_werte+k)->r);
      }
      *(L_coeff+j)=(*(L_coeff+j))*2.0/(n+1);
      *(B_coeff+j)=(*(B_coeff+j))*2.0/(n+1);
      *(R_coeff+j)=(*(R_coeff+j))*2.0/(n+1);
   }

   delete[] h_werte;
   delete[] f_werte;
}
