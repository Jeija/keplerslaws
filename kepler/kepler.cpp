/***************************************************************************
                          kepler.cpp  -  description
                             -------------------
    begin                : Sat Apr 15 2000
    copyright            : (C) 2000 by tkramer
    email                : tkramer@ph.tum.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

static char rcsid[] = "$Id: kepler.cpp,v 1.1 2000/04/16 08:35:23 tkramer Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "astromat.h"
#include "definiti.h"
#include "physdat.h"
#include "orbit.h"
#include "position.h"
#include "tscheby.h"
#include "gv_comm.h"
#include "kepler.h"

// #include <dmalloc.h>

char             PLANET_NAME[][8]={"Sonne",
                                   "Merkur",
                                   "Venus",
                                   "Erde",
                                   "Mars",
                                   "Jupiter",
                                 /* 1234567 */
                                   "Saturn",
                                   "Uranus",
                                   "Neptun",
                                   "Pluto",
				   "Mond"};

/* Erzeugt aus s den 'Blanklosen' String t */
char *strdel(char *s, char *t)
{  short i,j=0,l;

   l=strlen(s);
   for(i=0;i<l;i++)
      if (s[i]!=' ') t[j++]=s[i];
   t[j]='\0';
   return t;
}

/* Fuellt den String mit Leerzeichen auf len auf und setzt s[len]=0; */
char *strfill( char *s, short len)
{  short i=strlen(s);

   for(;i<len;i++)
      s[i]=' ';
   s[len]=0;
   return s;
}

/* Entfernt alle Blanks NACH dem letzen nicht-Blank Zeichen */
char *strcut( char *s )
{
   short i=strlen(s)-1;

   for(;(i>=0)&&(s[i]==' ');i--)
      s[i]=0;

   return s;
}

/* laden von EPHEM VERSION 4.0 Bahnelementdateien */
short elemente_laden( const char datei_name[], ELEMENT& e )
{
   FILE *f;
   int dummy;
   char string[100],a_string[20];

   f=fopen(datei_name,"rt");
   if (f==NULL)
   {
      fprintf(stderr,"Kann Datei %s nicht oeffnen!\n",datei_name);
      return 0;
   }
   fgets(string,100,f);
   if (strcmp("EPHEM4\n",string)!=0) // keine EPHEM 4.0 Bahnelementdatei
   {
      fclose(f);
      fprintf(stderr,"Die Datei %s ist keine EPHEM 4.x Bahnelementedatei.\n",datei_name);
      return 0;
   }

   fgets(string,100,f);
   sscanf(string,"%d",&dummy);
   e.komet=AUS_EIN(dummy % 2);
   fgets(string,100,f);
   strcut(string);
   string[strlen(string)-1]=0; /* '\n' rauswerfen */
   strcut(string);
   strcpy(e.name,string);
   fgets(string,100,f);
   sscanf(string,"%s",a_string);
   strdel(a_string,string);
   if (strcmp("des_datums",string)==0)
   {
      e.a_jd=0.0;
      e.a_mo=des_datums;
   }
   else if (strcmp("j2000",string)==0)
   {
      e.a_jd=J2000;
      e.a_mo=j2000;
   }
   else if (strcmp("b1950",string)==0)
   {
      e.a_jd=B1950;
      e.a_mo=b1950;
   }
   else
   {
      sscanf(string,"%lg",&e.a_jd);
      e.a_mo=anderes;
   }
   fgets(string,100,f); sscanf(string,"%lg",&e.e);
   fgets(string,100,f); sscanf(string,"%lg",&e.q);
   fgets(string,100,f); sscanf(string,"%lg",&e.M);
   fgets(string,100,f); sscanf(string,"%lg",&e.t);
   fgets(string,100,f); sscanf(string,"%lg",&e.om);
   fgets(string,100,f); sscanf(string,"%lg",&e.kl);
   fgets(string,100,f); sscanf(string,"%lg",&e.i);
   fgets(string,100,f); sscanf(string,"%lg",&e.mag0);
   fgets(string,100,f); sscanf(string,"%lg",&e.mag1);

   if (e.e<1.0)       /* Elliptische Bahn */
   {
      e.a=e.q/(1.0-e.e);
      e.n=M_2PI/sqrt(e.a*e.a*e.a*4.0*PI*PI/KGAUSS2);
   }
   else if (e.e>1.0)  /* Hyperbolische Bahn */
   {
      e.a=e.q/(1.0-e.e);
   }
   if (e.M==0.0)
      e.t0=e.t;

   e.M  *= DEG2RAD;
   e.om *= DEG2RAD;
   e.kl *= DEG2RAD;
   e.i  *= DEG2RAD;
   e.pi=mod(e.om+e.kl,M_2PI);

   fclose(f);

   return 1;
}
