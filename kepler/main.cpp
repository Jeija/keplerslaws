/***************************************************************************
                          main.cpp  -  Hauptprogramm
                             -------------------
    begin                : Sat Apr 15 13:49:34 CEST 2000
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
/* Fuer manche Systeme brauchen wir ein "usleep(Mikrosekunden)" Ersatz */
#ifdef WIN32
void usleep( long d ) { for(long i=0;i<d;i++,asin(i/d)); }
#else
#include <unistd.h>
#endif

#include "astromat.h"
#include "einaus.h"
#include "definiti.h"
#include "physdat.h"
#include "orbit.h"
#include "position.h"
#include "tscheby.h"
#include "monddat.h"
#include "gv_comm.h"
#include "texture.h"
#include "kepler.h"

static void verbinden(int id1, int id2);
static vektor pos_id( double JD, int i );
static int name2id( const char* name );
static void nameOnOff( const char* name, int on );
static void allesLoesen( void );
static void load_ImageText(char *filename);
static void load_ImageUntertitel(char *filename);
static void load_ImageWeiter(char *filename);
static void OgErstesGesetz( void );
static void OgZweitesGesetz( void );
static void OgDrittesGesetz( void );
static void OgGeozentrisch( void );
static void OgSchleife( void );
static void GlutGesetz11( int );
static void GlutGesetz12( int );
static void GlutGesetz31( int );
static void GlutGesetz32( int );
static void GlutIdleAnimate( void );
static void GlutIdleSchoner( void );
static void GlutTastatur( int key, int , int );
static void GlutTastatur2( int key, int , int );
static void GlutTastaturIdle( int key, int , int );
static void GlutSchoner( int value );
static void GlutImageWeiter( int value );
static void ShowText( char *str );
static void ShowUntertitel( char *str );
static void ShowGrafik( void );
static void GlutReshapeUntertitel(int w, int h);
static void GlutDisplayUntertitel( void );
static void GlutReshapeText(int w, int h);
static void GlutDisplayText( void );
static void GlutReshapeGrafik( int, int );
static void GlutDisplayGrafik( void );
static void SetGlutSpecialFuncAll( void (*func)(int key, int x,int y) );
static void SetGlutIdleFuncText( void (*func)(void) );
static void init2( void );
static void init( void );

#define MAX_VERBINDUNGEN 40

static char rcsid[] = "$Id: main.cpp,v 1.3 2000/04/16 08:41:09 tkramer Exp $";

static char *la_sprache_string[3]={"DE","EN","UN"};
enum LA_SPRACHE { LA_DE, LA_EN, LA_ENDE };

// enum ID_MODUS { ID_KAMPOS,ID_KAMDIR,ID_FIXPOS,ID_FIXDIR,ID_MODUS_ENDE };

static int LIST_SPEZIAL_SUB,LIST_SPEZIAL,LIST_STENCIL;
static PLANETEN planeten;
static ELEMENTE elemente;
static KAMERA_MODUS kamera_modus;
static OBJEKT_MODUS objekt_modus;
static int id_modus,id_kamdir,id_kampos,id_fixpos,id_fixdir,id_sonne,id_objekt1,id_objekt2;
static double jde_start,schrittweite;
static double jde;
static int f_animate;
static int GlutWinText,GlutWinUntertitel,GlutWinGrafik;
static int verbindungen[MAX_VERBINDUNGEN];
static double RADIUS;
static double AEQUINOKTIUM;
static AUS_EIN start,stencil;
static int IndexImageUntertitel;
static AUS_EIN schoner;
static time_t zeit;
static double rx,ry,rz;
static int n_F12;
static enum LA_SPRACHE la_sprache=LA_DE;

static unsigned *ImageText_data=NULL;
static int ImageText_width, ImageText_height, ImageText_components;
static GLenum ImageText_format = GL_RGBA;

static unsigned *ImageUntertitel_data=NULL;
static int ImageUntertitel_width, ImageUntertitel_height, ImageUntertitel_components;
static GLenum ImageUntertitel_format = GL_RGBA;

static unsigned *ImageWeiter_data=NULL;
static int ImageWeiter_width, ImageWeiter_height, ImageWeiter_components;
static GLenum ImageWeiter_format = GL_RGBA;

const int ZEIT_SCHONER=120; // Nach ZEIT_SCHONER Sekunden wird das Hauptmenu gezeigt
const int SCROLL_INTERVALL=100; // jede intervall Sekunden den Schonertext weiterscrollen

static GLfloat cube[8][3] =
{{-0.5,-0.5,-0.5},
 { 0.5,-0.5,-0.5},
 { 0.5,-0.5, 0.5},
 {-0.5,-0.5, 0.5},
 { 0.5, 0.5,-0.5},
 { 0.5, 0.5, 0.5},
 {-0.5, 0.5, 0.5},
 {-0.5, 0.5,-0.5}};

static int faceIndex[6][4] =
{{0, 1, 2, 3},
 {1, 4, 5, 2},
 {4, 7, 6, 5},
 {7, 0, 3, 6},
 {3, 2, 5, 6},
 {7, 4, 1, 0}};

static GLfloat colorIndex[6][3] =
{{0.800,0.098,0.098},
 {0.098,0.647,0.400},
 {0.098,0.098,0.800},
 {0.898,0.600,0.000},
 {0.000,0.600,0.800},
 {0.498,0.000,0.898}};

static void verbinden(int id1, int id2)
{
   int i,done=0;

   if (id1>id2) { int tmp; tmp=id2; id2=id1; id1=tmp; }
   for(i=0;i<MAX_VERBINDUNGEN;i+=2)
   {
      if ((verbindungen[i]==id1)&&(verbindungen[i+1]==id2))
      {
         verbindungen[i]=verbindungen[i+1]=-1;
	 done=1;
	 break;
      }
   }
   if (done!=1)
   {
      for(i=0;i<MAX_VERBINDUNGEN;i+=2)
      {
         if (verbindungen[i]==-1)
	 {
	    verbindungen[i]=id1;
	    verbindungen[i+1]=id2;
	    break;
	 }
      }
   }
}

static vektor pos_id( double JD, int i )
{
   vektor v;
   if (i<elemente.n)
   {
      vektor P,Q,vel;
      double d0;

      if (elemente.e[i].e<1.0)
         d0=elemente.e[i].t-elemente.e[i].M/elemente.e[i].n;
      else
         d0=elemente.e[i].t0;
      gaussvek(elemente.e[i].kl,elemente.e[i].i,elemente.e[i].om,P,Q);
      v=kepler(elemente.e[i].e,elemente.e[i].q,JD,d0,P,Q,vel);

      if (elemente.e[i].a_mo!=des_datums)
         v=m_praez_ekl(elemente.e[i].a_jd,AEQUINOKTIUM)*v;
      else
	 v=m_praez_ekl(JD,AEQUINOKTIUM)*v;
   }
   else if (i<elemente.n+planeten.n) // Planet
   {
      v=m_praez_ekl(JD,AEQUINOKTIUM)*planet_lo(planeten.p[i-elemente.n],JD);
   }
   else
   {
      v=vektor(0,1,0);
   }
   return v;
}

static int name2id( const char* name )
{
   int i;

   for(i=0;i<elemente.n;i++)
   {
      if (strcmp(elemente.e[i].name,name)==0)
         return i;
   }
   for(i=0;i<planeten.n;i++)
   {
      if (strcmp(PLANET_NAME[planeten.p[i]],name)==0)
         return i+elemente.n;
   }
   if (strcmp("(0,1,0)",name)==0)
      return elemente.n+planeten.n;
   return -1;
}

static void nameOnOff( const char* name, int on )
{
   int i;

   for(i=0;i<elemente.n;i++)
   {
      if (strcmp(elemente.e[i].name,name)==0)
         elemente.on[i]=on;
   }
   for(i=0;i<planeten.n;i++)
   {
      if (strcmp(PLANET_NAME[planeten.p[i]],name)==0)
         planeten.on[i]=on;
   }
}

static void allesLoesen( void )
{
   int i;

   for(i=0;i<MAX_VERBINDUNGEN;i++)
   {
      verbindungen[i]=-1;
   }
}

static void load_ImageText(char *filename)
{
  if (ImageText_data!=NULL)
     free(ImageText_data);

  if (filename) {
    ImageText_data = read_texture(filename, &ImageText_width, &ImageText_height, &ImageText_components);
    if (ImageText_data == NULL) {
      fprintf(stderr, "Error: Can't load ImageText file \"%s\".\n",
        filename);
      exit(1);
    }
    if (ImageText_components < 3 || ImageText_components > 4) {
      fprintf(stderr,"must be RGB or RGBA ImageText\n");
      exit(1);
    }
  }
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
}

static void load_ImageUntertitel(char *filename)
{
  if (ImageUntertitel_data!=NULL)
     free(ImageUntertitel_data);

  if (filename) {
    ImageUntertitel_data = read_texture(filename, &ImageUntertitel_width, &ImageUntertitel_height, &ImageUntertitel_components);
    if (ImageUntertitel_data == NULL) {
      fprintf(stderr, "Error: Can't load ImageUntertitel file \"%s\".\n",
        filename);
      exit(1);
    }
    if (ImageUntertitel_components < 3 || ImageUntertitel_components > 4) {
      fprintf(stderr,"must be RGB or RGBA ImageUntertitel\n");
      exit(1);
    }
  }
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
}

static void load_ImageWeiter(char *filename)
{
  if (ImageWeiter_data!=NULL)
     free(ImageWeiter_data);

  if (filename) {
    ImageWeiter_data = read_texture(filename, &ImageWeiter_width, &ImageWeiter_height, &ImageWeiter_components);
    if (ImageWeiter_data == NULL) {
      fprintf(stderr, "Error: Can't load ImageWeiter file \"%s\".\n",
        filename);
      exit(1);
    }
    if (ImageWeiter_components < 3 || ImageWeiter_components > 4) {
      fprintf(stderr,"must be RGB or RGBA ImageWeiter\n");
      exit(1);
    }
  }
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
}

static void OgErstesGesetz( void )
{
/* Erstes Keplersches Gesetz: Die Planeten bewegen sich auf Ellipsenbahnen
   um die Sonne, die im Brennpunkt steht.
*/
   int i,id;

   for(i=0;i<elemente.n;i++)
   {
      elemente.on[i]=1;
   }
   for(i=0;i<planeten.n;i++)
   {
      planeten.on[i]=1;
   }
   nameOnOff("Oterma",0);
   nameOnOff("Erde.",0);
   allesLoesen();
   rx=ry=rz=0.0;

   jde=juldat(1985,9,1,0.0);
   if ((id=name2id("Sonne"))  >=0) id_fixpos=id;
   if ((id=name2id("(0,1,0)"))>=0) id_fixdir=id;

   glNewList(LIST_SPEZIAL_SUB,GL_COMPILE);
   glEndList();

   glNewList(LIST_SPEZIAL,GL_COMPILE);
   for(i=0;i<elemente.n;i++)
   {
      if (elemente.on[i]==1)
      {
         if (i==name2id("Mars"))
	 {
            glColor3d(1.0,0.0,0.0);
	 }
	 else
	 {
            glColor3d(1.0,1.0,1.0);
	 }
         og_orbit(elemente.e[i],AEQUINOKTIUM);
      }
   }
   glEndList();

   stencil=AUS;
}

static void OgZweitesGesetz( void )
{
/* Zweites Keplersches Gesetz: Die Fahrstrahl Sonne-Planet ueberstreicht
   in gleichen Zeitintervallen gleiche Flaechen der Ellipsen
*/
   int i,id1,id2;

   for(i=0;i<elemente.n;i++)
   {
      elemente.on[i]=1;
   }
   for(i=0;i<planeten.n;i++)
   {
      planeten.on[i]=1;
   }
   nameOnOff("Oterma",0);
   nameOnOff("Erde.",0);
   allesLoesen();
   rx=ry=rz=0.0;
   jde=jde_start=juldat(1985,9,1,0.0);

   if ((id1=name2id("Sonne"))  >=0) id_fixpos=id1;
   if ((id1=name2id("(0,1,0)"))>=0) id_fixdir=id1;
   if (((id1=name2id("Sonne"))>=0)&&((id2=name2id("Halley"))>=0))
   {
      id_sonne=id1;
      id_objekt1=id2;
      verbinden(id1,id2);
   }

   glNewList(LIST_SPEZIAL_SUB,GL_COMPILE);
   glEndList();

   glNewList(LIST_SPEZIAL,GL_COMPILE);
   glColor3d(1.0,1.0,1.0);
   for(i=0;i<elemente.n;i++)
   {
      if (elemente.on[i]==1)
      {
         og_orbit(elemente.e[i],AEQUINOKTIUM);
      }
   }
   glEndList();

   stencil=AUS;
}

static void OgDrittesGesetz( void )
{
/* Drittes Keplersches Gesetz: Die Kuben der grossen Halbachsen verhalten
   sich wie die Quadrate der Umlaufzeiten
*/
   int i,id1,id2;

   for(i=0;i<elemente.n;i++)
   {
      elemente.on[i]=0;
   }
   for(i=0;i<planeten.n;i++)
   {
      planeten.on[i]=0;
   }
   nameOnOff("Sonne",1);
   nameOnOff("Oterma",1);
   nameOnOff("Erde.",1);
   allesLoesen();
   rx=ry=rz=0.0;

   jde=jde_start=juldat(2000,1,1,0.0);

   if ((id1=name2id("Sonne"))  >=0) id_fixpos=id1;
   if ((id1=name2id("(0,1,0)"))>=0) id_fixdir=id1;
   if (((id1=name2id("Sonne"))>=0)&&((id2=name2id("Oterma"))>=0))
   {
      id_sonne=id1;
      id_objekt1=id2;
      verbinden(id1,id2);
   }
   if (((id1=name2id("Sonne"))>=0)&&((id2=name2id("Erde."))>=0))
   {
      id_sonne=id1;
      id_objekt2=id2;
      verbinden(id1,id2);
   }

   glNewList(LIST_SPEZIAL_SUB,GL_COMPILE);
   glEndList();

   glNewList(LIST_SPEZIAL,GL_COMPILE);
   glColor3d(1.0,1.0,1.0);
   for(i=0;i<elemente.n;i++)
   {
      if (elemente.on[i]==1)
      {
         og_orbit(elemente.e[i],AEQUINOKTIUM);
      }
   }
   glEndList();

   stencil=AUS;
}

static void OgGeozentrisch( void )
{
/* Das Geozentrische Weltbild.
*/
   int i,id;

   for(i=0;i<elemente.n;i++)
   {
      elemente.on[i]=1;
   }
   for(i=0;i<planeten.n;i++)
   {
      planeten.on[i]=1;
   }
   nameOnOff("Halley",0);
   nameOnOff("Oterma",0);
   nameOnOff("Erde.",0);
   allesLoesen();
   rx=ry=rz=0.0;
   jde=juldat(1985,11,1,0.0);
   if ((id=name2id("Erde"))  >=0)  id_fixpos=id;
   if ((id=name2id("(0,1,0)"))>=0) id_fixdir=id;

   glNewList(LIST_SPEZIAL_SUB,GL_COMPILE);
   glEndList();

   glNewList(LIST_SPEZIAL,GL_COMPILE);
   glColor3d(1.0,1.0,1.0);
   for(i=0;i<elemente.n;i++)
   {
      if (elemente.on[i]==1)
      {
         if (i!=id_fixpos)
         {
	    og_orbit(elemente.e[i],AEQUINOKTIUM);
         }
      }
   }
   glEndList();

   glNewList(LIST_STENCIL,GL_COMPILE);
   glEndList();
   glClear(GL_STENCIL_BUFFER_BIT);

   stencil=AUS;
}

static void OgSchleife( void )
{
/* Schleifenbewegung der Planeten
*/
   int i,id;

   for(i=0;i<elemente.n;i++)
   {
      elemente.on[i]=1;
   }
   for(i=0;i<planeten.n;i++)
   {
      planeten.on[i]=1;
   }
   nameOnOff("Halley",0);
   nameOnOff("Oterma",0);
   nameOnOff("Erde.",0);
   allesLoesen();
   rx=  0.0;
   ry=-75.0;
   rz= 75.0;

   jde=jde_start=juldat(1992,8,1,0.0);
   if ((id=name2id("Sonne"))  >=0) id_fixpos=id;
   if ((id=name2id("(0,1,0)"))>=0) id_fixdir=id;
   if ((id=name2id("Erde"))   >=0) id_objekt1=id;
   if ((id=name2id("Mars"))   >=0) id_objekt2=id;

   glNewList(LIST_SPEZIAL_SUB,GL_COMPILE);
   glEndList();

   glNewList(LIST_SPEZIAL,GL_COMPILE);
   glColor3d(1.0,1.0,1.0);
   for(i=0;i<elemente.n;i++)
   {
      if (elemente.on[i]==1)
      {
         og_orbit(elemente.e[i],AEQUINOKTIUM);
      }
   }
   og_sterne(AUS,jde,AEQUINOKTIUM);
   glEndList();

   glNewList(LIST_STENCIL,GL_COMPILE);
   glEndList();
   glClear(GL_STENCIL_BUFFER_BIT);

   stencil=AUS;
}

static void GlutGesetz11( int )
{
   glutSetWindow(GlutWinGrafik);

   int i;

   glNewList(LIST_SPEZIAL_SUB,GL_COMPILE);
   glColor3d(1.0,0.0,0.0);
   glLineWidth(3);
   if ((i=name2id("Mars"))>=0)
   {
      og_focus(elemente.e[i],AEQUINOKTIUM);
   }
   glLineWidth(1);
   glEndList();

   SetGlutSpecialFuncAll(GlutTastaturIdle);
}

static void GlutGesetz12( int )
{
   glutSetWindow(GlutWinGrafik);

   glPushMatrix();
   for(ry=0.0;ry>=-90.0;ry-=1.0)
   {
      GlutDisplayGrafik();
   }
   glPopMatrix();

   SetGlutSpecialFuncAll(GlutTastaturIdle);
}


static void GlutGesetz31( int )
{
   SetGlutIdleFuncText(NULL);
   glutSetWindow(GlutWinGrafik);

   vektor v;

   glNewList(LIST_SPEZIAL_SUB,GL_COMPILE);

   v=vektor(0.5*elemente.e[id_objekt1].a,+0.5*elemente.e[id_objekt1].a+elemente.e[id_objekt1].q,0.5*elemente.e[id_objekt1].a);
   glPushMatrix();
   glTranslated(v[0],v[1],v[2]);
   glutWireCube(elemente.e[id_objekt1].a);
   glPopMatrix();

   v=vektor(0.5*elemente.e[id_objekt2].a,+0.5*elemente.e[id_objekt2].a-elemente.e[id_objekt2].q,0.5*elemente.e[id_objekt2].a);
   glPushMatrix();
   glTranslated(v[0],v[1],v[2]);
   glutWireCube(elemente.e[id_objekt2].a);
   glPopMatrix();

   glEndList();

   glPushMatrix();
   for(rx=0,ry=0,rz=0;ry<=360;rx+=1.0,ry+=2.0,rz+=1.0)
   {
      GlutDisplayGrafik();
   }
   glPopMatrix();
   ry=360.0;
   rx=180.0;
   rz=180.0;

   SetGlutSpecialFuncAll(GlutTastaturIdle);
}

static void GlutGesetz32( int )
{
   SetGlutIdleFuncText(NULL);
   glutSetWindow(GlutWinGrafik);

   vektor v,v0,v_delta;

   v0=vektor(0.5*elemente.e[id_objekt2].a,+0.5*elemente.e[id_objekt2].a-elemente.e[id_objekt2].q,0.5*elemente.e[id_objekt2].a);
   v_delta=-vektor(0.5*elemente.e[id_objekt2].a,+0.5*elemente.e[id_objekt2].a-elemente.e[id_objekt2].q,0.5*elemente.e[id_objekt2].a);
   v_delta+=vektor(0.5*elemente.e[id_objekt1].a,+0.5*elemente.e[id_objekt1].a+elemente.e[id_objekt1].q,0.5*elemente.e[id_objekt1].a);
   v_delta+=vektor(0.5*elemente.e[id_objekt2].a,+0.5*elemente.e[id_objekt2].a-2.0*elemente.e[id_objekt2].a,0.5*elemente.e[id_objekt2].a);


   // bewege kleinen Wuerfel
   int i,x,y,n;
   n=100;
   for(i=0;i<=n;i++)
   {
      glNewList(LIST_SPEZIAL_SUB,GL_COMPILE);
      glPushMatrix();
      glTranslated(v0[0]+v_delta[0]*i/n,v0[1]+v_delta[1]*i/n,v0[2]+v_delta[2]*i/n);
      glScaled(elemente.e[id_objekt2].a,elemente.e[id_objekt2].a,elemente.e[id_objekt2].a);
      for (x=0;x<6;x++)
      {
	 glColor3fv((GLfloat *)colorIndex[x]);
	 glBegin(GL_POLYGON);
	 for (y=0;y<4;y++)
	 {
            glVertex3fv((GLfloat *)cube[faceIndex[x][y]]);
	 }
	 glEnd();
      }
      glPopMatrix();
      v=vektor(0.5*elemente.e[id_objekt1].a,+0.5*elemente.e[id_objekt1].a+elemente.e[id_objekt1].q,0.5*elemente.e[id_objekt1].a);
      glPushMatrix();
      glTranslated(v[0],v[1],v[2]);
      glColor3d(1,1,1);
      glutWireCube(elemente.e[id_objekt1].a);
      glPopMatrix();
      glEndList();
      GlutDisplayGrafik();
   }

   {
      glNewList(LIST_SPEZIAL_SUB,GL_COMPILE);
      glPushMatrix();
      glTranslated(v0[0]+v_delta[0],v0[1]+v_delta[1],v0[2]+v_delta[2]);
      glScaled(elemente.e[id_objekt2].a,elemente.e[id_objekt2].a,elemente.e[id_objekt2].a);
      for (x=0;x<6;x++)
      {
	 glColor3fv((GLfloat *)colorIndex[x]);
	 glBegin(GL_POLYGON);
	 for (y=0;y<4;y++)
	 {
            glVertex3fv((GLfloat *)cube[faceIndex[x][y]]);
	 }
	 glEnd();
      }
      glPopMatrix();
      v=vektor(0.5*elemente.e[id_objekt1].a,+0.5*elemente.e[id_objekt1].a+elemente.e[id_objekt1].q,0.5*elemente.e[id_objekt1].a);
      glPushMatrix();
      glTranslated(v[0],v[1],v[2]);
      glColor3d(1,1,1);
      glutWireCube(elemente.e[id_objekt1].a);
      glPopMatrix();
      glPushMatrix();
      glTranslated(v[0],v[1],v[2]);
      glScaled(0.5,0.5,0.5);
      glBegin(GL_LINES);
      for(i=0;i<3;i++)
      {
	 for(x=-2;x<=2;x++)
	 {
	    for(y=-2;y<=2;y++)
	    {
		switch(i)
		{
	           case 0: glVertex3i(x,y,-2);
		           glVertex3i(x,y,+2);
		           break;
		   case 1: glVertex3i(-2,x,y);
		           glVertex3i(+2,x,y);
			   break;
		   case 2: glVertex3i(x,-2,y);
		           glVertex3i(x,+2,y);
		           break;
		}
	    }
	 }
      }
      glPopMatrix();
      glEnd();
      glEndList();
      GlutDisplayGrafik();
   }

   SetGlutSpecialFuncAll(GlutTastaturIdle);
}

static void GlutIdleAnimate( void )
{
   vektor v[5];
   double delta;
   static vektor v_alt=vektor(0,0,0);

   glutSetWindow(GlutWinGrafik);
   switch(f_animate)
   {
      case 2:
         delta=42.0;
	 if (fabs(mod(jde-jde_start,delta))==0.0)
         {
            int i,n;
	
	    n=int((jde-jde_start)/delta);
            if (n>20)
	    {
	       break;
	    }
	    glNewList(LIST_SPEZIAL_SUB,GL_COMPILE);
       	    for(i=0;i<=n;i++)
	    {
	       v[0]=pos_id(jde_start+i*delta,id_sonne);
	       v[1]=pos_id(jde_start+i*delta-0.0/3.0*delta/2.0,id_objekt1);
	       v[2]=pos_id(jde_start+i*delta-1.0/3.0*delta/2.0,id_objekt1);
	       v[3]=pos_id(jde_start+i*delta-2.0/3.0*delta/2.0,id_objekt1);
	       v[4]=pos_id(jde_start+i*delta-3.0/3.0*delta/2.0,id_objekt1);
	       og_polygon(5,v);
            }
	    glEndList();
         }
         break;
      case 3:
         if ((jde-jde_start)==0.0)
         {
	    glNewList(LIST_SPEZIAL_SUB,GL_COMPILE);
            og_strahl(pos_id(jde,id_sonne),pos_id(jde,id_objekt1));
            og_strahl(pos_id(jde,id_sonne),pos_id(jde,id_objekt2));
	    glEndList();
         }
         break;
      case 4:
         {
  	    stencil=EIN;
	    glNewList(LIST_STENCIL,GL_COMPILE);
	    glPointSize(1.0);
	    glBegin(GL_POINTS);
            int i;
	    for(i=0;i<elemente.n;i++)
	    {
	       if (elemente.on[i]==1)
	       {
		  v[0]=pos_id(jde,i);
		  glVertex3d(v[0][0],v[0][1],v[0][2]);
               }
	    }
	    for(i=0;i<planeten.n;i++)
	    {
	       if (planeten.on[i]==1)
	       {
		  v[0]=pos_id(jde,i+elemente.n);
		  glVertex3d(v[0][0],v[0][1],v[0][2]);
               }
	    }
	    glEnd();
	    glEndList();
	 }
	 break;
      case 5:
         {
	    v[0]=pos_id(jde,id_objekt2)-pos_id(jde,id_objekt1);
	    v[0]/=sqrt(v[0]|v[0]);
            v[0]=pos_id(jde,id_objekt1)+70.0*v[0];
            delta=1.0;
	    if (IndexImageUntertitel==0)
 	    {
	       glNewList(LIST_SPEZIAL_SUB,GL_COMPILE);
               glColor3d(1,1,0);
               og_strahl(pos_id(jde,id_objekt1),v[0]);
	       glEndList();
	       v_alt=vektor(0,0,0);
	       stencil=AUS;
	    }
	    else
	    {
	       glNewList(LIST_SPEZIAL_SUB,GL_COMPILE);
               glColor3d(1,0,0);
               og_strahl(pos_id(jde,id_objekt1),v[0]);
	       glEndList();
               stencil=EIN;
               if ((v_alt|v_alt)>0)
	       {
  		  glNewList(LIST_STENCIL,GL_COMPILE);
		  glBegin(GL_LINE_STRIP);
		  glVertex3d(v_alt[0],v_alt[1],v_alt[2]);
		  glVertex3d(v[0][0],v[0][1],v[0][2]);
		  glEnd();
		  glEndList();
	       }
  	       v_alt=v[0];
            }
	 }
         break;
   }
   jde+=schrittweite;
   GlutDisplayGrafik();
}

static void GlutIdleSchoner( void )
{
   static time_t talt=time(NULL);
   static int i=0;
   time_t t;
   int d;

   if (schoner==AUS)
      return;

   t=time(NULL);
   if (((t-zeit)%SCROLL_INTERVALL==0)&&(t!=talt))
   {
      talt=t;
      glutSetWindow(GlutWinText);
      glClear(GL_COLOR_BUFFER_BIT);
      if (ImageText_data!=NULL)
      {
         glRasterPos4d(0,glutGet((GLenum)GLUT_WINDOW_HEIGHT)-ImageText_height-i,0,1);
	 glRasterPos4d(0,i,0,1);
	 glDrawPixels(ImageText_width, ImageText_height,
	 ImageText_format, GL_UNSIGNED_BYTE, ImageText_data);
      }
      d=glutGet((GLenum)GLUT_WINDOW_HEIGHT)-ImageText_height;
      if (d>5)
      {
         i=(i+5)%(d);
      }
      else
      {
         fprintf(stderr,"GlutIdleSchoner: Bild zu gross, Programmende\n");
	 exit(1);
      }
   }
//++ CPU
   usleep(1000L);
}

// Tastatur, wenn Hauptauswahl oder Bildschirmschoner aktiv
static void GlutTastatur( int key, int , int )
{
   char filename[255];

   SetGlutSpecialFuncAll(NULL);

   zeit=time(NULL);
   glutTimerFunc(1000*ZEIT_SCHONER,GlutSchoner,0);

   if (schoner==EIN)
   {
//+++ Spracherkennung
      switch(key)
      {
         case GLUT_KEY_F11 : la_sprache=LA_EN; break;
         case GLUT_KEY_F12 : la_sprache=LA_DE; break;
      }
      sprintf(filename,"CropWeiter_%s.rgb",la_sprache_string[int(la_sprache)]);
      load_ImageWeiter(filename);

      schoner=AUS;
      sprintf(filename,"Text0_%s.rgb",la_sprache_string[int(la_sprache)]);
      ShowText(filename);
      SetGlutSpecialFuncAll(GlutTastatur);
      return;
   }
   switch(key)
   {
      case GLUT_KEY_F1 :
         f_animate=1;
         sprintf(filename,"Text1_%s.rgb",la_sprache_string[int(la_sprache)]);
         ShowText(filename);
         glutSetWindow(GlutWinGrafik);
         OgErstesGesetz();
	 SetGlutSpecialFuncAll(GlutTastatur2);
	 break;
      case GLUT_KEY_F2 :
         f_animate=2;
         sprintf(filename,"Text2_%s.rgb",la_sprache_string[int(la_sprache)]);
         ShowText(filename);
         glutSetWindow(GlutWinGrafik);
         OgZweitesGesetz();
	 SetGlutSpecialFuncAll(GlutTastatur2);
	 break;
      case GLUT_KEY_F3 :
         f_animate=3;
         sprintf(filename,"Text3_%s.rgb",la_sprache_string[int(la_sprache)]);
         ShowText(filename);
         glutSetWindow(GlutWinGrafik);
         OgDrittesGesetz();
	 SetGlutSpecialFuncAll(GlutTastatur2);
	 break;
      case GLUT_KEY_F4 :
         f_animate=4;
         sprintf(filename,"Text4_%s.rgb",la_sprache_string[int(la_sprache)]);
         ShowText(filename);
         glutSetWindow(GlutWinGrafik);
         OgGeozentrisch();
	 SetGlutSpecialFuncAll(GlutTastatur2);
	 break;
      case GLUT_KEY_F5 :
         f_animate=5;
         sprintf(filename,"Text5_%s.rgb",la_sprache_string[int(la_sprache)]);
         ShowText(filename);
         glutSetWindow(GlutWinGrafik);
         OgSchleife();
	 SetGlutSpecialFuncAll(GlutTastatur2);
	 break;
      case GLUT_KEY_F11 :
         la_sprache=LA_EN;
	 sprintf(filename,"CropWeiter_%s.rgb",la_sprache_string[int(la_sprache)]);
         load_ImageWeiter(filename);

         sprintf(filename,"Text0_%s.rgb",la_sprache_string[int(la_sprache)]);
         ShowText(filename);
         SetGlutSpecialFuncAll(GlutTastatur);
         return;
      case GLUT_KEY_F12 :
         la_sprache=LA_DE;
         sprintf(filename,"CropWeiter_%s.rgb",la_sprache_string[int(la_sprache)]);
         load_ImageWeiter(filename);

	 sprintf(filename,"Text0_%s.rgb",la_sprache_string[int(la_sprache)]);
         ShowText(filename);
         SetGlutSpecialFuncAll(GlutTastatur);
         return;
      default:
         SetGlutSpecialFuncAll(GlutTastatur);
	 break;
   }
}

// Tastatur nach dem ersten Infobildschirm:
// Ordne Fenster neu, zeige ersten Untertitel
static void GlutTastatur2( int key, int , int )
{
   char s[100];

   SetGlutSpecialFuncAll(NULL);
   zeit=time(NULL);
   glutTimerFunc(1000*ZEIT_SCHONER,GlutSchoner,0);
   IndexImageUntertitel=0;

   switch(key)
   {
      case GLUT_KEY_F12:
         break;
      default:
         SetGlutSpecialFuncAll(GlutTastatur2);
         return;
   }

   glutSetWindow(GlutWinGrafik);
   glutShowWindow();
   glutSetWindow(GlutWinUntertitel);
   glutShowWindow();

   if (f_animate<6)
   {
      sprintf(s,"Untertitel%d_%d_%s.rgb",
      f_animate,IndexImageUntertitel,la_sprache_string[int(la_sprache)]);
      ShowUntertitel(s);
   }
   SetGlutSpecialFuncAll(GlutTastaturIdle);
   glutTimerFunc(2000,GlutImageWeiter,0);
   n_F12=1;
   start=EIN;
}

// Tastatur, waehrend Animation laeuft
static void GlutTastaturIdle( int key, int , int )
{
   char s[100];
   FILE *f;

   SetGlutSpecialFuncAll(NULL);
   zeit=time(NULL);
   glutTimerFunc(1000*ZEIT_SCHONER,GlutSchoner,0);

   switch(key)
   {
      case GLUT_KEY_F12:
	 if (n_F12%2==0) // das Erste Mal F12
	 {
	    SetGlutIdleFuncText(NULL);
	    n_F12=(n_F12+1)%2;
	    IndexImageUntertitel++;
	    sprintf(s,"Untertitel%d_%d_%s.rgb",
	    f_animate,IndexImageUntertitel,la_sprache_string[int(la_sprache)]);

	    if ((f=fopen(s,"r"))!=NULL)
	    {
	       fclose(f);
               ShowUntertitel(s);
               glutTimerFunc(2000,GlutImageWeiter,0);
               SetGlutSpecialFuncAll(GlutTastaturIdle);
	       return;
	    }
	 }
	 else // das zweite Mal F12
	 {
            SetGlutIdleFuncText(GlutIdleAnimate);
	    n_F12=(n_F12+1)%2;
	    if (f_animate==3&&IndexImageUntertitel==1)
	    {
               glutTimerFunc(2000,GlutImageWeiter,0);
               glutTimerFunc(0,GlutGesetz31,0);
	       return;
	    }
	    else if (f_animate==3&&IndexImageUntertitel==2)
	    {
               glutTimerFunc(2000,GlutImageWeiter,0);
               glutTimerFunc(0,GlutGesetz32,0);
	       return;
	    }
	    else if (f_animate==1&&IndexImageUntertitel==1)
	    {
               glutTimerFunc(2000,GlutImageWeiter,0);
               glutTimerFunc(0,GlutGesetz11,0);
	       return;
	    }
	    else if (f_animate==1&&IndexImageUntertitel==2)
	    {
               glutTimerFunc(0,GlutGesetz12,0);
	       return;
	    }
	    else if (f_animate==5&&IndexImageUntertitel==1)
	    {
	       jde=jde_start=juldat(1992,8,1,0.0);
               glutTimerFunc(2000,GlutImageWeiter,0);
               SetGlutSpecialFuncAll(GlutTastaturIdle);
	       return;
            }
	    else
	    {
               glutTimerFunc(2000,GlutImageWeiter,0);
               SetGlutSpecialFuncAll(GlutTastaturIdle);
	       return;
	    }
	 }
	 break;
   }
   start=AUS;

   SetGlutSpecialFuncAll(NULL);
   SetGlutIdleFuncText(NULL);

   glutSetWindow(GlutWinGrafik);
   glutHideWindow();

   glutSetWindow(GlutWinUntertitel);
   glutHideWindow();

   glutSetWindow(GlutWinText);
   glutShowWindow();

   SetGlutSpecialFuncAll(GlutTastatur);

   sprintf(s,"Text0_%s.rgb",la_sprache_string[int(la_sprache)]);
   ShowText(s);
}

// System zuruecksetzen auf den ersten Bildschirm
static void GlutSchoner( int value )
{
   time_t t=time(NULL);

   if (int(t-zeit)<ZEIT_SCHONER)
   {
      return;
   }
   else if(schoner==EIN)
   {
      return;
   }
   else
   {
      start=AUS;
      schoner=EIN;

      SetGlutSpecialFuncAll(NULL);
      SetGlutIdleFuncText(NULL);

      glutSetWindow(GlutWinGrafik);
      glutHideWindow();

      glutSetWindow(GlutWinUntertitel);
      glutHideWindow();

      glutSetWindow(GlutWinText);
      glutShowWindow();

      SetGlutSpecialFuncAll(GlutTastatur);
      SetGlutIdleFuncText(GlutIdleSchoner);

      ShowText("CropTitel.rgb");
   }
}

static void GlutImageWeiter( int value )
{
   if (start==AUS)
      return;

   int window=glutGetWindow();

   glutSetWindow(GlutWinUntertitel);
   if (ImageWeiter_data!=NULL)
   {
      glRasterPos4d(glutGet((GLenum)GLUT_WINDOW_WIDTH)-ImageWeiter_width,0,0,1);
      glDrawPixels(ImageWeiter_width, ImageWeiter_height, ImageWeiter_format, GL_UNSIGNED_BYTE, ImageWeiter_data);
   }
   glutSetWindow(window);

   glFlush();
}

static void ShowText( char *str )
{
   load_ImageText(str);
   glutSetWindow(GlutWinText);
   glutShowWindow();
   GlutDisplayText();
}

static void ShowUntertitel( char *str )
{
   load_ImageUntertitel(str);
   glutSetWindow(GlutWinUntertitel);
   glutShowWindow();
   GlutDisplayUntertitel();
}

static void ShowGrafik( void )
{
   glutSetWindow(GlutWinGrafik);
   glutShowWindow();
   GlutDisplayGrafik();
}

static void GlutReshapeUntertitel(int w, int h)
{
   glViewport(0, 0, w, h);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluOrtho2D(0, w, 0, h);
   glMatrixMode(GL_MODELVIEW);
}

static void GlutDisplayUntertitel( void )
{
   glClear(GL_COLOR_BUFFER_BIT);
   if (ImageUntertitel_data!=NULL)
   {
      glRasterPos4d(0,0,0,1);
      glDrawPixels(ImageUntertitel_width, ImageUntertitel_height, ImageUntertitel_format, GL_UNSIGNED_BYTE, ImageUntertitel_data);
   }
   glFlush();
}

static void GlutReshapeText(int w, int h)
{
   glViewport(0, 0, w, h);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluOrtho2D(0, w, 0, h);
   glMatrixMode(GL_MODELVIEW);
}

static void GlutDisplayText( void )
{
   glClear(GL_COLOR_BUFFER_BIT);
   if (ImageText_data!=NULL)
   {
      glRasterPos4d(0,0,0,1);
      glDrawPixels(ImageText_width, ImageText_height, ImageText_format, GL_UNSIGNED_BYTE, ImageText_data);
   }
   glutSwapBuffers();
}

static void GlutReshapeGrafik( int, int )
{
   og_init(planeten,elemente);
}

static void GlutDisplayGrafik( void )
{
   char str[100];
   int i;

   og_camera(kamera_modus,pos_id(jde,id_kamdir),pos_id(jde,id_kampos),pos_id(jde,id_fixpos),pos_id(jde,id_fixdir));
   glRotated(rx,1,0,0);
   glRotated(ry,0,1,0);
   glRotated(rz,0,0,1);

   if (stencil==EIN) // wir haben was fuer den Stencilbuffer
   {
      glClear(GL_COLOR_BUFFER_BIT);

      glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
      glEnable(GL_STENCIL_TEST);
      glDepthMask(GL_FALSE);
      glDisable(GL_DEPTH_TEST);
      glStencilFunc(GL_ALWAYS,0,0);
      glStencilOp(GL_KEEP,GL_KEEP,GL_INCR);
      glPushMatrix();
      glCallList(LIST_STENCIL);
      glPopMatrix();

      // Male dort, wo Stencil==1
      glStencilFunc(GL_EQUAL, 1, 1);
      glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP);
      glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glLoadIdentity();
      gluOrtho2D(-5,5,-5,5);
      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      glLoadIdentity();
//++ Farbe: rot fuer Schleifen
      glColor3d(1,0,0);
      glBegin(GL_POLYGON);
      glVertex2d(-5,-5);
      glVertex2d(-5,+5);
      glVertex2d(+5,+5);
      glVertex2d(+5,-5);
      glEnd();
      glMatrixMode(GL_PROJECTION);
      glPopMatrix();
      glMatrixMode(GL_MODELVIEW);
      glPopMatrix();

      glDisable(GL_STENCIL_TEST);
      glDepthMask(GL_TRUE);
      glEnable(GL_DEPTH_TEST);
   }
   else
   {
      glClear(GL_COLOR_BUFFER_BIT);
      glClear(GL_STENCIL_BUFFER_BIT);
   }
   glClear(GL_DEPTH_BUFFER_BIT);

   glEnable(GL_LINE_SMOOTH);
   glEnable(GL_POINT_SMOOTH);
   glEnable(GL_POLYGON_SMOOTH);

   glPushMatrix();
   glCallList(LIST_SPEZIAL);
   glCallList(LIST_SPEZIAL_SUB);
   glPopMatrix();
   glDisable(GL_LINE_SMOOTH);
   glDisable(GL_POINT_SMOOTH);
   glDisable(GL_POLYGON_SMOOTH);

   og_planets(planeten,elemente,jde,RADIUS,AEQUINOKTIUM);

   for(i=0;i<MAX_VERBINDUNGEN;i+=2)
   {
      if (verbindungen[i]!=-1)
      {
         og_strahl(pos_id(jde,verbindungen[i]),pos_id(jde,verbindungen[i+1]));
      }
   }
   glutSwapBuffers();
}

static void SetGlutSpecialFuncAll( void (*func)(int key, int x,int y) )
{
   glutSetWindow(GlutWinText);
   glutSpecialFunc((*func));
   glutSetWindow(GlutWinGrafik);
   glutSpecialFunc((*func));
   glutSetWindow(GlutWinUntertitel);
   glutSpecialFunc((*func));
}

static void SetGlutIdleFuncText( void (*func)(void) )
{
   glutSetWindow(GlutWinText);
   glutIdleFunc((*func));
}

static void init2( void )
{

   og_sterne(EIN,jde,AEQUINOKTIUM);

   glutSetWindow(GlutWinGrafik);
   LIST_SPEZIAL    =glGenLists(1);
   LIST_SPEZIAL_SUB=glGenLists(1);
   LIST_STENCIL    =glGenLists(1);
   glNewList(LIST_SPEZIAL,GL_COMPILE);
   glEndList();
   glNewList(LIST_SPEZIAL_SUB,GL_COMPILE);
   glEndList();
   glNewList(LIST_STENCIL,GL_COMPILE);
   glEndList();
   glutHideWindow();

   glutSetWindow(GlutWinUntertitel);
//   load_ImageWeiter("CropWeiter.rgb");
   glutHideWindow();

   glutSetWindow(GlutWinText);
   load_ImageText("CropTitel.rgb");
   glutShowWindow();
   glutPostRedisplay();

   SetGlutSpecialFuncAll(GlutTastatur);

   schoner=EIN;
   SetGlutIdleFuncText(GlutIdleSchoner);
}

static void init( void )
{
   int i;

   AEQUINOKTIUM=J2000;
   RADIUS=0.05;
   jde=J2000;
   schrittweite=1.0;
   start=AUS;
   stencil=AUS;
   for(i=0;i<MAX_VERBINDUNGEN;i++)
   {
      verbindungen[i]=-1;
   }

   i=0;
   planeten.on[i]=1; planeten.p[i++]=so;
   planeten.n=i;

   i=0;
   elemente.on[i]=1; i+=elemente_laden("merkur.ele",elemente.e[i]);
   elemente.on[i]=1; i+=elemente_laden("venus.ele",elemente.e[i]);
   elemente.on[i]=1; i+=elemente_laden("erde.ele",elemente.e[i]);
   elemente.on[i]=1; i+=elemente_laden("mars.ele",elemente.e[i]);
   elemente.on[i]=1; i+=elemente_laden("halley.ele",elemente.e[i]);
   elemente.on[i]=1; i+=elemente_laden("planet1.ele",elemente.e[i]);
   elemente.on[i]=1; i+=elemente_laden("planet2.ele",elemente.e[i]);
   elemente.n=i;

   kamera_modus=FIXOBJ;
   id_modus=0;
   id_kamdir=1;
   id_kampos=0;
   id_fixpos=0;
   id_fixdir=1;
}

int main( int argc, char **argv )
{
   int w=800,h=600,uh=100;

   init();
   glutInit(&argc, argv);

   glutInitWindowSize(w,h);
   glutInitWindowPosition(0,0);
   glutInitDisplayMode(GLUT_RGB|GLUT_SINGLE);
   GlutWinText=glutCreateWindow(rcsid);
   glutDisplayFunc(GlutDisplayText);
   glutReshapeFunc(GlutReshapeText);
   glutSetCursor(GLUT_CURSOR_NONE);

   glutInitDisplayMode(GLUT_RGB|GLUT_SINGLE);
   GlutWinUntertitel=glutCreateSubWindow(GlutWinText,0,h-uh,w,uh);
   glutDisplayFunc(GlutDisplayUntertitel);
   glutReshapeFunc(GlutReshapeUntertitel);
   glutSetCursor(GLUT_CURSOR_NONE);

   glutInitWindowPosition(0,0);
   glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_STENCIL|GLUT_DEPTH);
   GlutWinGrafik=glutCreateSubWindow(GlutWinText,0,0,w,h-uh);
   glutDisplayFunc(GlutDisplayGrafik);
   glutReshapeFunc(GlutReshapeGrafik);
   glutSetCursor(GLUT_CURSOR_NONE);

   init2();

   glutMainLoop();

   return 0;
}
