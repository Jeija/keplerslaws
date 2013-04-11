/***************************************************************************
                          og_comm.cpp  -  description
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

static char rcsid[] = "$Id: og_comm.cpp,v 1.1 2000/04/16 08:37:55 tkramer Exp $";

#include <stdlib.h>
#include <string.h>

#include "astromat.h"
#include "definiti.h"
#include "physdat.h"
#include "orbit.h"
#include "position.h"
#include "tscheby.h"

#include "gv_comm.h"

// #include <dmalloc.h>

static const signed long C_L=10000L;                // Faktor bei long-Umwandlung
static const double VOR_MAG=1.0; // 0.3
static const double EXP_MAG=0.11; //0.095; // 0.11
static const double LIM_MAG=6.5; // 5.0

static int OG_LIST_OBJEKT;
static int Vega_SternZahl=0;
static VEGA_STERN *Vega_Stern=NULL;

void og_init( const PLANETEN& planet, const ELEMENTE& elem )
// Wir arbeiten in einem bereits bestehendem GLUT-Fenster
{
   GLfloat farbe_weiss[4] = {  1.0,1.0,1.0,1.0  };
   int n=10,m=10;

   // Clear
   glViewport(0, 0, glutGet((GLenum)GLUT_WINDOW_WIDTH), glutGet((GLenum)GLUT_WINDOW_HEIGHT));
   glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
   glShadeModel(GL_FLAT);
   glEnable(GL_DEPTH_TEST);
   glDisable(GL_DITHER);
   glDisable(GL_LIGHTING);
   // Projektion
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluPerspective(50.0,double(glutGet((GLenum)GLUT_WINDOW_WIDTH))/glutGet((GLenum)GLUT_WINDOW_HEIGHT),0.1,100);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   
   glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,  farbe_weiss);
   glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,  farbe_weiss);

   // Objekt
   OG_LIST_OBJEKT=glGenLists(1);
   glNewList(OG_LIST_OBJEKT,GL_COMPILE);
   glPushMatrix();
   {
      GLUquadricObj *q  = gluNewQuadric();
      gluSphere(q,1.0,n,m);
      gluDeleteQuadric(q);
   }
   glPopMatrix();
   glEndList();

/*
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glLightfv(GL_LIGHT0,GL_DIFFUSE,farbe_weiss);
   glLightfv(GL_LIGHT0,GL_AMBIENT,farbe_weiss);  
   glColor4fv(farbe_weiss);
   glLightModelfv(GL_LIGHT_MODEL_AMBIENT,farbe_weiss);
*/
}

void og_planets( const PLANETEN& planet,
                 const ELEMENTE& elem, 
                 double jd, double RADIUS, double AEQUINOKTIUM )
{
   vektor v;
   int i,j;
   
   GLdouble objx;
   GLdouble objy;
   GLdouble objz;
   GLdouble modelMatrix[16];
   GLdouble projMatrix[16];
   GLint viewport[4];
   GLdouble winx;
   GLdouble winy;
   GLdouble winz; 
   
   glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
   glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
   glGetIntegerv(GL_VIEWPORT,viewport);
   
   for(i=0;i<planet.n;i++)
   {  
      if (planet.on[i]==0)
      {
         continue;
      }
      if (planet.p[i]==mo)
      {
         v=m_praez_ekl(jd,AEQUINOKTIUM)*/*(*/planet_hi(er,jd);//+mond_hi(jd));
      }
      else
      {
         v=m_praez_ekl(jd,AEQUINOKTIUM)*planet_hi(planet.p[i],jd);
      }
      glPushMatrix();
      {
         glTranslated(v[0],v[1],v[2]);
	 {
	    if (planet.p[i]==so)
	    {
	       glColor3d(1,1,0);
  	       glScaled(2*RADIUS,2*RADIUS,2*RADIUS);
	       glCallList(OG_LIST_OBJEKT);
	    }
	    else
	    {
	       glColor3d(1,1,1);
  	       glScaled(RADIUS,RADIUS,RADIUS);
	       glCallList(OG_LIST_OBJEKT);
            }
	 }
      }
      glPopMatrix();
      objx=v[0];
      objy=v[1];
      objz=v[2];
      gluProject(objx,objy,objz,modelMatrix,projMatrix,viewport,&winx,&winy,&winz);
      winx-=double(strlen(elem.e[i].name))*glutBitmapWidth(GLUT_BITMAP_HELVETICA_12,'x')/2.0;
      winy+=14;
      gluUnProject(winx,winy,winz,modelMatrix,projMatrix,viewport,&objx,&objy,&objz);
      glRasterPos3d(objx,objy,objz);
      for(j=0;j<int(strlen(PLANET_NAME[planet.p[i]]));j++)
      {
	 glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,PLANET_NAME[planet.p[i]][j]);
      }
   }
   for(i=0;i<elem.n;i++)
   {  
      if (elem.on[i]==0)
      {
         continue;
      }
      vektor P,Q,vel;
      double d0;

      if (elem.e[i].e<1.0) d0=elem.e[i].t-elem.e[i].M/elem.e[i].n; else d0=elem.e[i].t0;
         gaussvek(elem.e[i].kl,elem.e[i].i,elem.e[i].om,P,Q);
      v=kepler(elem.e[i].e,elem.e[i].q,jd,d0,P,Q,vel);
      
      if (elem.e[i].a_mo!=des_datums)
         v=m_praez_ekl(elem.e[i].a_jd,AEQUINOKTIUM)*v;
      else
	 v=m_praez_ekl(jd,AEQUINOKTIUM)*v;
      
      glPushMatrix();
      {
         glTranslated(v[0],v[1],v[2]);
	 glColor3d(1,1,1);
	 glScaled(RADIUS,RADIUS,RADIUS);
	 glCallList(OG_LIST_OBJEKT);
      }
      glPopMatrix();
      objx=v[0];
      objy=v[1];
      objz=v[2];
      gluProject(objx,objy,objz,modelMatrix,projMatrix,viewport,&winx,&winy,&winz);
      winx-=double(strlen(elem.e[i].name))*glutBitmapWidth(GLUT_BITMAP_HELVETICA_12,'x')/2.0;
      winy+=14;
      gluUnProject(winx,winy,winz,modelMatrix,projMatrix,viewport,&objx,&objy,&objz);
      glRasterPos3d(objx,objy,objz);
      for(j=0;j<int(strlen(elem.e[i].name));j++)
      {
	 glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,elem.e[i].name[j]);
      }
   }
} 

void og_orbit( const ELEMENT& elem, double AEQUINOKTIUM )
{
   glBegin(GL_LINE_STRIP);
   for(int j=0;j<=360;j++)
   {
      vektor P,Q,v;
      double r,T,t;

      T=360.0/(elem.n*RAD2DEG);
      t=J2000+T*double(j)/360.0;

      gaussvek(elem.kl,elem.i,elem.om,P,Q);	 
      r=(elem.q*(1.0+elem.e))/(1.0+elem.e*cos(j*DEG2RAD));
      v=r*(cos(j*DEG2RAD)*P+sin(j*DEG2RAD)*Q);

      if (elem.a_mo!=des_datums)
	 v=m_praez_ekl(elem.a_jd,AEQUINOKTIUM)*v;
      else
	 v=m_praez_ekl(t,AEQUINOKTIUM)*v;
      glVertex3d(v[0],v[1],v[2]);
   }
   glEnd();
}

void og_apsiden( AUS_EIN a, const PLANETEN& planet, const ELEMENTE& elem, double RADIUS, double AEQUINOKTIUM )
{
#ifdef HUIBUH
   int i;
   char str[100];
   
   if(a==EIN)
   {
      for(i=0;i<elem.n;i++)
      {  
         if (elem.on[i]==0)
	 {
	    printf("(geometry apsiden.%s1 {})\n",elem.e[i].name);
            printf("(geometry apsiden.%s2 {})\n",elem.e[i].name);
            continue;
	 }
	 
	 vektor P,Q,v0,v1,v2,v3;
	 double r,T,t;
         T=360.0/(elem.e[i].n*RAD2DEG);
         
	 gaussvek(elem.e[i].kl,elem.e[i].i,elem.e[i].om,P,Q);	 
	 r=(elem.e[i].q*(1.0+elem.e[i].e))/(1.0+elem.e[i].e*cos(  0*DEG2RAD));
	 v0=r*(cos(  0*DEG2RAD)*P+sin(  0*DEG2RAD)*Q);
         r=(elem.e[i].q*(1.0+elem.e[i].e))/(1.0+elem.e[i].e*cos(180*DEG2RAD));
	 v1=r*(cos(180*DEG2RAD)*P+sin(180*DEG2RAD)*Q);
	 r=(elem.e[i].q*(1.0+elem.e[i].e))/(1.0+elem.e[i].e*cos( 90*DEG2RAD));
	 v2=r*(cos( 90*DEG2RAD)*P+sin( 90*DEG2RAD)*Q);
         r=(elem.e[i].q*(1.0+elem.e[i].e))/(1.0+elem.e[i].e*cos(270*DEG2RAD));
	 v3=r*(cos(270*DEG2RAD)*P+sin(270*DEG2RAD)*Q);

	 if (elem.e[i].a_mo!=des_datums)
	 {
            v0=m_praez_ekl(elem.e[i].a_jd,AEQUINOKTIUM)*v0;
            v1=m_praez_ekl(elem.e[i].a_jd,AEQUINOKTIUM)*v1;
            v2=m_praez_ekl(elem.e[i].a_jd,AEQUINOKTIUM)*v2;
            v3=m_praez_ekl(elem.e[i].a_jd,AEQUINOKTIUM)*v3;
	 }
	 else
	 {
	    t=J2000+T*double(  0)/360.0;
	    v0=m_praez_ekl(t,AEQUINOKTIUM)*v0;
            t=J2000+T*double(180)/360.0;
	    v1=m_praez_ekl(t,AEQUINOKTIUM)*v1;
	    t=J2000+T*double( 90)/360.0;
	    v2=m_praez_ekl(t,AEQUINOKTIUM)*v2;
            t=J2000+T*double(270)/360.0;
	    v3=m_praez_ekl(t,AEQUINOKTIUM)*v3;
          }
	  
         printf("(geometry apsiden.%s1 { VECT 1 2 0\t2\t0\t%.12f %.12f %.12f\t %.12f %.12f %.12f})",
         elem.e[i].name,v0[0],v0[1],v0[2],v1[0],v1[1],v1[2]);
         printf("(geometry apsiden.%s2 { VECT 1 2 0\t2\t0\t%.12f %.12f %.12f\t %.12f %.12f %.12f})",
         elem.e[i].name,v2[0],v2[1],v2[2],v3[0],v3[1],v3[2]);
      }
   }
   else
   {
      for(i=0;i<elem.n;i++)
      {  
	 printf("(geometry apsiden.%s1 {})\n",elem.e[i].name);
         printf("(geometry apsiden.%s2 {})\n",elem.e[i].name);
      }
   }
#endif
}

void og_focus( const ELEMENT& elem, double AEQUINOKTIUM )
{
   vektor P,Q,v_perihel,v_aphel,v_f1,v_f2,v_a;
   double r,T,t;

   T=360.0/(elem.n*RAD2DEG);
   gaussvek(elem.kl,elem.i,elem.om,P,Q);	 
   r=(elem.q*(1.0+elem.e))/(1.0+elem.e*cos(  0*DEG2RAD));
   v_perihel=r*(cos(  0*DEG2RAD)*P+sin(  0*DEG2RAD)*Q);
   r=(elem.q*(1.0+elem.e))/(1.0+elem.e*cos(180*DEG2RAD));
   v_aphel  =r*(cos(180*DEG2RAD)*P+sin(180*DEG2RAD)*Q);

   if (elem.a_mo!=des_datums)
   {
      v_aphel  =m_praez_ekl(elem.a_jd,AEQUINOKTIUM)*v_aphel;
      v_perihel=m_praez_ekl(elem.a_jd,AEQUINOKTIUM)*v_perihel;
   }
   else
   {
      t=J2000+T*double(  0)/360.0;
      v_perihel=m_praez_ekl(t,AEQUINOKTIUM)*v_perihel;
      t=J2000+T*double(180)/360.0;
      v_aphel  =m_praez_ekl(t,AEQUINOKTIUM)*v_aphel;
   }

   v_a=0.5*(v_perihel-v_aphel);
   v_f1=v_aphel+(1-elem.e)*v_a; 
   v_f2=v_aphel+(1+elem.e)*v_a;
   v_f1[2]+=0.01;
   v_f2[2]+=0.01;

   double d=0.15;

   glBegin(GL_LINES);
   glVertex3d(v_f1[0]+d,v_f1[1],v_f1[2]);
   glVertex3d(v_f1[0]-d,v_f1[1],v_f1[2]);
   glVertex3d(v_f1[0],v_f1[1]+d,v_f1[2]);
   glVertex3d(v_f1[0],v_f1[1]-d,v_f1[2]);
   glVertex3d(v_f2[0]+d,v_f2[1],v_f2[2]);
   glVertex3d(v_f2[0]-d,v_f2[1],v_f2[2]);
   glVertex3d(v_f2[0],v_f2[1]+d,v_f2[2]);
   glVertex3d(v_f2[0],v_f2[1]-d,v_f2[2]);
   glEnd();
}

void og_camera( KAMERA_MODUS mode, vektor target, vektor position, vektor fixobj, vektor fixdir )
{ 
   if (mode==FOLLOWOBJ) // Center: target, Spectator: position
   {
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      gluLookAt(
         position[0],position[1],position[2],
         target[0],target[1],target[2],
         0.0, 0.0, 1.0);

   }
   else if (mode==FIXOBJ) // Center: fixobj, Up: fixdir
   {
      double angle;
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
//++ Quadrant richtig ??
      angle=-acos((fixdir|vektor(1.0,0.0,0.0))/sqrt(fixdir|fixdir));
      glRotated(angle*RAD2DEG,0.0,0.0,1.0);
      glTranslated(-fixobj[0],-fixobj[1],-fixobj[2]);
      glTranslated(0,0,-5.0);
   }
}

#define MAX_KAT_ZEILE 300

/* Helligkeitsformel nach Sky & Telescope March 1990 p.311
*/
static double mag( double m )
{
   return VOR_MAG*pow(10.0,EXP_MAG*(LIM_MAG-m));
}

static double sdextract( const char *str, int x1, int x2 )
{
   char s[MAX_KAT_ZEILE];
   char *end;
   double x;

   strcpy(s,str);
   s[x2]=0;
   x=strtod(s+x1-1,&end);

   return x;
}

/* Gedacht als Vergleichsprozedur fuer qsort();
*/
static int vega_stern_vgl( VEGA_STERN *a, VEGA_STERN *b )
{
   double d=double((*a).mag-(*b).mag);

   return int(d/fabs(d));
}

/*
   "(/cds/V/53/) f2"
   lese Sterne aus Katalog in stern_sequ ein
*/
static VEGA_STERN *init_cdsV53katalog( const char *file, double jd, int *n, VEGA_STERN *stern_sequ )
{
   short i;
   char zeile[MAX_KAT_ZEILE];
   double h,m,s,d;
   FILE *Sternkatalog;

   stern_sequ=(VEGA_STERN *)calloc(1628,sizeof(VEGA_STERN));

   Sternkatalog=fopen(file,"rb");
   if (!Sternkatalog)
   {
      fprintf(stderr,"Kann Sternkatalog \"%s\"nicht <F6>ffnen!\nProgrammende.\n",file);
      exit(1);
   }

   for(i=0;i<1628;i++)
   {
      fgets(zeile,MAX_KAT_ZEILE-1,Sternkatalog);
      zeile[MAX_KAT_ZEILE-1]=0;
      
      h=sdextract(zeile,39,40);
      m=sdextract(zeile,42,43);
      s=sdextract(zeile,45,48);
      h=(h+m/60.0+s/3600.0)*H2RAD;

      d=sdextract(zeile,51,52);
      m=sdextract(zeile,54,55);
      s=sdextract(zeile,57,58);
      d=(d+m/60.0+s/3600.0)*DEG2RAD;
      if (zeile[50-1]=='-')
         d*=-1.0;

      m=sdextract(zeile,134,138);

      stern_sequ[i].lam=long(h*C_L);
      stern_sequ[i].bet=long(d*C_L);
      stern_sequ[i].mag=long(mag(m)*C_L);
      stern_sequ[i].nr=(int)sdextract(zeile,1,4);
   }
   fclose(Sternkatalog);

   *n=i;

   return stern_sequ;
}

/* Sterne: Aequator -> Ekliptik
   EINGABE:
   jd     : Julianisches Datum
   equ2ekl: Matrix equ -> ekl 
   n      : Anzahl der Sterne
   daten  : Sterne aequatorial J2000 (werden ueberschrieben!)
*/
static void stern_equ2ekl( double jd, const matrix& equ2ekl, int n, VEGA_STERN *daten )
{
   typedef int (*fcmp)(const void *, const void *);
   vektor v;
   sphaer s;
   short i;

   for(i=0;i<n;i++)
   {
      s.lambda=double(daten[i].lam)/C_L;
      s.beta  =double(daten[i].bet)/C_L;
      s.r     =1.0;
      v=vektor(s);
      s=sphaer(equ2ekl*v);
      daten[i].lam=long(s.lambda*C_L);
      daten[i].bet=long(s.beta  *C_L);
   }
   qsort(daten,n,sizeof(VEGA_STERN),(fcmp)vega_stern_vgl);
}

void og_sterne( AUS_EIN init, double JD, double AEQUINOKTIUM )
{   
   if (init==EIN)
   {
      if (Vega_Stern!=NULL)
      {
         free(Vega_Stern);
      }
      Vega_Stern=init_cdsV53katalog("f2",JD,&Vega_SternZahl,Vega_Stern);
      stern_equ2ekl(JD,m_equ2ekl(ekls_m(JD))*m_praez_equ(J2000,AEQUINOKTIUM),Vega_SternZahl,Vega_Stern);
   }
   else
   {
      vektor v;
      double r;
      int i;
      double r_alt;
      
      r_alt=floor(Vega_Stern[0].mag/C_L*10.0)/10.0;
      glPointSize(r_alt);
      glBegin(GL_POINTS);      
      for(i=0;i<Vega_SternZahl;i++)
      {
         v=vektor(sphaer(double(Vega_Stern[i].lam)/C_L,double(Vega_Stern[i].bet)/C_L,70.0));
         r=floor(Vega_Stern[i].mag/C_L*10.0)/10.0;
         if (r!=r_alt)
	 {
            r_alt=r;
            glEnd();      
            glPointSize(r);
            glBegin(GL_POINTS);      
	 }
	 glVertex3d(v[0],v[1],v[2]);
      }
      glEnd();      
   }
}

void og_polygon( int n, vektor v[] )
{
   int i;

   glBegin(GL_POLYGON);
   for(i=0;i<n;i++)
      glVertex3d(v[i][0],v[i][1],v[i][2]);
   glEnd();
}

void og_strahl( vektor start, vektor end )
{
   glBegin(GL_LINE_STRIP);   
   glVertex3d(start[0],start[1],start[2]);
   glVertex3d(end[0],end[1],end[2]);
   glEnd();
}

void og_erde( const vektor& v, double jdut, int list, int n, vektor *vz )
{
   double w,eps,d_lam,d_eps;

   eps=ekls_m(jdut);
   nutation_ekl(jdut,&d_lam,&d_eps);

   w=lmst(jdut,90.0*DEG2RAD)+d_lam*cos(eps);;
   glPushMatrix();
   glTranslated(v[0],v[1],v[2]);
   glRotated(180.0+w*RAD2DEG,0.0,0.0,1.0);
   if (n>0)
   {
      glDisable(GL_LIGHTING);
      glBegin(GL_LINE_STRIP);
      glColor3d(1,0,0);
      for(int i=0;i<n;i++)
      {
         glVertex3d(vz[i][0],vz[i][1],vz[i][2]);
      }
      glEnd();
      glEnable(GL_LIGHTING);
   }
   glCallList(list);
   glPopMatrix();
}

void og_mond( const vektor& v,  double jd, int list )
{
   vektor a;
   sphaer s;
   matrix Ry,Rz;
   double w,alpha,beta;

   rot_achse_null_mer(mo,jd,a,&w);
   s=sphaer(a);
   alpha=-s.beta;
   beta=s.lambda;
   
   glPushMatrix();
   glTranslated(v[0],v[1],v[2]);
   //++ keine Ahnung, ob das so stimmt...
   glRotated(alpha*RAD2DEG,0.0,0.0,1.0);
   glRotated(90.0-beta*RAD2DEG,0.0,1.0,0.0);
   glRotated(0.0+w*RAD2DEG,0.0,0.0,1.0);
   glCallList(list);
   glPopMatrix();
}

void og_kegel( const vektor& spitze, const vektor& basis, double basis_radius )
{
   double alpha,beta;
   vektor a=(basis-spitze);
   sphaer s=sphaer(a);
   vektor r=(a^vektor(0,0,1));

   alpha=PI+s.lambda;
   beta =-s.beta;

   glPushMatrix();
   // 1. Translation
   glTranslated(basis[0],basis[1],basis[2]);
   // 2. Rotation
   glRotated(alpha*RAD2DEG,0,0,1);
   glRotated(90.0-beta*RAD2DEG,0,1,0);
   glutSolidCone(basis_radius,sqrt(a|a),20,1); 
   glPopMatrix();
}
