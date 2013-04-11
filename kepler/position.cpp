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

/****************************************************************************/
/* POSITION.CC  Letzte Veraenderungen am ##.##.1993                         */
/*--------------------------------------------------------------------------*/
/* Berechnung von Planetenpositionen und von Positionen aus Bahnelementen   */
/*--------------------------------------------------------------------------*/
/* Abkuerzungen (Literatur):                                                */
/*                                                                          */
/* GdE   : Montenbruck: Grundlagen der Ephemeridenberechnung                */
/* AmdPC : Montenbruck,Pfleger: Astronomie mit dem Personal Computer        */
/* PAwyPC: Duffet-Smith: Practical Astronomy with your Personal Computer    */
/* AAyyyy: Astronomical Almanac of the Year yyyy                            */
/* AsAl  : Meeus: Astronomical Algorithms                                   */
/*****************************************************************************/

#include <assert.h>
#include <math.h>

#include "astromat.h"
#include "orbit.h"

#include "definiti.h"

#include "position.h"

sphaer planet_pos_hi( short planet, double jde );
sphaer planet_pos_xephem( PLANET planet, double jde );
sphaer pluto_hi( double jde );
sphaer mond_genau( double jde );

vektor element_gekl_geom( const ELEMENT& e, double jde, GENAUIGKEIT g )
/* Position aus Bahnelementen geozentrisch ekliptikal geometrisch (FK5)
   fuer das Aequinoktium des Datums

   EINGABE:
   e      : Bahnelemente
   jde    : Julianisches Datum (Ephemeridenzeit)
   g      : sph_hi: planet_hi() sonst: planet_lo()

   RETURN : s.o.
*/
{
   vektor P,Q,v_hekl_g,k_hekl_g;
   double d0;

   if (e.e<1.0) d0=e.t-e.M/e.n; else d0=e.t0;
   gaussvek(e.kl,e.i,e.om,P,Q);

   if (g==sph_hi)
      k_hekl_g=kepler(e.e,e.q,jde,d0,P,Q,v_hekl_g)-vektor(planet_hi(er,jde));
   else
      k_hekl_g=kepler(e.e,e.q,jde,d0,P,Q,v_hekl_g)-vektor(planet_lo(er,jde));

   if (e.a_mo!=des_datums) return m_praez_ekl(e.a_jd,jde)*k_hekl_g;
   else return k_hekl_g;
}

vektor element_gequ_app( const ELEMENT& e, double jde, GENAUIGKEIT g )
/* Position aus Bahnelementen geozentrisch aequatorial scheinbar (FK5)
   d.h. (a) Lichtlaufzeitkorrektur
        (b) stellare Aberrationskorrektur
        (c) Nutationskorrektur

   EINGABE:
   e      : Bahnelemente
   jde    : Julianisches Datum (Ephemeridenzeit)
   g      : sph_hi: planet_hi() sonst: planet_lo()
   RETURN : s.o.
*/
{
   double eps,d_lam,d_eps,d0,tau;
   vektor P,Q,v_hekl_g,k_hekl_g,k_gekl_s;
   sphaer s;
   matrix pmekl;

   eps=ekls_m(jde);

   if (e.e<1.0) d0=e.t-e.M/e.n; else d0=e.t0;
   gaussvek(e.kl,e.i,e.om,P,Q);
   k_hekl_g=kepler(e.e,e.q,jde,d0,P,Q,v_hekl_g);
   if (e.a_mo!=des_datums)
   {
      pmekl=m_praez_ekl(e.a_jd,jde);
      k_hekl_g=pmekl*k_hekl_g;
      v_hekl_g=pmekl*v_hekl_g;
   }

   if (g==sph_hi)
   {
      k_gekl_s=k_hekl_g-vektor(planet_hi(er,jde));
      tau=sqrt(k_gekl_s|k_gekl_s)*LZAE;
      k_gekl_s=(k_hekl_g-tau*v_hekl_g)-vektor(planet_hi(er,jde-tau));
   }
   else
   {
      k_gekl_s=k_hekl_g-vektor(planet_lo(er,jde));
      tau=sqrt(k_gekl_s|k_gekl_s)*LZAE;
      k_gekl_s=(k_hekl_g-tau*v_hekl_g)-vektor(planet_lo(er,jde-tau));
   }

   s=sphaer(k_gekl_s);
   nutation_ekl(jde,&d_lam,&d_eps);
   s.lambda+=d_lam;

   return m_ekl2equ(eps+d_eps)*vektor(sphaer(s.lambda,s.beta,s.r));
}

void element_pos( const ELEMENT& e,
                  double jde,
                  KOORDINATEN_TYP koor_typ,
                  AEQUINOKTIUM aqu_mo,
                  double aqu_jd,
                  GENAUIGKEIT g,

	          vektor& k_hekl_g, vektor& k_gekl_g, vektor& k_gequ_g,
                  vektor& s_gekl_g, vektor& k_gekl_k, vektor& k_gequ_k,
	          double* delta, double* tau, matrix& matrix_ekl2equ
                )
/* Diverse Positionen aus Bahnelementen

   EINGABE       :
   e             : Bahnelemente des Objekts
   jde           : Julianisches Ephemeriden Datum (ET)
   koor_typ      : (enum) geometrisch, astrometrisch oder scheinbar
   aqu_mo        : Aequinoktium (enum)
   aqu_jd        : Aequinoktium (JD)
   g             : sph_hi: planet_hi() sonst: planet_lo()

   AUSGABE       :
   k_hekl_g      : Komet hekl geometrisch Aequinoktium des Datums
   k_gekl_g      : Komet gekl geometrisch Aequinoktium des Datums
   k_gequ_g      : Komet gequ geometrisch Aequinoktium des Datums
   s_gekl_g      : Sonne gekl geometrisch Aequinoktium des Datums
   k_gekl_k      : Komet gekl nach koor_typ korrigiert Aequinoktium aqu_jd
   k_gequ_k      : Komet gequ nach koor_typ korrigiert Aequinoktium aqu_jd
   delta         : geometrische Entfernung
   tau           : geometrischer Lichtlaufweg
   matrix_ekl2equ: Matrix ekl -> equ (jde)
*/
{
   double eps,d_lam,d_eps,d0;
   vektor erde,P,Q,v_hekl_g;
   matrix pmekl;

   if (g==sph_hi)
      erde=(vektor)planet_hi(er,jde); /* geometrisch */
   else
      erde=(vektor)planet_lo(er,jde); /* geometrisch */

   s_gekl_g=-erde;
   eps=ekls_m(jde);
   matrix_ekl2equ=m_ekl2equ(eps);

   if (e.e<1.0) d0=e.t-e.M/e.n; else d0=e.t0;
   gaussvek(e.kl,e.i,e.om,P,Q);
   k_hekl_g=kepler(e.e,e.q,jde,d0,P,Q,v_hekl_g);
   if (e.a_mo!=des_datums)
   {
      pmekl=m_praez_ekl(e.a_jd,jde);
      k_hekl_g=pmekl*k_hekl_g;
      v_hekl_g=pmekl*v_hekl_g;
   }

   k_gekl_g=k_hekl_g-erde;
   k_gequ_g=matrix_ekl2equ*k_gekl_g;
   (*delta)=sqrt(k_gekl_g|k_gekl_g); /* geometrischer Lichtlaufweg */
   (*tau)  =(*delta)*LZAE;           /* geometrische Lichtlaufzeit */

   switch (koor_typ)
   {  case geometrisch  : if (aqu_mo!=des_datums)
                          {  k_gekl_k=m_praez_ekl(jde,aqu_jd)*k_gekl_g;
                             k_gequ_k=m_praez_equ(jde,aqu_jd)*k_gequ_g;
                          }
                          else
                          {  k_gekl_k=k_gekl_g;
                             k_gequ_k=k_gequ_g;
                          }
                          return;
      case astrometrisch: k_gekl_k=(k_hekl_g-(*tau)*v_hekl_g)-erde;
                          k_gequ_k=matrix_ekl2equ*k_gekl_k;
                          if (aqu_mo!=des_datums)
                          {  k_gekl_k=m_praez_ekl(jde,aqu_jd)*k_gekl_k;
                             k_gequ_k=m_praez_equ(jde,aqu_jd)*k_gequ_k;
                          }
                          return;
      case scheinbar    : if (g==sph_hi)
                             k_gekl_k=(k_hekl_g-(*tau)*v_hekl_g)
                                      -vektor(planet_hi(er,jde-(*tau)));
                          else
                             k_gekl_k=(k_hekl_g-(*tau)*v_hekl_g)
                                      -vektor(planet_lo(er,jde-(*tau)));
                          {
                             sphaer s=sphaer(k_gekl_k);
                             nutation_ekl(jde,&d_lam,&d_eps);
                             s.lambda+=d_lam;
                             k_gekl_k=vektor(sphaer(s.lambda,s.beta,s.r));
                             k_gequ_k=m_ekl2equ(eps+d_eps)*k_gekl_k;
                          }
                          return;
   }
}

double element_helligkeit( const ELEMENT& element, double i, double sp, double ep )
/* Helligkeit eines Kometen bzw. Kleinplaneten nach AsAl p. 217

   EINGABE:
   element: Bahnelemente des Objekts
   i      : Phasenwinkel
   sp     : Entfernung Sonne - Objekt in AE
   ep     : Entfernung Erde - Objekt in AE

   RETURN : s.o.
*/
{
   if (element.komet==EIN)
      return element.mag0+5.0*log10(ep)+2.5*element.mag1*log10(sp);
   else
   {  double Phi1,Phi2;

      Phi1=exp(-3.33*pow(tan(i*0.5),0.63));
      Phi2=exp(-1.87*pow(tan(i*0.5),1.22));
      return element.mag0+5.0*log10(ep*sp)
                         -2.5*log10((1.0-element.mag1)*Phi1+element.mag1*Phi2);
   }
}

sphaer planet_lo( PLANET planet, double jde )
/* Geometrische heliozentrisch ekliptikale Planetenposition (sphaer)
   fuer das Aequinoktium des Datums

   EINGABE:
   planet : (enum) Planet
   jde    : Julianisches Ephemeridendatum

   RETURN : s.o.
*/
{
   switch (planet)
   {  case so: return sphaer(0.0,0.0,0.0);
      case me: return mer(jde);
      case ve: return ven(jde);
      case er: {
                  sphaer s=son(jde);
                  return sphaer(s.lambda+PI,-s.beta,s.r);
               }
      case ma: return mar(jde);
      case ju: return jup(jde);
      case sa: return sat(jde);
      case ur: return ura(jde);
      case ne: return nep(jde);
      case pl: return plu(jde);
      default: assert(0); break;
   }
   return sphaer(0.0,0.0,0.0);
}

sphaer mond_lo( double jde )
/* Geometrische geozentrisch ekliptikale Mondposition (sphaer)
   fuer das Aequinoktium des Datums;
   Entfernung in AE

   EINGABE:
   jde    : Julianisches Ephemeridendatum

   RETURN : s.o.
*/
{  sphaer s;
   s=mon(jde);

   return sphaer(s.lambda,s.beta,s.r*RHO_AE);
}

sphaer planet_hi( PLANET planet, double jde )
/* Heliozentrisch ekliptikale Planetenposition (FK5) geometrisch (sphaer)
   fuer das Aequinoktium des Datums

   ACHTUNG: Die Files PLANETS.NDX und PLANETS.DAT muessen im aktuellen
            Verzeichnis stehen !

   EINGABE:
   planet : Planet (enum)
   jde    : Julianisches Ephemeridendatum

   RETURN : s.o.
*/
{
   switch (planet)
   {  case so: return sphaer(0.0,0.0,0.0);
      case me:
      case ve:
      case er:
      case ma:
      case ju:
      case sa:
      case ur:
      case ne: return planet_pos_xephem(planet,jde);
      case pl: return pluto_hi(jde);
      default: assert(0); break;
   }
   return sphaer(0.0,0.0,0.0);
}

sphaer mond_hi( double jde )
/* Mond geozentrisch ekliptikal ASTROMETRISCH (sphaer)
   fuer das Aequinoktium des Datums (allerdings fast kein Unterschied zu
   den geometrischen Koordinaten)
   Entfernung in km

   EINGABE:
   jde    : Julianisches Ephemeridendatum

   RETURN : s.o.

   Entfernung in AE
*/
{
   sphaer s;

   s=mond_genau(jde);
   s.r=s.r*1000.0/AE_M;

   return s;
}

vektor planet_gekl_geom( PLANET planet, double jde, GENAUIGKEIT g )
/* Planeten, Sonne und Mond geozentrisch ekliptikal geometrisch (FK5)
   geozentrisch ekliptikal geometrisch (vektor)
   fuer das Aequinoktium des Datums

   EINGABE:
   planet : Planet (enum)
   jde    : Julianisches Datum (Ephemeridenzeit)
   g      : sph_hi: planet_hi() sonst: planet_lo()

   RETURN : s.o.
*/
{
   switch (planet)
   {  case so: {
                  sphaer s;
                  if (g==sph_hi)
                     s=planet_hi(er,jde);
                  else
                     s=planet_lo(er,jde);

                  s.lambda+=PI;
                  s.beta   =-s.beta;
                  return vektor(sphaer(s.lambda,s.beta,s.r));
               }
      case er: return vektor(0.0,0.0,0.0);
      case mo: if (g==sph_hi)
                  return vektor(mond_hi(jde));
               else
                  return vektor(mond_lo(jde));
      default: if (g==sph_hi)
                  return vektor(planet_hi(planet,jde)-vektor(planet_hi(er,jde)));
               else
                  return vektor(planet_lo(planet,jde)-vektor(planet_lo(er,jde)));
   }
}

vektor planet_gequ_app( PLANET planet, double jde, GENAUIGKEIT g )
/* Planeten, Sonne und Mond geozentrisch aequatorial scheinbar (FK5) (vektor)
   d.h. (a) Lichtlaufzeitkorrektur
        (b) stellare Aberrationskorrektur
        (c) Nutationskorrektur

   EINGABE:
   planet : Planet (enum)
   jde    : Julianisches Datum (Ephemeridenzeit)
   g      : sph_hi: planet_hi() sonst: planet_lo()

   RETURN : s.o.
*/
{  double eps,tau,d_lam,d_eps;
   vektor p_gekl;

   switch (planet)
   {  case so: {
                  sphaer s;
                  if (g==sph_hi)
                     s=planet_hi(er,jde);
                  else
                     s=planet_lo(er,jde);

                  s.lambda+=(PI-20.4898*B2RAD/s.r);
                  s.beta   =-s.beta;
                  p_gekl=vektor(sphaer(s.lambda,s.beta,s.r));
               }
               break;
      case er: return vektor(0.0,0.0,0.0);
      case mo: if (g==sph_hi)
                  p_gekl=vektor(mond_hi(jde));
               else
                  p_gekl=vektor(mond_lo(jde));
               break;
      default: if (g==sph_hi)
               {  p_gekl=vektor(planet_hi(planet,jde)-vektor(planet_hi(er,jde)));
                  tau=sqrt(p_gekl|p_gekl)*LZAE;
                  p_gekl=vektor(planet_hi(planet,jde-tau))
                         -vektor(planet_hi(er,jde-tau));
               }
               else
               {  p_gekl=vektor(planet_lo(planet,jde)-vektor(planet_lo(er,jde)));
                  tau=sqrt(p_gekl|p_gekl)*LZAE;
                  p_gekl=vektor(planet_lo(planet,jde-tau))
                         -vektor(planet_lo(er,jde-tau));
               }
               break;
   }

   eps=ekls_m(jde);
   nutation_ekl(jde,&d_lam,&d_eps);

   {
      sphaer s=sphaer(p_gekl);
      s.lambda+=d_lam;
      p_gekl=vektor(sphaer(s.lambda,s.beta,s.r));
   }

   return m_ekl2equ(eps+d_eps)*p_gekl;
}

void planet_pos( PLANET planet,
                 double jde,
                 KOORDINATEN_TYP korr_typ,
                 AEQUINOKTIUM aqu_mo,
                 double aqu_jd,
                 GENAUIGKEIT g,

	         vektor& p_hekl_g, vektor& p_gekl_g, vektor& p_gequ_g,
                 vektor& s_gekl_g, vektor& p_gekl_k, vektor& p_gequ_k,
	         double* delta, double* tau, matrix& matrix_ekl2equ
               )
/* Planeten-, Sonne oder Mondposition (FK5)
   Beim Mond werden astrometrische und geometrische Koordinaten gleichgesetzt!

   EINGABE       :
   planet        : Planet
   jde           : Julianisches Ephemeriden Datum (ET)
   korr_typ      : (enum) geometrisch, astrometrisch oder scheinbar
   aqu_mo        : Aequinoktium (enum)
   aqu_jd        : Aequinoktium (JD)
   g             : sph_hi: planet_hi() sonst: planet_lo()

   AUSGABE       :
   p_hekl_g      : Planet hekl geometrisch Aequinoktium des Datums
   p_gekl_g      : Planet gekl geometrisch Aequinoktium des Datums
   p_gequ_g      : Planet gequ geometrisch Aequinoktium des Datums
   s_gekl_g      : Sonne gekl geometrisch Aequinoktium des Datums
   p_gekl_k      : Planet gekl nach korr_typ korrigiert Aequinoktium aqu_jd
   p_gequ_k      : Planet gequ nach korr_typ korrigiert Aequinoktium aqu_jd
   delta         : geometrische Entfernung
   tau           : geometrischer Lichtlaufweg
   matrix_ekl2equ: Matrix ekl -> equ (jde)
*/
{  double eps,d_lam,d_eps;
   vektor erde;

   if (g==sph_hi)
      erde=(vektor)planet_hi(er,jde); /* geometrisch */
   else
      erde=(vektor)planet_lo(er,jde); /* geometrisch */

   s_gekl_g=-erde;
   switch (planet)
   {  case so: p_hekl_g=vektor(0.0,0.0,0.0);
               p_gekl_g=s_gekl_g;
               break;
      case er: p_hekl_g=erde;
               p_gekl_g=vektor(0.0,0.0,0.0);
               break;
      case mo: if (g==sph_hi)
                  p_gekl_g=vektor(mond_hi(jde));
               else
                  p_gekl_g=vektor(mond_lo(jde));

               p_hekl_g=p_gekl_g+erde;
               break;
      default: if (g==sph_hi)
                  p_hekl_g=vektor(planet_hi(planet,jde));
               else
                  p_hekl_g=vektor(planet_lo(planet,jde));

               p_gekl_g=p_hekl_g-erde;
               break;
   }

   eps=ekls_m(jde);
   matrix_ekl2equ=m_ekl2equ(eps);
   p_gequ_g=matrix_ekl2equ*p_gekl_g;

   (*delta)=sqrt(p_gekl_g|p_gekl_g); /* geometrischer Lichtlaufweg */
   (*tau)  =(*delta)*LZAE;           /* geometrische Lichtlaufzeit */

   switch (korr_typ)
   {  case geometrisch  : if (aqu_mo!=des_datums)
                          {  p_gekl_k=m_praez_ekl(jde,aqu_jd)*p_gekl_g;
                             p_gequ_k=m_praez_equ(jde,aqu_jd)*p_gequ_g;
                          }
                          else
                          {  p_gekl_k=p_gekl_g;
                             p_gequ_k=p_gequ_g;
                          }
                          return;
      case astrometrisch: switch (planet)
                          {  case er: p_gekl_k=vektor(0.0,0.0,0.0);
                                      p_gequ_k=vektor(0.0,0.0,0.0);
                                      return;
                             case mo: /* Mond ist schon astrometrisch */
                                      p_gekl_k=p_gekl_g; break;
                             default: if (g==sph_hi)
                                         p_gekl_k=vektor(planet_hi(planet,jde-(*tau)))-erde;
                                      else
                                         p_gekl_k=vektor(planet_lo(planet,jde-(*tau)))-erde;
                                      break;
                          }
                          p_gequ_k=matrix_ekl2equ*p_gekl_k;
                          if (aqu_mo!=des_datums)
                          {  p_gekl_k=m_praez_ekl(jde,aqu_jd)*p_gekl_k;
                             p_gequ_k=m_praez_equ(jde,aqu_jd)*p_gequ_k;
                          }
                          return;
      case scheinbar    : switch (planet)
                          {  case so:{
                                        sphaer s=sphaer(p_gekl_g);
                                        nutation_ekl(jde,&d_lam,&d_eps);
                                        s.lambda+=(-20.4898*B2RAD/s.r+d_lam);
                                        p_gekl_k=
                                            vektor(sphaer(s.lambda,s.beta,s.r));
                                      }
                                      break;
                             case er: p_gekl_k=vektor(0.0,0.0,0.0);
                                      p_gequ_k=vektor(0.0,0.0,0.0);
                                      return;
                             case mo: {
                                         sphaer s=sphaer(p_gekl_g);
                                         nutation_ekl(jde,&d_lam,&d_eps);
                                         s.lambda+=d_lam;
                                         p_gekl_k=
                                            vektor(sphaer(s.lambda,s.beta,s.r));
                                      }
                                      break;
                             default: {
                                         sphaer s;
                                         if (g==sph_hi)
                                            p_gekl_k=vektor(planet_hi(planet,jde-(*tau))-vektor(planet_hi(er,jde-(*tau))));
                                         else
                                            p_gekl_k=vektor(planet_lo(planet,jde-(*tau))-vektor(planet_lo(er,jde-(*tau))));

                                         s=sphaer(p_gekl_k);
                                         nutation_ekl(jde,&d_lam,&d_eps);
                                         s.lambda+=d_lam;
                                         p_gekl_k=
                                            vektor(sphaer(s.lambda,s.beta,s.r));
                                      }
                                      break;
                          }
                          p_gequ_k=m_ekl2equ(eps+d_eps)*p_gekl_k;
                          return;
   }
}
