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

#if !defined (_ENUMS_H)
#define _ENUMS_H

/* Schalter */
enum AUS_EIN           { AUS,                     /* Option aus */
                         EIN                      /* Option ein */
                       };

/* Aequinoktium */
enum AEQUINOKTIUM      { des_datums,              /* Des Datums */
                         j2000,                   /* J 2000 */
                         b1950,                   /* B 1950 */
                         anderes                  /* Sonstiges */
                       };

/* Koordinaten-Typ */
enum KOORDINATEN_TYP   { geometrisch,             /* Geometrisch */
                         astrometrisch,           /* Astrometrisch */
                         scheinbar                /* Scheinbar */
                       };

/* Planeten (auch fuer physische Ephemeriden) */
enum PLANET            { so,                      /* Sonne */
                         me,                      /* Merkur */
                         ve,                      /* Venus */
                         er,                      /* Erde */
                         ma,                      /* Mars */
                         ju,                      /* Jupiter */
                         sa,                      /* Saturn */
                         ur,                      /* Uranus */
                         ne,                      /* Neptun */
                         pl,                      /* Pluto */
                         mo,                      /* Mond */
       	                 ju_i,                    /* Jupiter Rotationssystem I   */
                         ju_ii,                   /* Jupiter Rotationssystem II  */
                         ju_iii,                  /* Jupiter Rotationssystem III */
                         sa_i,                    /* Saturn Rotationssystem I   */
                         sa_iii,                  /* Saturn Rotationssystem III */
	                 ju_a,                    /* Jupiter (Aequator) */
                         ju_p,                    /* Jupiter (Pol) */
                         sa_a,                    /* Saturn (Aequator) */
                         sa_p,                    /* Saturn (Pol) */
                         sa_ri,                   /* Saturn (Ring innen) */
                         sa_ra                    /* Saturn (Ring aussen) */
                       };

/* Geraet fuer die Ausgabe */
enum AUSGABEGERAET     { bild,                    /* Bildschirm */
                         drucker,                 /* Drucker */
                         datei                    /* Datei */
                       };

/* Modus der Ephemeride */
enum EPHEMERIDEN_MODUS { koor_eph,                /* Koordinaten bevorzugt */
                         adu_eph                  /* ADU-Zeiten bevorzugt */
                       };

/* Genauigkeit (Ausgabe) */
enum GENAUIGKEIT       { sph_lo,                  /* Sphaerisch gering */
                         sph_no,                  /* Sphaerisch normal */
                         sph_hi                   /* Sphaerisch hoch */
                       };

/* Eingrenzung (Intervall) */
enum VERGLEICH         { GG_KG,                   /* a <= x <= b */
                         G_KG,                    /* a <  x <= b */
                         G_K,                     /* a <  x <  b */
                         GG_K                     /* a <= x <  b */
                       };

/* Mode bei der Dateieroeffnung */
enum DATEI_MODUS       { text,                    /* Textmodus */
                         binaer                   /* Binaermodus */
                       };

/* Bahnform (Bahnelemente) */
enum BAHN_TYP          { b_ellipse,               /* Ellipse */
                         parabel,                 /* Parabel */
                         hyperbel                 /* Hyperbel */
                       };

/* Auf-, Durch- und Untergang */
enum GANG              { A,                       /* Ausgang */
                         D,                       /* Durchgang */
                         U                        /* Untergang */
                       };

/* Formatierung der Ausgabe */
enum WINKEL            { HH,                      /* Stunden */
                         DD,                      /* Grad */
                         MM,                      /* Minuten */
                         SS                       /* Sekunden */
                       };

/* Jupitermond - Ereignisse */
enum JUPITER_EREIGNIS  { kein_ereignis,           /* kein Ereignis */
		         bedeckung,               /* Bedeckung */
		         durchgang,               /* Durchgang */
		         verfinsterung,           /* Verfinsterung */
		         schatten,                /* Schattendurchgang */
		         bedeckung_verfinsterung, /* Bedeckung+Verfinsterung */
		         durchgang_schatten       /* Durchgang+Schattendurchgang */
	               };

enum FINSTER_ART { total,
	           ring,
	           ring_total,
	           part,
	           halb
                 };


/* Intervall-Art (Sichtbarkeitsdiagramm) */
enum INTERVALL_ART {  weltzeit,   /* immer zur gleichen Weltzeit    */
                      sternzeit,  /* immer zur gleichen Sternzeit   */
                      sonnenhoehe /* immer zur gleichen Sonnenhoehe */
                   };

enum ZENTER_TEXT { LEFT,RIGHT,TOP,BOTTOM,CENTER };

enum DRUCKER { epson_fx, hp_laser };


#endif /* !defined (_ENUMS_H) */
