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

#if !defined (_MONDDAT_H)
#define _MONDDAT_H

const double SYN_MON=29.530588853;  /* Synodischer Monat */
const double ANO_MON=27.55454988;   /* Anomalistischer Monat */
const double NEU_MON=2451550.09765; /* Erster Neumond 2000 */
const double PER_MON=2451534.6698;

double mond_perigaeum      ( double jd_start );
double mond_apogaeum       ( double jd_start );
double mond_mittlere_phase ( double k, double T );
double mond_neu            ( double jd_start );
double mond_erstes_viertel ( double jd_start );
double mond_voll           ( double jd_start );
double mond_letztes_viertel( double jd_start );
double mond_alter          ( double jd );
void mond_apsiden_phasen   ( double jd_start, double jde[6], short art[6] );

#endif /* !defined (_MONDDAT_H) */
