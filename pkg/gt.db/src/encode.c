/*

  Copyright (C) 2009, 23andMe, Inc.
  
  Written by David A. Hinds <dhinds@sonic.net>
  
  This is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the license, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>
  
*/

#include <Rversion.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Arith.h>
#include <ctype.h>
#include <string.h>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("gt.db", String)
#else
#define _(String) (String)
#endif

/*---------------------------------------------------------------------*/

/*
  Parse a vector of space-delimited genotype strings, i.e. 'AA AC CC'
  as in HapMap files, into GT.DB notation, i.e. 'ahb', based on a
  vector of 'A/B' allele strings.
*/

SEXP do_encode_gt(SEXP sa, SEXP sg)
{
	int i, j, len;
	const char *a, *g;
	char *buf, a1, a2;
	SEXP ans;

	if (!isString(sa) || !isString(sg))
		error(_("arguments should be character vectors"));
	if (LENGTH(sa) != LENGTH(sg))
		error(_("argument lengths should be equal"));
	PROTECT(ans = allocVector(STRSXP, LENGTH(sa)));
	len = (LENGTH(STRING_ELT(sg,0))+1)/3;
	buf = (char *)R_alloc(len+1, sizeof(char));
	for (i = 0; i < LENGTH(sa); i++) {
		a = CHAR(STRING_ELT(sa,i));
		a1 = a[0]; a2 = a[2];
		g = CHAR(STRING_ELT(sg,i));
		for (j = 0; j < len && *g; j++) {
			buf[j] = '?';
			if (*g == a1) {
				if (g[1] == a1)
					buf[j] = 'a';
				else if (g[1] == a2)
					buf[j] = 'h';
			} else if (g[0] == a2) {
				if (g[1] == a2)
					buf[j] = 'b';
				else if (g[1] == a1)
					buf[j] = 'h';
			} else if (g[0] == 'N' && g[1] == 'N')
			       	buf[j] = 'n';
			if (buf[j] == '?')
				break;
			g += 2;
			while (*g == ' ') g++;
		}
		buf[j] = '\0';
		if (j < len)
			SET_STRING_ELT(ans, i, NA_STRING);
		else
			SET_STRING_ELT(ans, i, mkChar(buf));
	}
	UNPROTECT(1);
	return ans;
}
