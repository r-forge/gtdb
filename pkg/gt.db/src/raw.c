/*

  Copyright (C) 2009, Perlegen Sciences, Inc.
  Copyright (C) 2010, 23andMe, Inc.
  
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

#if (R_VERSION < R_Version(2,1,1))
#define RAW(x) ((Rbyte *)INTEGER(x))
#endif

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("gt.db", String)
#else
#define _(String) (String)
#endif

#define isRaw(x) (TYPEOF(x) == RAWSXP)

SEXP raw_to_hex(SEXP a)
{
	char *str;
	SEXP ans, dim;
	int nr, nc, i, j;

	if (!isRaw(a))
		error(_("argument should be a raw vector"));
	PROTECT(dim = getAttrib(a, R_DimSymbol));
	if (isMatrix(a)) {
		nr = INTEGER(dim)[0];
		nc = INTEGER(dim)[1];
	} else {
		nr = LENGTH(a);
		nc = 1;
	}
	PROTECT(ans = allocVector(STRSXP, nc));
	str = (char *)R_alloc(2*nr+1, sizeof(char));
	str[2*nr] = '\0';
	for (j = 0; j < nc; j++) {
		for (i = 0; i < nr; i++) {
			sprintf(str+i*2, "%02X", RAW(a)[i+j*nr]);
		}
		SET_STRING_ELT(ans, j, mkChar(str));
	}
	UNPROTECT(2);
	return ans;
}

SEXP hex_to_raw(SEXP s)
{
	SEXP ans;
	const char *str;
	int i, j, a, b, len, nc;

	if (!isString(s))
		error(_("argument should be a character vector"));
	nc = LENGTH(s);
	if (nc) {
		str = CHAR(STRING_ELT(s, 0));
		len = strlen(str);
	} else {
		len = 0;
	}
	if (len & 1)
		error(_("string length should be a multiple of 2"));
	for (i = 1; i < nc; i++)
		if (strlen(CHAR(STRING_ELT(s, i))) != len)
			error(_("string length mismatch"));
	PROTECT(ans = allocMatrix(RAWSXP, len>>1, nc));
	for (j = 0; j < nc; j++) {
		str = CHAR(STRING_ELT(s, j));
		for (i = 0; i < len; i += 2) {
			a = toupper(str[i]); b = toupper(str[i+1]);
			a = (a >= 'A') ? a - 'A' + 10 : a - '0';
			b = (b >= 'A') ? b - 'A' + 10 : b - '0';
			RAW(ans)[j*len + i>>1] = (a<<4)+b;
		}
	}
	UNPROTECT(1);
	return ans;
}

R_CallMethodDef callMethods[] = {
	{"raw_to_hex", (DL_FUNC)&raw_to_hex, 1},
	{"hex_to_raw", (DL_FUNC)&hex_to_raw, 1},
	{NULL, NULL, 0}
};

void R_init_raw(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
