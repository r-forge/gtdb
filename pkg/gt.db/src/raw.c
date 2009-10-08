/*

  Copyright (C) 2009, Perlegen Sciences, Inc.
  
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
	SEXP ans;
	int i, len = LENGTH(a);

	if (!isRaw(a))
		error(_("argument should be a raw vector"));
	str = (char *)R_alloc(2*len+1, sizeof(char));
	for (i = 0; i < len; i++) {
		sprintf(str+i*2, "%02X", RAW(a)[i]);
	}
	ans = mkString(str);
	return ans;
}

SEXP hex_to_raw(SEXP s)
{
	SEXP ans;
	const char *str;
	int i, a, b, len;

	if (!isString(s) || (LENGTH(s) != 1))
		error(_("argument should be a character vector of length 1"));
	str = CHAR(STRING_ELT(s, 0));
	len = strlen(str);
	if (len & 1)
		error(_("string length should be a multiple of 2"));
	PROTECT(ans = allocVector(RAWSXP, len>>1));
	for (i = 0; i < len; i += 2) {
		a = toupper(str[i]); b = toupper(str[i+1]);
		a = (a >= 'A') ? a - 'A' + 10 : a - '0';
		b = (b >= 'A') ? b - 'A' + 10 : b - '0';
		RAW(ans)[i>>1] = (a<<4)+b;
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
