/*

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
#include <string.h>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("gt.db", String)
#else
#define _(String) (String)
#endif

#define isRaw(x) (TYPEOF(x) == RAWSXP)

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
		if (a1 == '-') {
			/* codes for insertion/deletion polymorphisms */
			a1 = 'D';
			a2 = 'I';
		}
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

/*---------------------------------------------------------------------*/

/*
  Unpack a raw vector of packed 1/2/4-bit quantities into bytes
*/

SEXP do_unpack_bits(SEXP sr, SEXP sb)
{
	SEXP ans;
	const unsigned char *src;
	unsigned char *out, mask;
	int i, j, nb, len, bits;

	if (!isRaw(sr))
		error(_("first argument should be a raw vector"));

	bits = asInteger(sb);
	if ((bits != 1) && (bits != 2) && (bits != 4))
		error(_("invalid 'bits' (must be 1, 2, or 4)"));
	nb = 8/bits;
	mask = (1<<bits)-1;

	len = LENGTH(sr);
	PROTECT(ans = allocVector(RAWSXP, len*nb));
	src = RAW(sr);
	out = RAW(ans);

	for (i = 0; i < len; i++, src++) {
		unsigned char s = *src;
		for (j = 0; j < nb; j++, out++) {
			*out = s & mask;
			s >>= bits;
		}
	}
	UNPROTECT(1);
	return ans;
}

/*
  Unpack a raw vector of packed 1/2/4-bit quantities into bytes
*/

SEXP do_pack_bits(SEXP sr, SEXP sb)
{
	SEXP ans;
	const unsigned char *src;
	unsigned char *out, mask;
	int i, j, nb, len, bits, rem;

	if (!isRaw(sr))
		error(_("first argument should be a raw vector"));

	bits = asInteger(sb);
	if ((bits != 1) && (bits != 2) && (bits != 4))
		error(_("invalid 'bits' (must be 1, 2, or 4)"));
	nb = 8/bits;
	mask = (1<<bits)-1;

	len = LENGTH(sr) / nb;
	rem = LENGTH(sr) % nb;
	PROTECT(ans = allocVector(RAWSXP, len + (rem > 0)));
	src = RAW(sr);
	out = RAW(ans);

	for (i = 0; i < len; i++, out++) {
		unsigned char x = 0;
		for (j = 0; j < 8; j += nb, src++) {
			x |= (*src & mask) << j;
		}
		*out = x;
	}

	if (rem) {
		*out = 0;
		for (i = 0; i < rem; i++, src++) {
			*out |= (*src & mask) << (i*nb);
		}
	}

	UNPROTECT(1);
	return ans;
}
