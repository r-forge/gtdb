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
#include <string.h>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("gt.db", String)
#else
#define _(String) (String)
#endif

/*---------------------------------------------------------------------*/

SEXP do_mask_str(SEXP ss, SEXP sm, SEXP sc)
{
	int n, len, i, j, k;
	SEXP ans;
	const char *m;
	char *buf, ch;

	if (!isString(ss) || !isString(sm) || !isString(sc))
		error(_("arguments should be character vectors"));
	if ((LENGTH(sm) != 1) || (STRING_ELT(sm,0) == NA_STRING) ||
		(LENGTH(sc) != 1) || (STRING_ELT(sm,0) == NA_STRING))
		error(_("invalid argument(s)"));
	ch = CHAR(STRING_ELT(sc,0))[0];
	m = CHAR(STRING_ELT(sm,0));
	len = strlen(m);
	n = LENGTH(ss);
	for (i = 0; i < n; i++) {
		if ((STRING_ELT(ss,i) != NA_STRING) &&
			(strlen(CHAR(STRING_ELT(ss,i))) != len))
			error(_("string length mismatch"));
	}

	PROTECT(ans = allocVector(STRSXP, n));
	buf = (char *)R_alloc(len+1, sizeof(char));
	for (i = 0; i < n; i++) {
		if (STRING_ELT(ss,i) == NA_STRING) {
			SET_STRING_ELT(ans, i, NA_STRING);
		} else {
			const char *s = CHAR(STRING_ELT(ss,i));
			if (ch) {
				for (j = 0; j < len; j++) {
					buf[j] = (m[j] == 'F') ? ch : s[j];
				}
			} else {
				for (j = k = 0; k < len; k++) {
					if (m[k] != 'F') buf[j++] = s[k];
				}
			}
			buf[j] = '\0';
			SET_STRING_ELT(ans, i, mkChar(buf));
		}
	}
	UNPROTECT(1);
	return ans;
}

/*---------------------------------------------------------------------*/

SEXP do_index_str(SEXP ss, SEXP si)
{
	int i, j, nc, ns, ni;
	char *buf;
    SEXP ans;

	if (!isString(ss))
		error(_("first argument should be a character vector"));
	if (!isInteger(si))
		error(_("second argument should be an integer vector"));

	ns = LENGTH(ss);
    ni = LENGTH(si);
	nc = strlen(CHAR(STRING_ELT(ss,0)));
	for (i = 0; i < ns; i++) {
		if ((STRING_ELT(ss,i) != NA_STRING) &&
			(strlen(CHAR(STRING_ELT(ss,i))) != nc))
			error(_("string length mismatch"));
	}
	for (i = 0; i < ni; i++) {
		if (INTEGER(si)[i] > nc)
			error(_("index out of range"));
	}

	PROTECT(ans = allocVector(STRSXP, ns));
	buf = (char *)R_alloc(ni+1, sizeof(char));
	buf[ni] = '\0';
	for (i = 0; i < ns; i++) {
		if (STRING_ELT(ss,i) == NA_STRING) {
			SET_STRING_ELT(ans, i, NA_STRING);
		} else {
			const char *s = CHAR(STRING_ELT(ss,i));
			for (j = 0; j < ni; j++) {
				buf[j] = s[INTEGER(si)[j]-1];
			}
			SET_STRING_ELT(ans, i, mkChar(buf));
		}
	}
	UNPROTECT(1);
	return ans;
}

/*---------------------------------------------------------------------*/

SEXP do_nsubstr(SEXP ss, SEXP sc)
{
	SEXP ans;
	int i, j, k, nc, ns;
	const char *s, *c;

	if (!isString(ss) || !isString(sc))
		error(_("arguments should be character vectors"));
	c = CHAR(STRING_ELT(sc, 0));
	nc = length(STRING_ELT(sc, 0));
	if (nc < 1)
		error(_("substr length should be > 0"));

	PROTECT(ans = allocVector(INTSXP, LENGTH(ss)));

	for (i = 0; i < LENGTH(ss); i++) {
		if (STRING_ELT(ss, i) != NA_STRING) {
			s = CHAR(STRING_ELT(ss, i));
			ns = length(STRING_ELT(ss, i));
			if (nc == 1)
				for (j = k = 0; k < ns; k++)
					j += (s[k] == *c);
			else
				for (j = k = 0; k < ns-nc+1; k++)
					j += !memcmp(s+k, c, nc);
			INTEGER(ans)[i] = j;
		} else {
			INTEGER(ans)[i] = NA_INTEGER;
		}
	}
	UNPROTECT(1);
	return ans;
}

/*---------------------------------------------------------------------*/

SEXP do_ch_table(SEXP ss, SEXP sch)
{
	int i, j, len, n[256], nc, ns, *ians;
	char ch[256];
	const char *s;
	SEXP ans, dimnames;

	if (!isString(ss) || !isString(sch))
		error(_("arguments should be character vectors"));

	/* get list of characters to be counted */
	nc = LENGTH(sch);
	for (j = 0; j < nc; j++) {
		if (STRING_ELT(sch,j) == NA_STRING)
			error(_("chars cannot include NA"));
		ch[j] = *CHAR(STRING_ELT(sch,j));
	}

	ns = LENGTH(ss);
	PROTECT(ans = allocMatrix(INTSXP, ns, nc));
	ians = INTEGER(ans);

	for (i = 0; i < ns; i++) {
		if (STRING_ELT(ss,i) == NA_STRING) {
			for (j = 0; j < nc; j++)
				ians[i+ns*j] = NA_INTEGER;
		} else {
			s = CHAR(STRING_ELT(ss,i));
			len = strlen(s);
			memset(n, 0, sizeof n);
			for (j = 0; j < len; j++)
				n[(int)s[j]]++;
			for (j = 0; j < nc; j++)
				ians[i+ns*j] = n[(int)ch[j]];
		}
	}

	PROTECT(dimnames = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dimnames, 1, sch);
	setAttrib(ans, R_DimNamesSymbol, dimnames);

	UNPROTECT(2);
	return(ans);
}

/*---------------------------------------------------------------------*/

SEXP do_ch_table_2(SEXP ss1, SEXP ss2, SEXP sch)
{
	int i, i1, i2, j, len, nc, n1, n2, ns, *ians;
	int m[256], *n;
	char ch;
	const char *s1, *s2;
	SEXP ans, dim, dimnames;

	if (!isString(ss1) || !isString(ss2) || !isString(sch))
		error(_("arguments should be character vectors"));

	/* set up mapping from character code to index */
	nc = LENGTH(sch);
	memset(m, 0, sizeof m);
	for (j = 0; j < nc; j++) {
		if (STRING_ELT(sch,j) == NA_STRING)
			error(_("chars cannot include NA"));
		ch = *CHAR(STRING_ELT(sch,j));
		m[(int)ch] = j+1;
	}
	n = (int *)R_alloc((nc+1)*(nc+1), sizeof(int));

	n1 = LENGTH(ss1);
	n2 = LENGTH(ss2);
	ns = (n1 && n2) ? ((n1 > n2) ? n1 : n2) : 0;
	PROTECT(ans = allocVector(INTSXP, ns*nc*nc));
	ians = INTEGER(ans);

	for (i = i1 = i2 = 0; i < ns; i++) {
		if (STRING_ELT(ss1, i1) == NA_STRING ||
			STRING_ELT(ss2, i2) == NA_STRING) {
			for (j = 0; j < nc*nc; j++)
				ians[i+ns*j] = NA_INTEGER;
		} else {
			s1 = CHAR(STRING_ELT(ss1,i1));
			s2 = CHAR(STRING_ELT(ss2,i2));
			len = strlen(s1);
			if (len != strlen(s2))
				error(_("string length mismatch"));
			memset(n, 0, (nc+1)*(nc+1)*sizeof(int));
			for (j = 0; j < len; j++)
				n[m[(int)s1[j]] + (nc+1)*m[(int)s2[j]]]++;
			for (j = 0; j < nc*nc; j++)
				ians[i+ns*j] = n[j+nc+2+j/nc];
		}
		if (++i1 >= n1) i1 = 0;
		if (++i2 >= n2) i2 = 0;
	}
	if (i1 || i2)
		warning(_("longer object length is not a multiple"
				  " of shorter object length"));

	PROTECT(dim = allocVector(INTSXP, 3));
	INTEGER(dim)[0] = ns;
	INTEGER(dim)[1] = nc;
	INTEGER(dim)[2] = nc;
	setAttrib(ans, R_DimSymbol, dim);

	PROTECT(dimnames = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(dimnames, 1, sch);
	SET_VECTOR_ELT(dimnames, 2, sch);
	setAttrib(ans, R_DimNamesSymbol, dimnames);

	UNPROTECT(3);
	return(ans);
}
