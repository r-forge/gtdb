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
#include <string.h>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("gt.db", String)
#else
#define _(String) (String)
#endif

/*
  There is one 4x4 mapping table for each supported binary operator
*/

static int map_nr[] =
{ 0, 0, 0, 0,
  0, 1, 1, 1,
  0, 1, 1, 1,
  0, 1, 1, 1 };

static int map_eq[] =
{ 0, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1 };

static int map_ne[] =
{ 0, 0, 0, 0,
  0, 0, 1, 1,
  0, 1, 0, 1,
  0, 1, 1, 0 };

static int map_s0[] =
{ 0, 0, 0, 0,
  0, 0, 0, 1,
  0, 0, 0, 0,
  0, 1, 0, 0 };

static int map_s1[] =
{ 0, 0, 0, 0,
  0, 0, 1, 0,
  0, 1, 0, 1,
  0, 0, 1, 0 };

static int map_sn[] =
{ 0, 0, 0, 0,
  0, 2, 1, 0,
  0, 1, 2, 1,
  0, 0, 1, 2 };

SEXP do_gt_dist(SEXP gt1, SEXP gt2, SEXP op, SEXP ag, SEXP na)
{
	SEXP ans;
	int *ia, nc1, nc2, i, m, n, *map = map_nr;
	const char *ops, *ags;
	char *g1, *g2;
	enum { SUM, MIN, MAX } agi = SUM;

	if (!isString(gt1) || (length(gt1)<1))
		error(_("'gt1' argument should be a character vector"));
	if (!isString(gt2) || (length(gt2)<1))
		error(_("'gt2' argument should be a character vector"));
	if (length(gt1) != length(gt2))
		error(_("'gt1' and 'gt2' should be the same length"));
	if (!isString(op) || (length(op)!=1))
		error(_("'op' argument should be a single character string"));
	PROTECT(na = coerceVector(na, INTSXP));
	if (!isInteger(na) || (length(na)!=1))
		error(_("'na' argument should be a scalar"));

	nc1 = length(STRING_ELT(gt1, 0));
	nc2 = length(STRING_ELT(gt2, 0));
	for (i = 1; i < length(gt1); i++) {
		if ((nc1 != length(STRING_ELT(gt1, i))) ||
			(nc2 != length(STRING_ELT(gt2, i))))
			error(_("'gt' strings are not all the same length"));
	}

	ops = CHAR(STRING_ELT(op, 0));
	if (strcmp(ops, "n") == 0) {
		map = map_nr;
	} else if ((strcmp(ops, "==") == 0) || 
			   (strcmp(ops, "ibs==2") == 0)) {
		map = map_eq;
	} else if (strcmp(ops, "!=") == 0) {
		map = map_ne;
	} else if (strcmp(ops, "ibs==0") == 0) {
		map = map_s0;
	} else if (strcmp(ops, "ibs==1") == 0) {
		map = map_s1;
	} else if (strcmp(ops, "ibs") == 0) {
		map = map_sn;
	} else {
		error(_("unknown operation"));
	}

	ags = CHAR(STRING_ELT(ag, 0));
	if (strcmp(ags, "sum") == 0) {
		agi = SUM;
	} else if (strcmp(ags, "min") == 0) {
		agi = MIN;
	} else if (strcmp(ags, "max") == 0) {
		agi = MAX;
	} else {
		error(_("unknown aggregation function"));
	}

	for (i = 0; i < 4; i++)
		map[i*4] = map[i] = INTEGER(na)[0];

	PROTECT(ans = allocMatrix(INTSXP, nc1, nc2));
	ia = INTEGER(ans);
	for (m = 0; m < nc1*nc2; m++, ia++)
		*ia = NA_INTEGER;

	g1 = (char *)R_alloc(nc1, sizeof(char));
	g2 = (char *)R_alloc(nc2, sizeof(char));
	for (i = 0; i < LENGTH(gt1); i++) {

		memcpy(g1, CHAR(STRING_ELT(gt1, i)), nc1);
		memcpy(g2, CHAR(STRING_ELT(gt2, i)), nc2);
		for (m = 0; m < nc1; m++) {
			switch (g1[m]) {
			case 'a': g1[m] = 1; break;
			case 'h': g1[m] = 2; break;
			case 'b': g1[m] = 3; break;
			default:  g1[m] = 0; break;
			}
		}
		for (m = 0; m < nc2; m++) {
			switch (g2[m]) {
			case 'a': g2[m] = 1; break;
			case 'h': g2[m] = 2; break;
			case 'b': g2[m] = 3; break;
			default:  g2[m] = 0; break;
			}
		}

		ia = INTEGER(ans);
		for (m = 0; m < nc2; m++) {
			int *gm = map + g2[m]*4, nai = NA_INTEGER;
			switch (agi) {
			case SUM:
				for (n = 0; n < nc1; n++, ia++) {
					int v = gm[(int)g1[n]], p = *ia;
					if (v != nai)
						*ia = (p != nai) ? p + v : v;
				}
				break;
			case MIN:
				for (n = 0; n < nc1; n++, ia++) {
					int v = gm[(int)g1[n]], p = *ia;
					if (v != nai)
						if ((p == nai) || (p > v)) *ia = v;
				}
				break;
			case MAX:
				for (n = 0; n < nc1; n++, ia++) {
					int v = gm[(int)g1[n]], p = *ia;
					if (v != nai)
						if ((p == nai) || (p < v)) *ia = v;
				}
				break;
			}
		}
	}

	UNPROTECT(2);
	return ans;
}


