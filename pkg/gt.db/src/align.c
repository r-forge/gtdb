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
   Global and ends-free sequence alignment by dynamic programming

   Needleman, S.B. & Wunsch C.D.  A general method applicable to the
   search for similarities in the amino acid sequence of two proteins.
   J Mol Biol 48: 443-453 (1970).
*/

SEXP do_nw_align(SEXP gt1, SEXP gt2, SEXP sparam, SEXP sends)
{
	int n1, n2, i, j, n, k;
	const char *s1, *s2;
	char *a1, *a2, *eq;
	double gap, pm, mm, **score;
	int ends[4], **link;
	SEXP ans, ssc, sal, end;
	enum gap { GAP_NO = 0, GAP_S1, GAP_S2 };

	if (!isString(gt1) || (length(gt1) != 1))
		error(_("'gt1' argument should be a character string"));
	if (!isString(gt2) || (length(gt2) != 1))
		error(_("'gt2' argument should be a character string"));
	if (!isReal(sparam) || LENGTH(sparam) != 3)
		error(_("'param' should a numeric vector of length 3"));
	if (!isLogical(sends) || (length(sends) != 4))
		error(_("'ends' argument should be a logical vector of length 4"));
	s1 = CHAR(STRING_ELT(gt1,0));
	s2 = CHAR(STRING_ELT(gt2,0));
	n1 = strlen(s1)+1;
	n2 = strlen(s2)+1;
	gap = REAL(sparam)[0];
	pm = REAL(sparam)[1];
	mm = REAL(sparam)[2];
	memcpy(ends, LOGICAL(sends), 4*sizeof(int));

	score = (double **)R_alloc(n1, sizeof(double *));
	score[0] = (double *)R_alloc(n2*n1, sizeof(double));
	link = (int **)R_alloc(n1, sizeof(int *));
	link[0] = (int *)R_alloc(n2*n1, sizeof(int));
	for (i = 1; i < n1; i++) {
		score[i] = score[0] + i*n2;
		link[i] = link[0] + i*n2;
	}

	/* initialize first rows of score, link matrices */
	if (ends[2]) {
		memset(score[0], 0, n2*sizeof(double));
	} else {
		score[0][0] = 0;
		for (j = 1; j < n2; j++)
			score[0][j] = j * gap;
	}
	link[0][0] = GAP_NO;
	for (j = 1; j < n2; j++)
		link[0][j] = GAP_S1;
	
	/* dynamic programming loops */
	for (i = 1; i < n1; i++) {
		score[i][0] = ends[0] ? 0 : (i * gap);
		link[i][0] = GAP_S2;
		for (j = 1; j < n2; j++) {
			double sc1, sc2, scm;
			sc1 = score[i][j-1] + gap;
			sc2 = score[i-1][j] + gap;
			scm = score[i-1][j-1] + (s1[i-1] == s2[j-1] ? pm : mm);
			if (scm >= sc1) {
				if (scm > sc2) {
					link[i][j] = GAP_NO;
					score[i][j] = scm;
				} else {
					link[i][j] = GAP_S2;
					score[i][j] = sc2;
				}
			} else {
				if (sc1 >= sc2) {
					link[i][j] = GAP_S1;
					score[i][j] = sc1;
				} else {
					link[i][j] = GAP_S2;
					score[i][j] = sc2;
				}
			}
		}
	}

	/* identify ends of best alignment */
	i = n1-1 ; j = n2-1;
	if (ends[3]) {
		for (n = n2-1; n > 0; n--) {
			if (score[i][n] >= score[i][j])
				j = n;
		}
	}
	if (ends[1]) {
		for (n = n1-1; n > 0; n--) {
			if (score[n][n2-1] >= score[i][j]) {
				i = n;
				j = n2-1;
			}
		}
	}
	PROTECT(end = allocVector(INTSXP,4));
	INTEGER(end)[1] = i;
	INTEGER(end)[3] = j;

	PROTECT(ssc = ScalarReal(score[i][j]));

	/* these will receive the final alignment */
	k = n1+n2;
	a1 = (char *)R_alloc(k+1, sizeof(char));
	a2 = (char *)R_alloc(k+1, sizeof(char));
	eq = (char *)R_alloc(k+1, sizeof(char));
	a1[k] = a2[k] = eq[k] = '\0';

	/* backtrack to extract best alignment */
	while ((i && j) || (i && !ends[0]) || (j && !ends[2])) {
		n = link[i][j];
		a1[--k] = (n == GAP_S1) ? '-' : toupper(s1[--i]);
		a2[k]   = (n == GAP_S2) ? '-' : toupper(s2[--j]);
		eq[k]   = (a1[k] == a2[k]) ? '|' : ' ';
	}
	INTEGER(end)[0] = i+1;
	INTEGER(end)[2] = j+1;

	PROTECT(sal = allocVector(STRSXP,3));
	SET_STRING_ELT(sal,0,mkChar(a1+k));
	SET_STRING_ELT(sal,1,mkChar(eq+k));
	SET_STRING_ELT(sal,2,mkChar(a2+k));

	PROTECT(ans = allocVector(VECSXP,3));
	SET_VECTOR_ELT(ans,0,ssc);
	SET_VECTOR_ELT(ans,1,end);
	SET_VECTOR_ELT(ans,2,sal);
	UNPROTECT(4);
	return ans;
}

/*
   Optimized version of alignment code that only computes the score
*/

SEXP do_nw_score(SEXP gt1, SEXP gt2, SEXP sparam, SEXP sends)
{
	int n1, n2, i, j, n;
	const char *s1, *s2;
	double gap, pm, mm, m, *score;
	int ends[4];

	if (!isString(gt1) || (length(gt1) != 1))
		error(_("'gt1' argument should be a character string"));
	if (!isString(gt2) || (length(gt2) != 1))
		error(_("'gt2' argument should be a character string"));
	if (!isReal(sparam) || LENGTH(sparam) != 3)
		error(_("'param' should a numeric vector of length 3"));
	if (!isLogical(sends) || (length(sends) != 4))
		error(_("'ends' argument should be a logical vector of length 4"));
	s1 = CHAR(STRING_ELT(gt1,0));
	s2 = CHAR(STRING_ELT(gt2,0));
	n1 = strlen(s1)+1;
	n2 = strlen(s2)+1;
	gap = REAL(sparam)[0];
	pm = REAL(sparam)[1];
	mm = REAL(sparam)[2];
	memcpy(ends, LOGICAL(sends), 4*sizeof(int));

	/* initialize first rows of score, link vectors */
	score = (double *)R_alloc(n2, sizeof(double));
	if (ends[2]) {
		memset(score, 0, n2*sizeof(double));
	} else {
		score[0] = 0;
		for (j = 1; j < n2; j++)
			score[j] = j * gap;
	}
	
	/* dynamic programming loops */
	m = 0.0;
	for (i = 1; i < n1; i++) {
		double sc0, sc1, sc2, scm;
		sc0 = score[0];
		sc1 = score[0] = ends[0] ? 0 : (i * gap);
		for (j = 1; j < n2; j++) {
			scm = sc0 + (s1[i-1] == s2[j-1] ? pm : mm);
			sc0 = score[j];
			sc2 = sc0 + gap;
			if (scm > sc1) {
				score[j] = (scm > sc2) ? scm : sc2;
			} else {
				score[j] = (sc1 > sc2) ? sc1 : sc2;
			}
			sc1 = score[j] + gap;
		}
		if (ends[1]) {
			if (score[n2-1] > m) m = score[n2-1];
		}
	}

	/* identify end of best alignment */
	if (ends[3]) {
		for (n = n2-1; n > 0; n--)
			if (score[n] > m)
				m = score[n];
	} else {
		if (score[n2-1] > m) m = score[n2-1];
	}
	return ScalarReal(m);
}
