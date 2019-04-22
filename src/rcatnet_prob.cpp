/*
 *  catnet : categorical Bayesian network inference
 *  Copyright (C) 2009--2019  Nikolay Balov
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.gnu.org/licenses/gpl-2.0.html
 */

/*
 * rcatnet_prob.cpp
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

#include "utils.h"
#include "rcatnet.h"

char *gen_prob_string(int node, SEXP parlist, int paridx, SEXP catlist, SEXP problist, char *str) {
	int j, npar;
	SEXP parprobs, pcats;
	char *newstr, *aux, *aux2, *aux3;

	//char ss[128];

	if(!str) {
		str = (char*)CATNET_MALLOC(1);
		str[0] = 0;
	}

	if(paridx >= length(parlist)) {
		pcats = VECTOR_ELT(catlist, node);
		newstr = (char*)CATNET_MALLOC(((strlen(str)+1+32)*length(pcats))*sizeof(char));
		if (!newstr)
			return str;
		newstr[0] = 0;
		for(j = 0; j < length(pcats); j++) {
			sprintf(newstr, "%s%s%s %f\n", newstr, str, CHAR(STRING_ELT(pcats, j)), NUMERIC_POINTER(problist)[j]);
		}

		CATNET_FREE(str);
		str = newstr;
		return str;
	}

	npar = INTEGER_POINTER(parlist)[paridx] - 1;
	pcats = VECTOR_ELT(catlist, npar);

	newstr = (char*)CATNET_MALLOC(sizeof(char));
	if (!newstr)
		return str;
	newstr[0] = 0;
	for(j = 0; j < length(pcats); j++) {
		parprobs = VECTOR_ELT(problist, j);

		aux = (char*)CATNET_MALLOC((strlen(str)+1+8)*sizeof(char));
		if (!aux)
			break;
		sprintf(aux, "%s%s", str, CHAR(STRING_ELT(pcats, j)));
		aux2 = gen_prob_string(node, parlist, paridx + 1, catlist, parprobs, aux);
		aux3 = (char*)CATNET_MALLOC((strlen(newstr)+strlen(aux2)+2)*sizeof(char));
		if (!aux2 || !aux3)
			break;
		sprintf(aux3, "%s%s", newstr, aux2);
		CATNET_FREE(newstr);
		newstr = aux3;

		CATNET_FREE(aux2);
	}
	CATNET_FREE(str);
	str = newstr;

	return str;
}

SEXP prob_string(SEXP rnodes, SEXP rparents, SEXP rcatlist, SEXP rproblist) {

	int node;
	SEXP nodepars, nodeproblist, pstr;
	char *str = NULL, *newstr, *aux;

	PROTECT(rnodes    = AS_LIST(rnodes));
	PROTECT(rparents  = AS_LIST(rparents));
	PROTECT(rcatlist  = AS_LIST(rcatlist));
	PROTECT(rproblist = AS_LIST(rproblist));

	for(node = 0; node < length(rnodes); node++) {
		nodepars = VECTOR_ELT(rparents, node);
		nodeproblist = VECTOR_ELT(rproblist, node);
		newstr = gen_prob_string(node, nodepars, 0, rcatlist, nodeproblist, NULL);
		if(str) {
			aux = (char*)CATNET_MALLOC((strlen(str)+strlen(newstr)+1+16)*sizeof(char));
			if (!aux)
				break;
			sprintf(aux, "%sNode [%d]:\n%s", str, node, newstr);
			CATNET_FREE(str);
			CATNET_FREE(newstr);
			str = aux;
		}
		else
			str = newstr;
	}

	PROTECT(pstr = allocVector(STRSXP, 1));
	SET_STRING_ELT(pstr, 0, mkChar(str));

	UNPROTECT(5);
	return pstr;
}


void gen_prob_vector(int node, SEXP parlist, int paridx, SEXP catlist, SEXP problist, double *&pvec, int &nvec) {
	int j, npar;
	SEXP parprobs, pcats;
	double *newvec;

	if(!pvec) {
		pvec = (double*)CATNET_MALLOC(sizeof(double));
		nvec = 0;
	}

	if(paridx >= length(parlist)) {
		pcats = VECTOR_ELT(catlist, node);
		if (length(problist) != length(pcats)) {
			Rprintf("gen_prob_vector: %d:  %d, %d\n", node, length(problist), length(pcats));
			error("Wrong probability table");
			return;
		}
		newvec = (double*)CATNET_MALLOC((nvec + length(pcats))*sizeof(double));
		if (!newvec)
			return;
		memcpy(newvec, pvec, nvec*sizeof(double));
		for(j = 0; j < length(pcats); j++) {
			newvec[nvec+j] = NUMERIC_POINTER(problist)[j];
		}
		CATNET_FREE(pvec);
		pvec = newvec;
		nvec += length(pcats);
		return;
	}

	npar = INTEGER_POINTER(parlist)[paridx] - 1;
	pcats = VECTOR_ELT(catlist, npar);
	if (length(problist) != length(pcats)) {
		Rprintf("gen_prob_vector: %d:  %d, %d\n", node, length(problist), length(pcats));
		error("Wrong probability table");
		return;
	}
	for(j = 0; j < length(pcats); j++) {
		parprobs = VECTOR_ELT(problist, j);
		gen_prob_vector(node, parlist, paridx + 1, catlist, parprobs, pvec, nvec);
	}
}

SEXP prob_vector(SEXP rnodes, SEXP rparents, SEXP rcatlist, SEXP rproblist) {

	int node, nvec;
	SEXP nodepars, nodeproblist, pstr;
	SEXP rvec;
	double *pvec, *prvec;

	PROTECT(rnodes = AS_LIST(rnodes));
	PROTECT(rparents = AS_LIST(rparents));
	PROTECT(rcatlist = AS_LIST(rcatlist));
	PROTECT(rproblist = AS_LIST(rproblist));

	PROTECT(pstr = allocVector(VECSXP, length(rnodes)));

	for(node = 0; node < length(rnodes); node++) {
		nodepars     = VECTOR_ELT(rparents, node);
		nodeproblist = VECTOR_ELT(rproblist, node);
		pvec = 0;
		nvec = 0;
		gen_prob_vector(node, nodepars, 0, rcatlist, nodeproblist, pvec, nvec);
		PROTECT(rvec = NEW_NUMERIC(nvec));
		prvec = NUMERIC_POINTER(rvec);
		if (prvec && pvec)
			memcpy(prvec, pvec, nvec*sizeof(double));
		CATNET_FREE(pvec);
		SET_VECTOR_ELT(pstr, node, rvec);
		UNPROTECT(1);
	}

	UNPROTECT(5);
	return pstr;
}


