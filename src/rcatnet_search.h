/*
 *  catnet : categorical Bayesian network inference
 *  Copyright (C) 2009--2011  Nikolay Balov
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
 * rcatnet_search.h
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

#ifndef RCATNET_SEARCH_H
#define RCATNET_SEARCH_H

#define DISCRETE_SAMPLE
#include "catnet_search.h"
#undef DISCRETE_SAMPLE
#define PROB_SAMPLE
#include "catnet_search.h"
#undef PROB_SAMPLE

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

class RCatnetSearchD: public CATNETD_SEARCH<double> {
protected:
	int *m_pRorder, *m_pRorderInverse;
public:
	SEARCH_PARAMETERS *m_pSearchParams;

public:
	RCatnetSearchD();
	~RCatnetSearchD();

	SEXP estimate(SEXP rSamples, SEXP rPerturbations, SEXP rClasses, SEXP rClsdist, 
		SEXP rMaxParents, SEXP rParentSizes, SEXP rMaxComplexity, SEXP rOrder, SEXP rNodeCats, 
		SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rMatEdgeLiks, SEXP rUseCache, SEXP rEcho);

};

class RCatnetSearchP: public CATNETP_SEARCH<double> {
protected:
	int *m_pRorder, *m_pRorderInverse;
public:
	SEARCH_PARAMETERS *m_pSearchParams;

public:
	RCatnetSearchP();
	~RCatnetSearchP();

	SEXP estimate(SEXP rSamples, SEXP rPerturbations, SEXP rClasses, SEXP rClsdist, 
		SEXP rMaxParents, SEXP rParentSizes, SEXP rMaxComplexity, SEXP rOrder, SEXP rNodeCats, 
		SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rMatEdgeLiks, SEXP rUseCache, SEXP rEcho);

};

class RDagSearch {
protected:
	int *m_pRorder;
	int m_numNodes, m_numSamples;
public:
	SEARCH_PARAMETERS *m_pSearchParams;

public:
	RDagSearch();
	~RDagSearch();

	SEXP estimate(SEXP rSamples, SEXP rPerturbations, SEXP rClasses, SEXP rClsdist, 
		SEXP rMaxParents, SEXP rParentSizes, SEXP rMaxComplexity, SEXP rOrder, SEXP rNodeCats, 
		SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rMatEdgeLiks, SEXP rUseCache, SEXP rEcho, int bIntSample);

};

#endif /* RCATNET_SEARCH_H */
