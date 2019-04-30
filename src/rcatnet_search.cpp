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
 * rcatnet_search.cpp
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

#include "utils.h"
#include "rcatnet.h"
#include "rcatnet_search.h"

#define DISCRETE_SAMPLE
#include "dag_search.h"
#include "dag_search_dc.h"
#undef DISCRETE_SAMPLE
#define PROB_SAMPLE
#include "dag_search.h"
#include "dag_search_dc.h"
#undef PROB_SAMPLE

RCatnetSearchD::RCatnetSearchD() {
	m_pRorder = 0;
	m_pRorderInverse = 0;
	m_pSearchParams = 0;
}

RCatnetSearchD::~RCatnetSearchD() {
	if(m_pRorder)
		CATNET_FREE(m_pRorder);
	m_pRorder = 0;
	if(m_pRorderInverse)
		CATNET_FREE(m_pRorderInverse);
	m_pRorderInverse = 0;
	if(m_pSearchParams)
		delete m_pSearchParams;
	m_pSearchParams = 0;
}

SEXP RCatnetSearchD::estimate(SEXP rSamples, SEXP rPerturbations, SEXP rClasses, SEXP rClsdist, 
			SEXP rMaxParents, SEXP rParentSizes, SEXP rMaxComplexity, SEXP rOrder, SEXP rNodeCats, 
			SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rMatEdgeLiks, SEXP rUseCache, SEXP rEcho) {

	int i, j, k, len, bUseCache, maxParentSet, maxComplexity, numnets, inet, echo, klmode;
 	int *pRsamples, *pRperturbations, *pSamples, *pPerturbations, 
		hasClasses, *pRclasses, *pClasses, 
		**parentsPool, **fixedParentsPool, *pPool, *pParentSizes;
	double *matEdgeLiks, *pMatEdgeLiks;

	RCatnet rcatnet;
	SEXP dim, rnodecat, rparpool, cnetlist, cnetnode;

	if(!isMatrix(rSamples))
		error("Data is not a matrix");

	PROTECT(rMaxParents = AS_INTEGER(rMaxParents));
	maxParentSet = INTEGER_POINTER(rMaxParents)[0];
	UNPROTECT(1);

	PROTECT(rMaxComplexity = AS_INTEGER(rMaxComplexity));
	maxComplexity = INTEGER_POINTER(rMaxComplexity)[0];
	UNPROTECT(1);

	PROTECT(rUseCache = AS_LOGICAL(rUseCache));
	bUseCache = LOGICAL(rUseCache)[0];
	//Rprintf("bUseCache = %d\n", bUseCache);
	UNPROTECT(1);

	PROTECT(rEcho = AS_LOGICAL(rEcho));
	echo = LOGICAL(rEcho)[0];
	UNPROTECT(1);

	klmode = 0;
	PROTECT(rClsdist = AS_INTEGER(rClsdist));
	klmode = INTEGER_POINTER(rClsdist)[0];
	UNPROTECT(1);

	hasClasses = 0;
	if(!isNull(rClasses) && isInteger(rClasses))
		hasClasses = 1;

	PROTECT(rSamples = AS_INTEGER(rSamples));
	dim = GET_DIM(rSamples);
	m_numNodes = INTEGER(dim)[0];
	m_numSamples = INTEGER(dim)[1]; 

	if(m_pRorder)
		CATNET_FREE(m_pRorder);
	m_pRorder = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
	if(m_pRorderInverse)
		CATNET_FREE(m_pRorderInverse);
	m_pRorderInverse = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
	if (!m_pRorder || !m_pRorderInverse) {
		CATNET_MEM_ERR();
	}

	PROTECT(rOrder = AS_INTEGER(rOrder));
	if(length(rOrder) < m_numNodes) {
		warning("Invalid nodeOrder parameter - reset to default node order.");
		for(i = 0; i < m_numNodes; i++)
			m_pRorder[i] = i + 1;
	}
	else {
		memcpy(m_pRorder, INTEGER(rOrder), m_numNodes*sizeof(int));
	}
	for(i = 0; i < m_numNodes; i++) {
		if(m_pRorder[i] <= 0 || m_pRorder[i] > m_numNodes) {
			error("Invalid nodeOrder parameter");		
		}	
		else
			m_pRorderInverse[m_pRorder[i]-1] = i + 1;
	}
	UNPROTECT(1);

	if(m_pSearchParams)
		delete m_pSearchParams;
	m_pSearchParams = new SEARCH_PARAMETERS(
		m_numNodes, m_numSamples, 
		maxParentSet, maxComplexity, echo, 
		!isNull(rNodeCats), 
		!isNull(rParentSizes), !isNull(rPerturbations), 
		!isNull(rParentsPool), !isNull(rFixedParentsPool), 
		!isNull(rMatEdgeLiks), 0, 
		NULL, this, 0, 0, hasClasses, klmode);
	if (!m_pSearchParams) {
		CATNET_MEM_ERR();
	}

	if(!isNull(rParentSizes)) {
		pParentSizes = m_pSearchParams->m_pParentSizes;
		PROTECT(rParentSizes = AS_INTEGER(rParentSizes));
		if(length(rParentSizes) == m_numNodes) { 
			for(i = 0; i < m_numNodes; i++)
				pParentSizes[i] = INTEGER(rParentSizes)[m_pRorder[i] - 1];
		}
		UNPROTECT(1);
	}

	pSamples = (int*)m_pSearchParams->m_pSamples;
	pRsamples = INTEGER(rSamples);
	for(j = 0; j < m_numSamples; j++) {
		for(i = 0; i < m_numNodes; i++) {
			pSamples[j*m_numNodes + i] = pRsamples[j*m_numNodes + m_pRorder[i] - 1];
			if(R_IsNA(pSamples[j*m_numNodes + i]) || pSamples[j*m_numNodes + i] < 1) {
				pSamples[j*m_numNodes + i] = CATNET_NAN;
			}
		}
	}
	UNPROTECT(1); // rSamples

	pPerturbations = 0;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
		pPerturbations = m_pSearchParams->m_pPerturbations;
		pRperturbations = INTEGER_POINTER(rPerturbations);
		for(j = 0; j < m_numSamples; j++) {
			for(i = 0; i < m_numNodes; i++) {
				pPerturbations[j*m_numNodes + i] = pRperturbations[j*m_numNodes + m_pRorder[i] - 1];
			}
		}
		UNPROTECT(1);
	}

	if(hasClasses) {
		pClasses = (int*)m_pSearchParams->m_pClasses;
		PROTECT(rClasses = AS_INTEGER(rClasses));
		pRclasses = INTEGER(rClasses);
		if (pClasses && pRclasses)
			memcpy(pClasses, pRclasses, m_numSamples*sizeof(int));
		UNPROTECT(1); // rClasses
	}

	if(!isNull(rNodeCats)) {
		PROTECT(rNodeCats = AS_LIST(rNodeCats));
		for(i = 0; i < m_numNodes; i++) {
			rnodecat = AS_INTEGER(VECTOR_ELT(rNodeCats, (int)(m_pRorder[i] - 1)));
			len = length(rnodecat);
			if(isVector(rnodecat) && len > 0) {
				m_pSearchParams->m_pNodeNumCats[i] = len;
				m_pSearchParams->m_pNodeCats[i] = (int*)CATNET_MALLOC(len*sizeof(int));
				if (m_pSearchParams->m_pNodeCats[i]) {
					for(j = 0; j < len; j++)
						m_pSearchParams->m_pNodeCats[i][j] = INTEGER(rnodecat)[j];
				}
			}
		}
		UNPROTECT(1);
	}

	parentsPool = 0;
	if(!isNull(rParentsPool)) {
		PROTECT(rParentsPool = AS_LIST(rParentsPool));
		parentsPool = m_pSearchParams->m_parentsPool;
		for(i = 0; i < m_numNodes; i++) {
			rparpool = AS_INTEGER(VECTOR_ELT(rParentsPool, (int)(m_pRorder[i] - 1)));
			len = length(rparpool);
			if(isVector(rparpool) && len > 0 && len <= m_numNodes) {
				parentsPool[i] = (int*)CATNET_MALLOC((len+1)*sizeof(int));
				pPool = INTEGER(rparpool);
				if (parentsPool[i] && pPool) {
					for(j = 0; j < len; j++) {
						if(pPool[j] > 0 && pPool[j] <= m_numNodes) {
							for(k = 0; k < m_numNodes; k++)
								if(pPool[j] == m_pRorder[k])
									break;
							if(k < m_numNodes)
								parentsPool[i][j] = k;
							else
								parentsPool[i][j] = -1;
						}
					}
					parentsPool[i][len] = -1;
				}
				if(m_pSearchParams->m_maxParentsPool < len)
					m_pSearchParams->m_maxParentsPool = len;
			}
		}
		UNPROTECT(1);
	}

	fixedParentsPool = 0;
	if(!isNull(rFixedParentsPool)) {
		PROTECT(rFixedParentsPool = AS_LIST(rFixedParentsPool));
		fixedParentsPool = m_pSearchParams->m_fixedParentsPool;
		for(i = 0; i < m_numNodes; i++) {
			rparpool = AS_INTEGER(VECTOR_ELT(rFixedParentsPool, (int)(m_pRorder[i] - 1)));
			len = length(rparpool);
			if(isVector(rparpool) && len > 0 && len <= m_numNodes) {
				fixedParentsPool[i] = (int*)CATNET_MALLOC((len+1)*sizeof(int));
			 	if(maxParentSet < len)
			    		maxParentSet = len;
				pPool = INTEGER(rparpool);
				if (fixedParentsPool[i] && pPool) {
					for(j = 0; j < len; j++) {
						if(pPool[j] > 0 && pPool[j] <= m_numNodes) {
							for(k = 0; k < m_numNodes; k++)
								if(pPool[j] == m_pRorder[k])
									break;
							if(k < m_numNodes)
								fixedParentsPool[i][j] = k;
							else
								fixedParentsPool[i][j] = -1;
						}
					}
					fixedParentsPool[i][len] = -1;
				}
				if(m_pSearchParams->m_maxParentsPool < len)
					m_pSearchParams->m_maxParentsPool = len;
			}
		}
		UNPROTECT(1);
	}

	if(!isNull(rMatEdgeLiks) && m_pSearchParams->m_matEdgeLiks) {
		PROTECT(rMatEdgeLiks = AS_NUMERIC(rMatEdgeLiks));
		matEdgeLiks  = m_pSearchParams->m_matEdgeLiks;
		pMatEdgeLiks = REAL(rMatEdgeLiks);
		for(j = 0; j < m_numNodes; j++) {
			for(i = 0; i < m_numNodes; i++) {
				matEdgeLiks[j*m_numNodes + i] = pMatEdgeLiks[(m_pRorder[j] - 1)*m_numNodes + m_pRorder[i] - 1];
			}
		}
		UNPROTECT(1);
	}

	if(bUseCache)
		setCacheParams(m_numNodes, maxParentSet, m_pRorder, m_pRorderInverse);

	search(m_pSearchParams);

	if(m_pSearchParams)
		delete m_pSearchParams;
	m_pSearchParams = 0;

	if(!m_nCatnets || !m_pCatnets) {
		warning("No networks are found");
		return R_NilValue;
	}

	// create a R-list of catNetworks
	numnets = 0;
	for(i = 0; i < m_nCatnets; i++) {
		if(m_pCatnets[i]) {
			m_pCatnets[i]->setNodesOrder(m_pRorder);
			numnets++;
		}
	}

	PROTECT(cnetlist = allocVector(VECSXP, numnets));

	inet = 0;
	for(i = 0; i < m_nCatnets; i++) {
		if(!m_pCatnets[i])
			continue;

		rcatnet = *m_pCatnets[i];

		//rcatnet.setCategoryIndices(m_pNodeNumCats, m_pNodeCats);

		PROTECT(cnetnode = rcatnet.genRcatnet("catNetwork"));

		SET_VECTOR_ELT(cnetlist, inet, cnetnode);
		UNPROTECT(1);
		inet++;
	}

	UNPROTECT(1);

	if(m_pRorder)
		CATNET_FREE(m_pRorder);
	m_pRorder = 0;
	if(m_pRorderInverse)
		CATNET_FREE(m_pRorderInverse);
	m_pRorderInverse = 0;

	return cnetlist;
}

/*
	RDagSearch class members
*/

RDagSearch::RDagSearch() {
	m_pRorder = 0;
	m_pSearchParams = 0;
}

RDagSearch::~RDagSearch() {
	if(m_pRorder)
		CATNET_FREE(m_pRorder);
	m_pRorder = 0;
	if(m_pSearchParams)
		delete m_pSearchParams;
	m_pSearchParams = 0;
}

SEXP RDagSearch::estimate(SEXP rSamples, SEXP rPerturbations, SEXP rClasses, SEXP rClsdist, 
			SEXP rMaxParents, SEXP rParentSizes, SEXP rMaxComplexity, SEXP rOrder, SEXP rNodeCats, 
			SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rMatEdgeLiks, SEXP rUseCache, SEXP rEcho, int bIntSample = 0) {

	int i, j, k, len, maxParentSet, maxCategories, maxComplexity, bEqualCategories, node, echo, klmode;
 	int *pRperturbations, *pPerturbations, **parentsPool, **fixedParentsPool, *pPool, *pParentSizes;
	double *matEdgeLiks, *pMatEdgeLiks;
	SEXP dim, rnodecat, rparpool;

	int sampleline, *pNodeOffsets;
	int *pRsamples, *pSamples;
	double *pfRsamples, *pfSamples;
	DAG_LIST<double, int> *pDagList;
	int hasClasses, *pRclasses, *pClasses;

	if(!isMatrix(rSamples))
		error("Data is not a matrix");

	PROTECT(rMaxParents = AS_INTEGER(rMaxParents));
	maxParentSet = INTEGER_POINTER(rMaxParents)[0];
	UNPROTECT(1);

	PROTECT(rMaxComplexity = AS_INTEGER(rMaxComplexity));
	maxComplexity = INTEGER_POINTER(rMaxComplexity)[0];
	UNPROTECT(1);

	PROTECT(rEcho = AS_LOGICAL(rEcho));
	echo = LOGICAL(rEcho)[0];
	UNPROTECT(1);

	klmode = 0;
	PROTECT(rClsdist = AS_INTEGER(rClsdist));
	klmode = INTEGER_POINTER(rClsdist)[0];
	UNPROTECT(1);

	hasClasses = 0;
	if(!isNull(rClasses) && isInteger(rClasses))
		hasClasses = 1;

	sampleline = 0;
	if(bIntSample) {
		dim = GET_DIM(rSamples);
		m_numNodes = INTEGER(dim)[0];
		m_numSamples = INTEGER(dim)[1]; 
	}
	else {
		dim = GET_DIM(rSamples);
		sampleline = INTEGER(dim)[0];
		m_numSamples = INTEGER(dim)[1]; 
		if(isNull(rNodeCats)) 
			error("Node categories must be specified");
		m_numNodes = length(rNodeCats);
	}

	if(m_pRorder)
		CATNET_FREE(m_pRorder);
	m_pRorder = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
	if (!m_pRorder) {
		CATNET_MEM_ERR();
	}

	PROTECT(rOrder = AS_INTEGER(rOrder));
	if(length(rOrder) < m_numNodes) {
		warning("Invalid nodeOrder parameter - reset to default node order.");
		for(i = 0; i < m_numNodes; i++)
			m_pRorder[i] = i + 1;
	}
	else {
		memcpy(m_pRorder, INTEGER(rOrder), m_numNodes*sizeof(int));
	}
	UNPROTECT(1);

	if(m_pSearchParams)
		delete m_pSearchParams;
	m_pSearchParams = new SEARCH_PARAMETERS(
		m_numNodes, m_numSamples, 
		maxParentSet, maxComplexity, echo, 
		!isNull(rNodeCats), 
		!isNull(rParentSizes), !isNull(rPerturbations), 
		!isNull(rParentsPool), !isNull(rFixedParentsPool), 
		!isNull(rMatEdgeLiks), 0, 
		NULL, this, sampleline, 0, hasClasses, klmode);
	if (!m_pSearchParams) {
		CATNET_MEM_ERR();
	}

	pPerturbations = 0;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
		pPerturbations = m_pSearchParams->m_pPerturbations;
		pRperturbations = INTEGER_POINTER(rPerturbations);
		for(j = 0; j < m_numSamples; j++) {
			for(i = 0; i < m_numNodes; i++) {
				pPerturbations[j*m_numNodes + i] = pRperturbations[j*m_numNodes + m_pRorder[i] - 1];
			}
		}
		UNPROTECT(1);
	}

	if(hasClasses) {
		pClasses = (int*)m_pSearchParams->m_pClasses;
		PROTECT(rClasses = AS_INTEGER(rClasses));
		pRclasses = INTEGER(rClasses);
		memcpy(pClasses, pRclasses, m_numSamples*sizeof(int));
		UNPROTECT(1); // rClasses
	}

	parentsPool = 0;
	if(!isNull(rParentsPool)) {
		PROTECT(rParentsPool = AS_LIST(rParentsPool));
		parentsPool = m_pSearchParams->m_parentsPool;
		for(i = 0; i < m_numNodes; i++) {
			rparpool = AS_INTEGER(VECTOR_ELT(rParentsPool, (int)(m_pRorder[i] - 1)));
			len = length(rparpool);
			if(isVector(rparpool) && len > 0 && len <= m_numNodes) {
				parentsPool[i] = (int*)CATNET_MALLOC((len+1)*sizeof(int));
				pPool = INTEGER(rparpool);
				if (parentsPool[i] && pPool) {
					for(j = 0; j < len; j++) {
						if(pPool[j] > 0 && pPool[j] <= m_numNodes) {
							for(k = 0; k < m_numNodes; k++)
								if(pPool[j] == m_pRorder[k])
									break;
							if(k < m_numNodes)
								parentsPool[i][j] = k;
							else
								parentsPool[i][j] = -1;
						}
					}
					parentsPool[i][len] = -1;
				}
				if(m_pSearchParams->m_maxParentsPool < len)
					m_pSearchParams->m_maxParentsPool = len;
			}
		}
		UNPROTECT(1);
	}

	fixedParentsPool = 0;
	if(!isNull(rFixedParentsPool)) {
		PROTECT(rFixedParentsPool = AS_LIST(rFixedParentsPool));
		fixedParentsPool = m_pSearchParams->m_fixedParentsPool;
		for(i = 0; i < m_numNodes; i++) {
			rparpool = AS_INTEGER(VECTOR_ELT(rFixedParentsPool, (int)(m_pRorder[i] - 1)));
			len = length(rparpool);
			if(isVector(rparpool) && len > 0 && len <= m_numNodes) {
				fixedParentsPool[i] = (int*)CATNET_MALLOC((len+1)*sizeof(int));
			 	if(maxParentSet < len)
			    		maxParentSet = len;
				pPool = INTEGER(rparpool);
				if (fixedParentsPool[i] && pPool) {
					for(j = 0; j < len; j++) {
						if(pPool[j] > 0 && pPool[j] <= m_numNodes) {
							for(k = 0; k < m_numNodes; k++)
								if(pPool[j] == m_pRorder[k])
									break;
							if(k < m_numNodes)
								fixedParentsPool[i][j] = k;
							else
								fixedParentsPool[i][j] = -1;
						}
					}
					fixedParentsPool[i][len] = -1;
				}
				if(m_pSearchParams->m_maxParentsPool < len)
					m_pSearchParams->m_maxParentsPool = len;
			}
		}
		UNPROTECT(1);
	}

	if(!isNull(rMatEdgeLiks) && m_pSearchParams->m_matEdgeLiks) {
		PROTECT(rMatEdgeLiks = AS_NUMERIC(rMatEdgeLiks));
		matEdgeLiks = m_pSearchParams->m_matEdgeLiks;
		pMatEdgeLiks = REAL(rMatEdgeLiks);
		for(j = 0; j < m_numNodes; j++) {
			for(i = 0; i < m_numNodes; i++) {
				matEdgeLiks[j*m_numNodes + i] = pMatEdgeLiks[(m_pRorder[j] - 1)*m_numNodes + m_pRorder[i] - 1];
			}
		}
		UNPROTECT(1);
	}

	if(!isNull(rParentSizes)) {
		pParentSizes = m_pSearchParams->m_pParentSizes;
		PROTECT(rParentSizes = AS_INTEGER(rParentSizes));
		if(length(rParentSizes) == m_numNodes) { 
			for(i = 0; i < m_numNodes; i++)
				pParentSizes[i] = INTEGER(rParentSizes)[m_pRorder[i] - 1];
		}
		UNPROTECT(1);
	}

	pDagList = 0;

	if(bIntSample) {
		PROTECT(rSamples = AS_INTEGER(rSamples));
		pSamples = (int*)m_pSearchParams->m_pSamples;
		pRsamples = INTEGER(rSamples);
		for(j = 0; j < m_numSamples; j++) {
			for(i = 0; i < m_numNodes; i++) {
				pSamples[j*m_numNodes + i] = pRsamples[j*m_numNodes + m_pRorder[i] - 1];
				if(R_IsNA(pSamples[j*m_numNodes + i]) || pSamples[j*m_numNodes + i] < 1) {
					pSamples[j*m_numNodes + i] = CATNET_NAN;
				}
			}
		}
		UNPROTECT(1); // rSamples
		
		maxCategories = 0;
		if(!isNull(rNodeCats)) {
			PROTECT(rNodeCats = AS_LIST(rNodeCats));
			for(i = 0; i < m_numNodes; i++) {
				rnodecat = AS_INTEGER(VECTOR_ELT(rNodeCats, (int)(m_pRorder[i] - 1)));
				len = length(rnodecat);
				if(maxCategories < len)
					maxCategories = len;
				//if(maxCategories > 0 && maxCategories != len)
				//	CATNET_ERR("Nodes should have equal number of categories");
				if(isVector(rnodecat) && len > 0) {
					m_pSearchParams->m_pNodeNumCats[i] = len;
					m_pSearchParams->m_pNodeCats[i] = (int*)CATNET_MALLOC(len*sizeof(int));
					if (!m_pSearchParams->m_pNodeCats[i]) {
						CATNET_MEM_ERR();
					}
					for(j = 0; j < len; j++)
						m_pSearchParams->m_pNodeCats[i][j] = INTEGER(rnodecat)[j];
				}
			}
			UNPROTECT(1);
		}
		
		bEqualCategories = 1;
		for(i = 0; i < m_numNodes; i++) 
			if(i > 1 && m_pSearchParams->m_pNodeNumCats[i] != m_pSearchParams->m_pNodeNumCats[0])
				bEqualCategories = 0;

		if(bEqualCategories) { 

		switch(maxParentSet) {
		case 1: switch(maxCategories) {
			case 2: pDagList = new DAGD_SEARCH<double, int, int, 1, 2>; break;
			case 3: pDagList = new DAGD_SEARCH<double, int, int, 1, 3>; break;
			case 4: pDagList = new DAGD_SEARCH<double, int, int, 1, 4>; break;
			default: CATNET_NOTSUPP_ERR();break;
			}
		break;
		case 2: switch(maxCategories) {
			case 2: pDagList = new DAGD_SEARCH<double, int, int, 2, 2>; break;
			case 3: pDagList = new DAGD_SEARCH<double, int, int, 2, 3>; break;
			case 4: pDagList = new DAGD_SEARCH<double, int, int, 2, 4>; break;
			default: CATNET_NOTSUPP_ERR();break;
			}
		break;
		case 3: switch(maxCategories) {
			case 2: pDagList = new DAGD_SEARCH<double, int, int, 3, 2>; break;
			case 3: pDagList = new DAGD_SEARCH<double, int, int, 3, 3>; break;
			case 4: pDagList = new DAGD_SEARCH<double, int, int, 3, 4>; break;
			default: CATNET_NOTSUPP_ERR();break;
			}
		break;
		case 4: switch(maxCategories) {
			case 2: pDagList = new DAGD_SEARCH<double, int, int, 4, 2>; break;
			case 3: pDagList = new DAGD_SEARCH<double, int, int, 4, 3>; break;
			case 4: pDagList = new DAGD_SEARCH<double, int, int, 4, 4>; break;
			default: CATNET_NOTSUPP_ERR();break;
			}
		break;
		default: CATNET_NOTSUPP_ERR();break;
		}
		} /* bEqualCategories */
		else {
			switch(maxParentSet) {
			case 1: 
				pDagList = new DAGD_SEARCH_DC<double, int, int, 1>; break;
			case 2: 
				pDagList = new DAGD_SEARCH_DC<double, int, int, 2>; break;
			case 3: 
				pDagList = new DAGD_SEARCH_DC<double, int, int, 3>; break;
			case 4: 
				pDagList = new DAGD_SEARCH_DC<double, int, int, 4>; break;
			default: CATNET_NOTSUPP_ERR();break;
			}
		} /* !bEqualCategories */

	}
	else /* !bIntSample */ {
		pNodeOffsets = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
		if (!pNodeOffsets) {
			CATNET_MEM_ERR();
		}
		memset(pNodeOffsets, 0, m_numNodes*sizeof(int));

		maxCategories = 0;
		PROTECT(rNodeCats = AS_LIST(rNodeCats));
		for(i = 0; i < m_numNodes; i++) {
			//rnodecat = AS_INTEGER(VECTOR_ELT(rNodeCats, (int)(m_pRorder[i] - 1)));
			rnodecat = AS_INTEGER(VECTOR_ELT(rNodeCats, i));
			len = length(rnodecat);
			if(maxCategories < len)
				maxCategories = len;
			//if(maxCategories > 0 && maxCategories != len)
			//	CATNET_ERR("Nodes should have equal number of categories");
			pNodeOffsets[i] = len;
			if(i > 0)
				pNodeOffsets[i] = pNodeOffsets[i-1] + len;
			if(isVector(rnodecat) && len > 0) {
				m_pSearchParams->m_pNodeNumCats[i] = len;
				m_pSearchParams->m_pNodeCats[i] = (int*)CATNET_MALLOC(len*sizeof(int));
				if (m_pSearchParams->m_pNodeCats[i]) {
					for(j = 0; j < len; j++)
						m_pSearchParams->m_pNodeCats[i][j] = INTEGER(rnodecat)[j];
				}
			}
		}
		for(i = m_numNodes - 1; i > 0; i--) 
			pNodeOffsets[i] = pNodeOffsets[i-1];
		pNodeOffsets[0] = 0;
		UNPROTECT(1);

		PROTECT(rSamples = AS_NUMERIC(rSamples));
		pfSamples = (double*)m_pSearchParams->m_pSamples;
		pfRsamples = REAL(rSamples);
		int ii = 0;
		if (pfSamples && pfRsamples) {
			for(i = 0; i < m_numNodes; i++) {
				for(j = 0; j < m_numSamples; j++) {
					memcpy(pfSamples+j*sampleline + ii, 
						pfRsamples+j*sampleline + pNodeOffsets[m_pRorder[i] - 1], 
						m_pSearchParams->m_pNodeNumCats[i]*sizeof(double));
					if(R_IsNA(pfSamples[j*sampleline + ii]) || pfSamples[j*sampleline + ii] < 0) {
						pfSamples[j*sampleline + ii] = CATNET_NAN;
					}
				}
				ii += m_pSearchParams->m_pNodeNumCats[i];
			}
		}
		UNPROTECT(1); // rSamples

		CATNET_FREE(pNodeOffsets);
		pNodeOffsets = 0;

		bEqualCategories = 1;
		for(i = 0; i < m_numNodes; i++) 
			if(i > 1 && m_pSearchParams->m_pNodeNumCats[i] != m_pSearchParams->m_pNodeNumCats[0])
				bEqualCategories = 0;

		if(bEqualCategories) {

		switch(maxParentSet) {
		case 1: switch(maxCategories) {
			case 2: pDagList = new DAGP_SEARCH<double, int, 1, 2>; break;
			case 3: pDagList = new DAGP_SEARCH<double, int, 1, 3>; break;
			case 4: pDagList = new DAGP_SEARCH<double, int, 1, 4>; break;
			default: CATNET_NOTSUPP_ERR();break;
			}
		break;
		case 2: switch(maxCategories) {
			case 2: pDagList = new DAGP_SEARCH<double, int, 2, 2>; break;
			case 3: pDagList = new DAGP_SEARCH<double, int, 2, 3>; break;
			case 4: pDagList = new DAGP_SEARCH<double, int, 2, 4>; break;
			default: CATNET_NOTSUPP_ERR();break;
			}
		break;
		case 3: switch(maxCategories) {
			case 2: pDagList = new DAGP_SEARCH<double, int, 3, 2>; break;
			case 3: pDagList = new DAGP_SEARCH<double, int, 3, 3>; break;
			case 4: pDagList = new DAGP_SEARCH<double, int, 3, 4>; break;
			default: CATNET_NOTSUPP_ERR();break;
			}
		break;
		case 4: switch(maxCategories) {
			case 2: pDagList = new DAGP_SEARCH<double, int, 4, 2>; break;
			case 3: pDagList = new DAGP_SEARCH<double, int, 4, 3>; break;
			case 4: pDagList = new DAGP_SEARCH<double, int, 4, 4>; break;
			default: CATNET_NOTSUPP_ERR();break;
			}
		break;
		default: CATNET_NOTSUPP_ERR();break;
		}
		} /* bEqualCategories */
		else {
			switch(maxParentSet) {
			case 1: 
				pDagList = new DAGP_SEARCH_DC<double, int, 1>; break;
			case 2: 
				pDagList = new DAGP_SEARCH_DC<double, int, 2>; break;
			case 3: 
				pDagList = new DAGP_SEARCH_DC<double, int, 3>; break;
			case 4: 
				pDagList = new DAGP_SEARCH_DC<double, int, 4>; break;
			default: CATNET_NOTSUPP_ERR();break;
			}
		} /* !bEqualCategories */
	}

	if(!pDagList) 
		CATNET_MEM_ERR();

	pDagList->search(m_pSearchParams);

	if(m_pSearchParams)
		delete m_pSearchParams;
	m_pSearchParams = 0;

	if(!pDagList->m_dagPars || pDagList->m_numDags < 1) {
		warning("No networks are found");
		return R_NilValue;
	}

	int *pn;
	SEXP plist, pint, ppars, pLoglik, pComplx;
	SEXP daglist = PROTECT(NEW_OBJECT(MAKE_CLASS("dagEvaluate")));

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_numNodes;
	SET_SLOT(daglist, install("numnodes"), pint);
	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_numSamples;
	SET_SLOT(daglist, install("numsamples"), pint);
	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = maxCategories;
	SET_SLOT(daglist, install("maxcats"), pint);
	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = maxParentSet;
	SET_SLOT(daglist, install("maxpars"), pint);
	UNPROTECT(1);

	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	for(k = 0; k < m_numNodes; k++) {
		node = m_pRorder[k]-1;
		if(!pDagList->m_parSlots[k] || pDagList->m_numParSlots[k] <= 0) {
			SET_VECTOR_ELT(plist, node, R_NilValue);
			continue;
		}
		PROTECT(ppars = NEW_INTEGER(pDagList->m_numParSlots[k]/*maxParentSet*/*maxParentSet));
		pn = INTEGER_POINTER(ppars);
		for(j = 0; j < pDagList->m_numParSlots[k]/*maxParentSet*/; j++) {
			i = 0;
			while(i < maxParentSet && pDagList->m_parSlots[k][j*maxParentSet+i] >= 0) {
				pn[j*maxParentSet+i] = 
					m_pRorder[pDagList->m_parSlots[k][j*maxParentSet+i]];
				i++;
			}
			for(; i < maxParentSet; i++)
				pn[j*maxParentSet+i] = 0;
		}
		SET_VECTOR_ELT(plist, node, ppars);
		UNPROTECT(1);
	}
	SET_SLOT(daglist, install("parSlots"), plist);
	UNPROTECT(1);

	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	for(k = 0; k < m_numNodes; k++) {
		node = m_pRorder[k]-1;
		if(!pDagList->m_parLogliks[k] || pDagList->m_numParSlots[k] <= 0) {
			SET_VECTOR_ELT(plist, node, R_NilValue);
			continue;
		}
		PROTECT(ppars = NEW_NUMERIC(pDagList->m_numParSlots[k]));
		memcpy(NUMERIC_POINTER(ppars), pDagList->m_parLogliks[k], pDagList->m_numParSlots[k]*sizeof(double));
		SET_VECTOR_ELT(plist, node, ppars);
		UNPROTECT(1);
	}
	SET_SLOT(daglist, install("parLogliks"), plist);
	UNPROTECT(1);

	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	for(k = 0; k < m_numNodes; k++) {
		node = m_pRorder[k]-1;
		if(!pDagList->m_parComplx[k] || pDagList->m_numParSlots[k] <= 0) {
			SET_VECTOR_ELT(plist, node, R_NilValue);
			continue;
		}
		PROTECT(ppars = NEW_INTEGER(pDagList->m_numParSlots[k]));
		memcpy(INTEGER_POINTER(ppars), pDagList->m_parComplx[k], pDagList->m_numParSlots[k]*sizeof(int));
		SET_VECTOR_ELT(plist, node, ppars);
		UNPROTECT(1);
	}
	SET_SLOT(daglist, install("parComplx"), plist);
	UNPROTECT(1);

	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	for(k = 0; k < m_numNodes; k++) {
		node = m_pRorder[k]-1;
		if(!pDagList->m_parSampleSize[k] || pDagList->m_numParSlots[k] <= 0) {
			SET_VECTOR_ELT(plist, node, R_NilValue);
			continue;
		}
		PROTECT(ppars = NEW_INTEGER(pDagList->m_numParSlots[k]));
		memcpy(INTEGER_POINTER(ppars), pDagList->m_parSampleSize[k], pDagList->m_numParSlots[k]*sizeof(int));
		SET_VECTOR_ELT(plist, node, ppars);
		UNPROTECT(1);
	}
	SET_SLOT(daglist, install("parSampleSize"), plist);
	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = pDagList->m_numDags;
	SET_SLOT(daglist, install("numDags"), pint);
	UNPROTECT(1);

	PROTECT(plist =  allocVector(VECSXP, pDagList->m_numDags));
	PROTECT(pLoglik = NEW_NUMERIC(pDagList->m_numDags));
	PROTECT(pComplx = NEW_INTEGER(pDagList->m_numDags));
	DAG_PARS<double> *pDags = pDagList->m_dagPars;
	char *pParBuff = (char*)CATNET_MALLOC((m_numNodes+1)*sizeof(int));
	int  *pIntBuff =  (int*)CATNET_MALLOC((m_numNodes+1)*sizeof(int));
	int nParBuff;
	if (!pParBuff || !pIntBuff) {
		CATNET_MEM_ERR();
	}
	for(k = 0; k < pDagList->m_numDags && pDags; k++) {
		NUMERIC_POINTER(pLoglik)[k] = pDags->loglik;
		INTEGER_POINTER(pComplx)[k] = pDags->complx;
		if(pDags->numPars == 0) {
			SET_VECTOR_ELT(plist, k, R_NilValue);
			continue;
		}
		nParBuff = m_numNodes;
		if(pDags->compressNumPars(pIntBuff, pParBuff, nParBuff, m_pRorder) <= 0) {
			SET_VECTOR_ELT(plist, k, R_NilValue);
			continue;
		}
		nParBuff = 1 + (int)((nParBuff*sizeof(char))/sizeof(int));
		PROTECT(ppars = NEW_INTEGER(nParBuff));
		memcpy(INTEGER_POINTER(ppars), pParBuff, nParBuff*sizeof(int));
		SET_VECTOR_ELT(plist, k, ppars);
		UNPROTECT(1);
		pDags = pDags->next;
	}

	CATNET_FREE(pParBuff);
	CATNET_FREE(pIntBuff);
	SET_SLOT(daglist, install("numPars"), plist);
	SET_SLOT(daglist, install("loglik"), pLoglik);
	SET_SLOT(daglist, install("complx"), pComplx);
	UNPROTECT(3);

	UNPROTECT(1); // cnet

	delete pDagList;
	pDagList = 0;

	if(m_pRorder)
		CATNET_FREE(m_pRorder);
	m_pRorder = 0;

	return daglist;
}


/*
	RCatnetSearchP class members
*/

RCatnetSearchP::RCatnetSearchP() {
	m_pRorder = 0;
	m_pRorderInverse = 0;
	m_pSearchParams = 0;
}

RCatnetSearchP::~RCatnetSearchP() {
	if(m_pRorder)
		CATNET_FREE(m_pRorder);
	m_pRorder = 0;
	if(m_pRorderInverse)
		CATNET_FREE(m_pRorderInverse);
	m_pRorderInverse = 0;
	if(m_pSearchParams)
		delete m_pSearchParams;
	m_pSearchParams = 0;
}

SEXP RCatnetSearchP::estimate(SEXP rSamples, SEXP rPerturbations, SEXP rClasses, SEXP rClsdist, 
			SEXP rMaxParents, SEXP rParentSizes, SEXP rMaxComplexity, SEXP rOrder, SEXP rNodeCats, 
			SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rMatEdgeLiks, SEXP rUseCache, SEXP rEcho) {

	int i, ii, j, k, len, sampleline, bUseCache, maxParentSet, maxComplexity, numnets, inet, echo, klmode;
 	int *pRperturbations, *pPerturbations, *pNodeOffsets, 
		**parentsPool, **fixedParentsPool, *pPool, *pParentSizes, 
		hasClasses, *pRclasses, *pClasses;
	double *pRsamples, *pSamples, *matEdgeLiks, *pMatEdgeLiks;

	RCatnetP rcatnet;
	SEXP dim, rnodecat, rparpool, cnetlist, cnetnode;

	if(!isMatrix(rSamples))
		error("Data is not a matrix");
Rprintf("RCatnetSearchP\n");
	PROTECT(rMaxParents = AS_INTEGER(rMaxParents));
	maxParentSet = INTEGER_POINTER(rMaxParents)[0];
	UNPROTECT(1);

	PROTECT(rMaxComplexity = AS_INTEGER(rMaxComplexity));
	maxComplexity = INTEGER_POINTER(rMaxComplexity)[0];
	UNPROTECT(1);

	PROTECT(rUseCache = AS_LOGICAL(rUseCache));
	bUseCache = LOGICAL(rUseCache)[0];
	//Rprintf("bUseCache = %d\n", bUseCache);
	UNPROTECT(1);

	PROTECT(rEcho = AS_LOGICAL(rEcho));
	echo = LOGICAL(rEcho)[0];
	UNPROTECT(1);

	klmode = 0;
	PROTECT(rClsdist = AS_INTEGER(rClsdist));
	klmode = INTEGER_POINTER(rClsdist)[0];
	UNPROTECT(1);

	hasClasses = 0;
	if(!isNull(rClasses) && isInteger(rClasses))
		hasClasses = 1;

	dim = GET_DIM(rSamples);
	sampleline = INTEGER(dim)[0];
	m_numSamples = INTEGER(dim)[1]; 

	if(isNull(rNodeCats)) 
		error("Node categories must be specified");
	m_numNodes = length(rNodeCats);

	if(m_pSearchParams)
		delete m_pSearchParams;
	m_pSearchParams = new SEARCH_PARAMETERS(
		m_numNodes, m_numSamples, 
		maxParentSet, maxComplexity, echo, 
		!isNull(rNodeCats), 
		!isNull(rParentSizes), !isNull(rPerturbations), 
		!isNull(rParentsPool), !isNull(rFixedParentsPool), 
		!isNull(rMatEdgeLiks), 0, 
		NULL, this, sampleline, 0, hasClasses, klmode);
	if (!m_pSearchParams) {
		CATNET_MEM_ERR();
	}

	if(m_pRorder)
		CATNET_FREE(m_pRorder);
	m_pRorder = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
	if(m_pRorderInverse)
		CATNET_FREE(m_pRorderInverse);
	m_pRorderInverse = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
	if (!m_pRorder || !m_pRorderInverse) {
		CATNET_MEM_ERR();
	}

	PROTECT(rOrder = AS_INTEGER(rOrder));
	if(length(rOrder) < m_numNodes) {
		warning("Invalid nodeOrder parameter - reset to default node order.");
		for(i = 0; i < m_numNodes; i++)
			m_pRorder[i] = i + 1;
	}
	else {
		memcpy(m_pRorder, INTEGER(rOrder), m_numNodes*sizeof(int));
	}
	for(i = 0; i < m_numNodes; i++) {
		if(m_pRorder[i] <= 0 || m_pRorder[i] > m_numNodes) {
			error("Invalid nodeOrder parameter");		
		}	
		else
			m_pRorderInverse[m_pRorder[i]-1] = i + 1;
	}
	UNPROTECT(1);

	pNodeOffsets = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
	if (!pNodeOffsets) {
		CATNET_MEM_ERR();
	}
	memset(pNodeOffsets, 0, m_numNodes*sizeof(int));

	PROTECT(rNodeCats = AS_LIST(rNodeCats));
	for(i = 0; i < m_numNodes; i++) {
		rnodecat = AS_INTEGER(VECTOR_ELT(rNodeCats, i));
		len = length(rnodecat);
		pNodeOffsets[i] = len;
		if(i > 0)
			pNodeOffsets[i] = pNodeOffsets[i-1] + len;
		if(isVector(rnodecat) && len > 0) {
			m_pSearchParams->m_pNodeNumCats[i] = len;
			m_pSearchParams->m_pNodeCats[i] = (int*)CATNET_MALLOC(len*sizeof(int));
			if (m_pSearchParams->m_pNodeCats[i]) {
				for(j = 0; j < len; j++)
					m_pSearchParams->m_pNodeCats[i][j] = INTEGER(rnodecat)[j];
			}
		}
	}
	for(i = m_numNodes - 1; i > 0; i--) 
		pNodeOffsets[i] = pNodeOffsets[i-1];
	pNodeOffsets[0] = 0;
	UNPROTECT(1);

	if(!isNull(rParentSizes)) {
		pParentSizes = m_pSearchParams->m_pParentSizes;
		PROTECT(rParentSizes = AS_INTEGER(rParentSizes));
		if(length(rParentSizes) == m_numNodes) { 
			for(i = 0; i < m_numNodes; i++)
				pParentSizes[i] = INTEGER(rParentSizes)[m_pRorder[i] - 1];
		}
		UNPROTECT(1);
	}
	
	PROTECT(rSamples = AS_NUMERIC(rSamples));
	pSamples  = (double*)m_pSearchParams->m_pSamples;
	pRsamples = REAL(rSamples);
	if (pSamples && pRsamples) {
		ii = 0;
		for(i = 0; i < m_numNodes; i++) {
			for(j = 0; j < m_numSamples; j++) {
				memcpy(pSamples+j*sampleline + ii, 
					pRsamples+j*sampleline + pNodeOffsets[m_pRorder[i] - 1], 
					m_pSearchParams->m_pNodeNumCats[i]*sizeof(double));
				if(R_IsNA(pSamples[j*sampleline + ii]) || pSamples[j*sampleline + ii] < 0) {
					pSamples[j*sampleline + ii] = CATNET_NAN;
				}
			}
			ii += m_pSearchParams->m_pNodeNumCats[i];
		}
	}
	UNPROTECT(1); // rSamples

	CATNET_FREE(pNodeOffsets);
	pNodeOffsets = 0;

	pPerturbations = 0;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
		pPerturbations = m_pSearchParams->m_pPerturbations;
		pRperturbations = INTEGER_POINTER(rPerturbations);
		for(j = 0; j < m_numSamples; j++) {
			for(i = 0; i < m_numNodes; i++) {
				pPerturbations[j*m_numNodes + i] = pRperturbations[j*m_numNodes + m_pRorder[i] - 1];
			}
		}
		UNPROTECT(1);
	}

	if(hasClasses) {
		PROTECT(rClasses = AS_INTEGER(rClasses));
		pClasses = (int*)m_pSearchParams->m_pClasses;
		pRclasses = INTEGER(rClasses);
		if (pClasses && pRclasses)
			memcpy(pClasses, pRclasses, m_numSamples*sizeof(int));
		UNPROTECT(1); // rClasses
	}

	parentsPool = 0;
	if(!isNull(rParentsPool)) {
		PROTECT(rParentsPool = AS_LIST(rParentsPool));
		parentsPool = m_pSearchParams->m_parentsPool;
		for(i = 0; i < m_numNodes; i++) {
			rparpool = AS_INTEGER(VECTOR_ELT(rParentsPool, (int)(m_pRorder[i] - 1)));
			len = length(rparpool);
			if(isVector(rparpool) && len > 0 && len <= m_numNodes) {
				parentsPool[i] = (int*)CATNET_MALLOC((len+1)*sizeof(int));
				pPool = INTEGER(rparpool);
				if (parentsPool[i] && pPool) {
					for(j = 0; j < len; j++) {
						if(pPool[j] > 0 && pPool[j] <= m_numNodes) {
							for(k = 0; k < m_numNodes; k++)
								if(pPool[j] == m_pRorder[k])
									break;
							if(k < m_numNodes)
								parentsPool[i][j] = k;
							else
								parentsPool[i][j] = -1;
						}
					}
					parentsPool[i][len] = -1;
				}
				if(m_pSearchParams->m_maxParentsPool < len)
					m_pSearchParams->m_maxParentsPool = len;
			}
		}
		UNPROTECT(1);
	}

	fixedParentsPool = 0;
	if(!isNull(rFixedParentsPool)) {
		PROTECT(rFixedParentsPool = AS_LIST(rFixedParentsPool));
		fixedParentsPool = m_pSearchParams->m_fixedParentsPool;
		for(i = 0; i < m_numNodes; i++) {
			rparpool = AS_INTEGER(VECTOR_ELT(rFixedParentsPool, (int)(m_pRorder[i] - 1)));
			len = length(rparpool);
			if(isVector(rparpool) && len > 0 && len <= m_numNodes) {
				fixedParentsPool[i] = (int*)CATNET_MALLOC((len+1)*sizeof(int));
			 	if(maxParentSet < len)
			    		maxParentSet = len;
				pPool = INTEGER(rparpool);
				if (fixedParentsPool[i] && pPool) {
					for(j = 0; j < len; j++) {
						if(pPool[j] > 0 && pPool[j] <= m_numNodes) {
							for(k = 0; k < m_numNodes; k++)
								if(pPool[j] == m_pRorder[k])
									break;
							if(k < m_numNodes)
								fixedParentsPool[i][j] = k;
							else
								fixedParentsPool[i][j] = -1;
						}
					}
				}
				fixedParentsPool[i][len] = -1;
				if(m_pSearchParams->m_maxParentsPool < len)
					m_pSearchParams->m_maxParentsPool = len;
			}
		}
		UNPROTECT(1);
	}

	if(!isNull(rMatEdgeLiks) && m_pSearchParams->m_matEdgeLiks) {
		PROTECT(rMatEdgeLiks = AS_NUMERIC(rMatEdgeLiks));
		matEdgeLiks = m_pSearchParams->m_matEdgeLiks;
		pMatEdgeLiks = REAL(rMatEdgeLiks);
		for(j = 0; j < m_numNodes; j++) {
			for(i = 0; i < m_numNodes; i++) {
				matEdgeLiks[j*m_numNodes + i] = pMatEdgeLiks[(m_pRorder[j] - 1)*m_numNodes + m_pRorder[i] - 1];
			}
		}
		UNPROTECT(1);
	}

	if(bUseCache)
		setCacheParams(m_numNodes, maxParentSet, m_pRorder, m_pRorderInverse);

	search(m_pSearchParams);

	if(m_pSearchParams)
		delete m_pSearchParams;
	m_pSearchParams = 0;

	if(!m_nCatnets || !m_pCatnets) {
		warning("No networks are found");
		return R_NilValue;
	}

	// create a R-list of catNetworks
	numnets = 0;
	for(i = 0; i < m_nCatnets; i++) {
		if(m_pCatnets[i]) {
			m_pCatnets[i]->setNodesOrder(m_pRorder);
			numnets++;
		}
	}

	PROTECT(cnetlist = allocVector(VECSXP, numnets));

	inet = 0;
	for(i = 0; i < m_nCatnets; i++) {
		if(!m_pCatnets[i])
			continue;

		rcatnet = *m_pCatnets[i];

		PROTECT(cnetnode = rcatnet.genRcatnet("catNetwork"));

		SET_VECTOR_ELT(cnetlist, inet, cnetnode);
		UNPROTECT(1);
		inet++;
	}

	UNPROTECT(1);

	if(m_pRorder)
		CATNET_FREE(m_pRorder);
	m_pRorder = 0;
	if(m_pRorderInverse)
		CATNET_FREE(m_pRorderInverse);
	m_pRorderInverse = 0;
Rprintf("estimate exit");
	return cnetlist;
}


