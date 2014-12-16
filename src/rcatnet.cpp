/*
 *  catnet : categorical Bayesian network inference
 *  Copyright (C) 2009--2010  Nikolay Balov
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
 * rcatnet.cpp
 *
 *  Created on: Sep 21, 2011
 *      Author: nbalov
 */

#include "utils.h"
#include "rcatnet.h"

#ifdef RCATNET_CLASS

RCATNET_CLASS::RCATNET_CLASS(SEXP cnet) {

	double *pvec;
	int nnode, i, nvec, *pn;
	char const *pstr;

	SEXP rname, rnodes, rparents, rcatlist, rproblist, pf, nodepars, nodeproblist, pint, rnodeprob;

	if(!isS4(cnet))
		return;

	rname     = GET_SLOT(cnet, install("objectName"));
	rnodes    = GET_SLOT(cnet, install("nodes"));
	rparents  = GET_SLOT(cnet, install("pars"));
	rcatlist  = GET_SLOT(cnet, install("cats"));
	rproblist = GET_SLOT(cnet, install("probs"));

	if(rnodes == R_NilValue || rparents == R_NilValue || 
	   rcatlist == R_NilValue || rproblist == R_NilValue) {
		return;
	}

	PROTECT(rname     = AS_CHARACTER(rname));
	PROTECT(rnodes    = AS_LIST(rnodes));
	PROTECT(rparents  = AS_LIST(rparents));
	PROTECT(rcatlist  = AS_LIST(rcatlist));
	PROTECT(rproblist = AS_LIST(rproblist));

	pint       = GET_SLOT(cnet, install("numnodes"));
	m_numNodes = INTEGER_POINTER(pint)[0];

	pint = GET_SLOT(cnet, install("maxpars"));
	m_maxParents = INTEGER_POINTER(pint)[0];
	pint = GET_SLOT(cnet, install("maxcats"));
	m_maxCategories = INTEGER_POINTER(pint)[0];
	pint = GET_SLOT(cnet, install("complx"));
	m_complexity = INTEGER_POINTER(pint)[0];
	pint = GET_SLOT(cnet, install("loglik"));
	m_loglik = NUMERIC_POINTER(pint)[0];

	if (length(rproblist) != m_numNodes) {
		UNPROTECT(5);
		warning("length(rproblist) != m_numNodes");
		return;
	}

	m_nodeNames     = (char**)CATNET_MALLOC(m_numNodes * sizeof(char*));
	m_numParents    =   (int*)CATNET_MALLOC(m_numNodes * sizeof(int));
	m_parents       =  (int**)CATNET_MALLOC(m_numNodes * sizeof(int*));	
	m_numCategories =   (int*)CATNET_MALLOC(m_numNodes * sizeof(int));

	m_pProbLists = (PROB_LIST<double>**) CATNET_MALLOC(m_numNodes * sizeof(PROB_LIST<double>*));

	if (!m_nodeNames || !m_numParents || !m_parents || 
	    !m_numCategories || !m_pProbLists) {
		if (m_nodeNames)
			CATNET_FREE(m_nodeNames);
		m_nodeNames = 0;
		if (m_numParents)
			CATNET_FREE(m_numParents);
		m_numParents = 0;
		if (m_parents)
			CATNET_FREE(m_parents);
		m_parents = 0;
		if (m_numCategories)
			CATNET_FREE(m_numCategories);
		m_numCategories = 0;
		if (m_pProbLists)
			CATNET_FREE(m_pProbLists);
		m_pProbLists = 0;
		UNPROTECT(5);
		return;
	}

	memset(m_nodeNames,     0, m_numNodes * sizeof(char*));
	memset(m_numParents,    0, m_numNodes * sizeof(int));
	memset(m_parents,       0, m_numNodes * sizeof(int*));
	memset(m_numCategories, 0, m_numNodes * sizeof(int));

	for (i = 0; i < m_numNodes; i++)
		m_pProbLists[i] = 0;

	for(nnode = 0; nnode < m_numNodes; nnode++) {
		pf = VECTOR_ELT(rnodes, nnode);
		m_nodeNames[nnode] = 0;
		if(IS_VECTOR(pf)) {
			pstr = CHAR(asChar(pf));
			if (pstr) {
				m_nodeNames[nnode] = (char*) CATNET_MALLOC((strlen(pstr)+1) * sizeof(char));
				if (m_nodeNames[nnode])
					strcpy(m_nodeNames[nnode], pstr);
			}
		}

		pf = VECTOR_ELT(rparents, nnode);
		m_numParents[nnode] = 0;
		m_parents[nnode]    = 0;
		if (IS_VECTOR(pf)) {
			m_numParents[nnode] = length(pf);
			pn = INTEGER_POINTER(pf);
			m_parents[nnode] = (int*) CATNET_MALLOC(m_numParents[nnode] * sizeof(int));
			if (m_parents[nnode]) {
				for(i = 0; i < m_numParents[nnode]; i++) {
					m_parents[nnode][i] = pn[i] - 1;
				}
			}
		}

		pf = VECTOR_ELT(rcatlist, nnode);
		m_numCategories[nnode] = length(pf);
	}

	// get probabilities
	if(!strcmp(CHAR(asChar(rname)), "catNetworkC")) {
		for (nnode = 0; nnode < m_numNodes; nnode++) {
			rnodeprob = VECTOR_ELT(rproblist, nnode);
			setCondProb(nnode, NUMERIC_POINTER(rnodeprob), length(rnodeprob));   
		}
	}
	else {
		for (nnode = 0; nnode < m_numNodes; nnode++) {
			nodepars = VECTOR_ELT(rparents, nnode);
			nodeproblist = VECTOR_ELT(rproblist, nnode);
			pvec = 0;
			nvec = 0;
			gen_prob_vector(nnode, nodepars, 0, rcatlist, nodeproblist, pvec, nvec);
			setCondProb(nnode, pvec, nvec);
			CATNET_FREE(pvec);
		}
	}

	UNPROTECT(5);
}

SEXP RCATNET_CLASS::genRcatnet(const char * objectName = (const char*)"catNetwork") {

	char str[256];
	int node, i, *pslotcats, *pn;
	double *pf, floglik;
	SEXP plist, ppars, pcats, pnodeprobs, rNodeNames, pint;

	if(strcmp(objectName, "catNetwork") && strcmp(objectName, "catNetworkC") ) 
		return R_NilValue;

	SEXP cnet = PROTECT(NEW_OBJECT(MAKE_CLASS(objectName)));

	PROTECT(plist = allocVector(STRSXP, 1));
	SET_STRING_ELT(plist, 0, mkChar(objectName));
	SET_SLOT(cnet, install("objectName"), plist);
	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_numNodes;
	SET_SLOT(cnet, install("numnodes"), pint);
	UNPROTECT(1);

	PROTECT(rNodeNames = allocVector(STRSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		if(m_nodeNames && m_nodeNames[node]) {
			SET_STRING_ELT(rNodeNames, node, mkChar(m_nodeNames[node]));
		}
		else {
			sprintf(str, "N%d", node+1);
			SET_STRING_ELT(rNodeNames, node, mkChar(str));
		}
	}
	SET_SLOT(cnet, install("nodes"), rNodeNames);
 	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_maxParents;
	SET_SLOT(cnet, install("maxpars"), pint);
	UNPROTECT(1);

	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		if(m_numParents[node] <= 0) {
			SET_VECTOR_ELT(plist, node, R_NilValue);
			continue;
		}
		PROTECT(ppars = NEW_INTEGER(m_numParents[node]));
		pn = INTEGER_POINTER(ppars);
		for(i = 0; i < m_numParents[node]; i++) {
			// remember to increase the index by 1
			pn[i] = m_parents[node][i] + 1;
		}
		SET_VECTOR_ELT(plist, node, ppars);
		UNPROTECT(1);
	}
	SET_SLOT(cnet, install("pars"), plist);
	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_maxCategories;
	SET_SLOT(cnet, install("maxcats"), pint);
	UNPROTECT(1);

	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		if(m_numCategories[node] > m_maxCategories)
			break;
		PROTECT(pcats = allocVector(STRSXP, m_numCategories[node]));
		for(i = 0; i < m_numCategories[node]; i++) {
			if(m_catIndices && m_catIndices[node])
				sprintf(str, "%d", m_catIndices[node][i]);
			else
				sprintf(str, "C%d", i+1);
			SET_STRING_ELT(pcats, i, mkChar(str));
		}
		SET_VECTOR_ELT(plist, node, pcats);
		UNPROTECT(1);
	}
	SET_SLOT(cnet, install("cats"), plist);
	UNPROTECT(1);

	pslotcats = (int*)CATNET_MALLOC(m_maxParents*sizeof(int));
	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	if(!strcmp(objectName, "catNetworkC")) {	
		for(node = 0; node < m_numNodes; node++) {
			if(!m_pProbLists[node])
				continue;
			PROTECT(pnodeprobs = NEW_NUMERIC(m_pProbLists[node]->nProbSize));
			pf = NUMERIC_POINTER(pnodeprobs);
			if (pf && m_pProbLists[node]->pProbs)
				memcpy(pf, m_pProbLists[node]->pProbs, 
					m_pProbLists[node]->nProbSize * sizeof(double));
			SET_VECTOR_ELT(plist, node, pnodeprobs);
			UNPROTECT(1);
		}
	}
	else { 
		for(node = 0; node < m_numNodes; node++) {
			if (pslotcats) {
				memset(pslotcats, 0, m_maxParents*sizeof(int));
			}
			pnodeprobs = genProbList(node, 0, pslotcats);
			SET_VECTOR_ELT(plist, node, pnodeprobs);
			if(pnodeprobs != R_NilValue)
				UNPROTECT(1);
		}
	}
	SET_SLOT(cnet, install("probs"), plist);
	UNPROTECT(1); // plist

	if (pslotcats)
		CATNET_FREE(pslotcats);

	SET_SLOT(cnet, install("meta"), mkString("catNetwork object"));

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = complexity();
	SET_SLOT(cnet, install("complx"), pint);
	UNPROTECT(1); // pint
	
	PROTECT(pint = NEW_NUMERIC(1));
	floglik = loglik();
	if(floglik > -FLT_MAX)
		NUMERIC_POINTER(pint)[0] = floglik;
	else
		NUMERIC_POINTER(pint)[0] = R_NegInf;
	SET_SLOT(cnet, install("loglik"), pint);
	UNPROTECT(1); // pint

	PROTECT(pint = NEW_INTEGER(m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		if(m_pProbLists[node])
			INTEGER_POINTER(pint)[node] = m_pProbLists[node]->sampleSize;
	}
	SET_SLOT(cnet, install("nodeSampleSizes"), pint);
	UNPROTECT(1); // pint

	UNPROTECT(1); // cnet

	return cnet;
}

SEXP RCATNET_CLASS::genRcatnetFromDagEvaluate(SEXP rDagEval, SEXP rDagIndex) {

	char str[256];
	int numDags, nIndex, node, i, n, k, nbits;
	int nIntBuff, *pIntBuff, *pParBuff, *pParIndex;
	char *pByteBuff;
	const char *pstr;
	SEXP rslot, rnodes, plist, poldlist, ppars, pcats, poldcats, rNodeNames, pint;

	if(!isS4(rDagEval))
		CATNET_PARAM_ERR();
	rslot = GET_SLOT(rDagEval, install("numDags"));
	if(rslot == R_NilValue)
		CATNET_PARAM_ERR();

	PROTECT(rslot = AS_INTEGER(rslot));
	numDags = INTEGER_POINTER(rslot)[0];
	UNPROTECT(1);

	PROTECT(rDagIndex = AS_INTEGER(rDagIndex));
	nIndex = INTEGER_POINTER(rDagIndex)[0] - 1;
	UNPROTECT(1);
	if(numDags < 1 || nIndex < 0 || nIndex >= numDags)
		CATNET_PARAM_ERR();

	pint = GET_SLOT(rDagEval, install("numnodes"));
	m_numNodes = INTEGER_POINTER(pint)[0];
	pint = GET_SLOT(rDagEval, install("maxpars"));
	m_maxParents = INTEGER_POINTER(pint)[0];
	pint = GET_SLOT(rDagEval, install("maxcats"));
	m_maxCategories = INTEGER_POINTER(pint)[0];

	if(m_numNodes < 0 || m_maxParents < 0 || m_maxCategories < 0 || 
	   m_numNodes > 1e16 || m_maxParents > 1e4 || m_maxCategories > 1e4) {
		CATNET_PARAM_ERR();
	}

	/* initialize the CATNET attributes */
	init(m_numNodes, m_maxParents, m_maxCategories);

	rnodes = GET_SLOT(rDagEval, install("nodes"));
	PROTECT(rnodes = AS_LIST(rnodes));
	for(node = 0; node < m_numNodes; node++) {
		rNodeNames = VECTOR_ELT(rnodes, node);
		if(m_nodeNames[node])
			CATNET_FREE(m_nodeNames[node]);
		m_nodeNames[node] = 0;
		if(length(rNodeNames) > 0) {
			pstr = CHAR(asChar(rNodeNames));
			m_nodeNames[node] = (char*)CATNET_MALLOC((strlen(pstr)+1) * sizeof(char));
			if (m_nodeNames[node] && pstr)
				strcpy(m_nodeNames[node], pstr);
		}
	}
	UNPROTECT(1);

	pint = GET_SLOT(rDagEval, install("complx"));
	if(length(pint) != numDags)
		CATNET_ERR("Wrong complx slot");
	m_complexity = INTEGER_POINTER(pint)[nIndex];
	pint = GET_SLOT(rDagEval, install("loglik"));
	if(length(pint) != numDags)
		CATNET_ERR("Wrong loglik slot");
	m_loglik = NUMERIC_POINTER(pint)[nIndex];

	plist = GET_SLOT(rDagEval, install("numPars"));
	PROTECT(plist = AS_LIST(plist));
	if(length(plist) != numDags)
		CATNET_ERR("Wrong numPars slot");
	ppars = VECTOR_ELT(plist, nIndex);
	PROTECT(ppars = AS_INTEGER(ppars));
	nIntBuff = length(ppars);

	pIntBuff  = (int*)CATNET_MALLOC(nIntBuff*sizeof(int));
	pParIndex = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));

	if (!pIntBuff || !pParIndex) {
		CATNET_PARAM_ERR();
	}

	memcpy(pIntBuff, INTEGER_POINTER(ppars), nIntBuff*sizeof(int));
	memset(pParIndex, 0, m_numNodes*sizeof(int));

	UNPROTECT(2); // ppars & plist

	nbits = *((int*)pIntBuff);
	if(nbits <= 8) {
		pByteBuff = (char*)pIntBuff;
		i = sizeof(int);
		node = 0;
		while(i < nIntBuff*(int)sizeof(int) && node < m_numNodes) {
			n = (int)pByteBuff[i++];
			if(n < 0) {
				k = -n;
				n = (int)pByteBuff[i++];
				while(k-- > 0)
					pParIndex[node++] = n;
			}
			else 
				pParIndex[node++] = n;
		}
		if(node < m_numNodes) {
			sprintf(str, "Wrong numPars[%d] slot", nIndex+1);
			CATNET_ERR(str);
		}
	}
	else {
		i = 1;
		node = 0;
		while(i < nIntBuff*(int)sizeof(int) && node < m_numNodes) {
			n = (int)pIntBuff[i++];
			if(n < 0) {
				k = -n;
				n = (int)pIntBuff[i++];
				while(k-- > 0)
					pParIndex[node++] = n;
			}
			else 
				pParIndex[node++] = n;
		}
		if(node < m_numNodes) {
			sprintf(str, "Wrong numPars[%d] slot", nIndex+1);
			CATNET_ERR(str);
		}
	}

	pParBuff = (int*)CATNET_MALLOC(m_maxParents*sizeof(int));
	if (!pParBuff) {
		CATNET_FREE(pIntBuff);
		CATNET_FREE(pParIndex);
		CATNET_PARAM_ERR();
	}

 	plist = GET_SLOT(rDagEval, install("parSlots"));
	if(length(plist) != m_numNodes)
		CATNET_ERR("Wrong parSlots slot");
	PROTECT(plist = AS_LIST(plist));
	for(node = 0; node < m_numNodes; node++) {
		PROTECT(ppars = AS_INTEGER(ppars));
		ppars = VECTOR_ELT(plist, node);
		if(pParIndex[node] < 0 || length(ppars) < (pParIndex[node]+1)*m_maxParents) {
			sprintf(str, "Wrong parSlots for node %d", node+1);
			CATNET_ERR(str);
		}
		if (m_maxParents > 0)
			memcpy(pParBuff, INTEGER_POINTER(ppars) + pParIndex[node]*m_maxParents, m_maxParents*sizeof(int));
		UNPROTECT(1);
		i = 0;
		while(pParBuff[i] != 0 && i < m_maxParents) 
			i++;
		m_numParents[node] = i;

		if(m_numParents[node] <= 0) 
			continue;
		m_parents[node] = (int*)CATNET_MALLOC(m_numParents[node] * sizeof(int));
		if (m_parents[node] && m_numParents[node] > 0)
			memcpy(m_parents[node], pParBuff, m_numParents[node]*sizeof(int));
	}
	UNPROTECT(1); // plist

	CATNET_FREE(pParBuff);
	pParBuff = 0;
	CATNET_FREE(pIntBuff);
	pIntBuff = 0;

	SEXP cnet = PROTECT(NEW_OBJECT(MAKE_CLASS("catNetwork")));

	PROTECT(plist = allocVector(STRSXP, 1));
	SET_STRING_ELT(plist, 0, mkChar("catNetwork"));
	SET_SLOT(cnet, install("objectName"), plist);
	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_numNodes;
	SET_SLOT(cnet, install("numnodes"), pint);
	UNPROTECT(1);

	PROTECT(rNodeNames = allocVector(STRSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		if(m_nodeNames && m_nodeNames[node]) {
			SET_STRING_ELT(rNodeNames, node, mkChar(m_nodeNames[node]));
		}
		else {
			sprintf(str, "N%d", node+1);
			SET_STRING_ELT(rNodeNames, node, mkChar(str));
		}
	}
	SET_SLOT(cnet, install("nodes"), rNodeNames);
 	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_maxParents;
	SET_SLOT(cnet, install("maxpars"), pint);
	UNPROTECT(1);

	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		if(m_numParents[node] <= 0) {
			SET_VECTOR_ELT(plist, node, R_NilValue);
			continue;
		}
		PROTECT(ppars = NEW_INTEGER(m_numParents[node]));
		if (m_parents[node] && m_numParents[node] > 0)
			memcpy(INTEGER_POINTER(ppars), m_parents[node], m_numParents[node]*sizeof(int));
		SET_VECTOR_ELT(plist, node, ppars);
		UNPROTECT(1);
	}
	setAttrib(plist, R_NamesSymbol, rNodeNames);
	SET_SLOT(cnet, install("pars"), plist);
	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_maxCategories;
	SET_SLOT(cnet, install("maxcats"), pint);
	UNPROTECT(1);

	poldlist = GET_SLOT(rDagEval, install("cats"));
	if(length(poldlist) != m_numNodes)
		CATNET_ERR("Wrong categories slot");
	PROTECT(poldlist = AS_LIST(poldlist));
	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		poldcats = AS_LIST(VECTOR_ELT(poldlist, node));
		m_numCategories[node] = length(poldcats);
		if(m_numCategories[node] > m_maxCategories)
			CATNET_ERR("Wrong categories slot");
		PROTECT(pcats = allocVector(STRSXP, m_numCategories[node]));
		for(i = 0; i < m_numCategories[node]; i++) {
			SET_STRING_ELT(pcats, i, asChar(VECTOR_ELT(poldcats, i)));
		}
		SET_VECTOR_ELT(plist, node, pcats);
		UNPROTECT(1);
	}
	setAttrib(plist, R_NamesSymbol, rNodeNames);
	SET_SLOT(cnet, install("cats"), plist);
	UNPROTECT(2);

	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	SET_SLOT(cnet, install("probs"), plist);
	UNPROTECT(1);

	SET_SLOT(cnet, install("meta"), mkString("dagEvaluate object"));

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_complexity;
	SET_SLOT(cnet, install("complx"), pint);
	UNPROTECT(1);
	
	PROTECT(pint = NEW_NUMERIC(1));
	NUMERIC_POINTER(pint)[0] = m_loglik;
	SET_SLOT(cnet, install("loglik"), pint);
	UNPROTECT(1);

	plist = GET_SLOT(rDagEval, install("parSampleSize"));
	if(length(plist) != m_numNodes)
		CATNET_ERR("Wrong parSampleSize slot");
	PROTECT(plist = AS_LIST(plist));
	PROTECT(pint = NEW_INTEGER(m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		ppars = VECTOR_ELT(plist, node);
		INTEGER_POINTER(pint)[node] = INTEGER_POINTER(ppars)[pParIndex[node]];
	}
	SET_SLOT(cnet, install("nodeSampleSizes"), pint);
	UNPROTECT(2);

	CATNET_FREE(pParIndex);

	UNPROTECT(1); // cnet
	return cnet;
}

SEXP RCATNET_CLASS::genProbList(int node, int paridx, int *pcats) {
	int j, npar;
	SEXP problist;
	double *pslot, *pp;

	if(m_pProbLists == 0 || m_pProbLists[node] == 0 || paridx < 0)
		return R_NilValue;

	if(paridx >= m_numParents[node]) {
		pslot = m_pProbLists[node]->find_slot(0, pcats, 0);
		PROTECT(problist = NEW_NUMERIC(m_numCategories[node]));
		pp = NUMERIC_POINTER(problist);
		if (pp && pslot && m_numCategories[node] > 0)
			memcpy(pp, pslot, m_numCategories[node]*sizeof(double));
		return problist;
	}

	npar = m_parents[node][paridx];
	PROTECT(problist = allocVector(VECSXP, m_numCategories[npar]));
	for(j = 0; j < m_numCategories[npar]; j++) {
		pcats[paridx] = j;
		SET_VECTOR_ELT(problist, j, genProbList(node, paridx + 1, pcats));
		UNPROTECT(1);
	}

	return problist;
}

SEXP RCATNET_CLASS::genSamples(SEXP rNumSamples, SEXP rPerturbations, SEXP rNaRate) {

	SEXP rsamples = R_NilValue;
	int numsamples;
	int *pSamples, *pRsamples, *pPerturbations;
	double fNaRate;
	int i, j, k, nnode, *pnodepars, *pnodesample, *porder;
	double u, v, *pnodeprob;
	PROB_LIST<double>* pProbList; 
	
	PROTECT(rNumSamples = AS_INTEGER(rNumSamples));
	numsamples = INTEGER_POINTER(rNumSamples)[0];
	UNPROTECT(1);

	PROTECT(rNaRate = AS_NUMERIC(rNaRate));
	fNaRate = NUMERIC_POINTER(rNaRate)[0];
	UNPROTECT(1);

	pPerturbations = 0;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
		pPerturbations = INTEGER(rPerturbations);
	}

	porder = getOrder();
	if(!porder)
		return R_NilValue;

	pSamples = (int*)CATNET_MALLOC(m_numNodes*numsamples*sizeof(int));
	if(!pSamples) {
		CATNET_FREE(porder);
		return R_NilValue;
	}

	pnodesample = 0;
	if(m_maxParents > 0) {
		pnodesample = (int*)CATNET_MALLOC(m_maxParents * sizeof(int));
		if(!pnodesample) {
			CATNET_FREE(pSamples);
			CATNET_FREE(porder);
			return R_NilValue;
		}
	}

	GetRNGstate();
	for(k = 0; k < m_numNodes; k++) {
		nnode = porder[k];
		pnodepars = m_parents[nnode];
		pProbList = (PROB_LIST<double>*)getNodeProb(nnode);
		if (!pProbList) {
			if (pnodesample)
				CATNET_FREE(pnodesample);
			CATNET_FREE(pSamples);
			CATNET_FREE(porder);
			return R_NilValue;
		}

		for (j = 0; j < numsamples; j++) {

			if(pPerturbations) {
				if(!R_IsNA(pPerturbations[j*m_numNodes + nnode]) && 
					pPerturbations[j*m_numNodes + nnode] >= 1 && pPerturbations[j*m_numNodes + nnode] <= m_numCategories[nnode]) {
					pSamples[j * m_numNodes + nnode] = pPerturbations[j * m_numNodes + nnode];
					continue;
				}
			}

			for (i = 0; i < m_numParents[nnode]; i++) {
				if (pnodepars[i] < 0 || pnodepars[i] >= m_numNodes)
					break;
				pnodesample[i] = (int)(pSamples[j * m_numNodes + pnodepars[i]] - 1);
			}
			pnodeprob = pProbList->find_slot(0, pnodesample, 0);

			u = (double)unif_rand();
			v = 0;
			for(i = 0; i < m_numCategories[nnode]; i++) {
				v += pnodeprob[i];
				if(u <= v)
					break;
			}
			pSamples[j * m_numNodes + nnode] = i + 1;
		}
	}

	k = (int)(fNaRate*m_numNodes);
	if(k > 0 && k < m_numNodes) {
		int ii, fmax, *paux = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
		if (!paux) {
			if (pnodesample)
				CATNET_FREE(pnodesample);
			CATNET_FREE(pSamples);
			CATNET_FREE(porder);
			return R_NilValue;
		}
		for (j = 0; j < numsamples; j++) {
			for(ii = 0; ii < m_numNodes; ii++)
				paux[ii] = (int)(RAND_MAX*unif_rand());
			for(i = 0; i < k; i++) {
				nnode = 0;
				fmax = -(int)RAND_MAX;
				for(ii = 0; ii < m_numNodes; ii++) {
					if(fmax < paux[ii]) {
						fmax = paux[ii];
						nnode = ii;
					}
				}
				if(nnode >= 0 && nnode < m_numNodes) {
					paux[nnode] = -(int)RAND_MAX;
					pSamples[j * m_numNodes + nnode] = R_NaInt;
				}
			}
		}
		CATNET_FREE(paux);
	}
	PutRNGstate();

	if(!isNull(rPerturbations))
		UNPROTECT(1);

	if(pnodesample)
		CATNET_FREE(pnodesample);		
	if(porder)
		CATNET_FREE(porder);

	// output the new matrix
	PROTECT(rsamples = NEW_INTEGER(m_numNodes * numsamples));
	pRsamples = (int*)INTEGER_POINTER(rsamples);
	if (pRsamples && pSamples)
		memcpy(pRsamples, pSamples, m_numNodes*numsamples*sizeof(int));
	UNPROTECT(1);

	CATNET_FREE(pSamples);

	return rsamples;
}

#endif

