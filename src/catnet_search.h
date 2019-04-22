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
 * catnet_search.h
 * includes soft edge constraints in addition to the original CATNET_SEARCH
 *
 *  Created on: Sep 25, 2009
 *      Author: nbalov
 */

#include "utils.h"
#include "catnetd.h"
#include "thread.h"
#include "cache.h"
#include "search_params.h"

#if (defined(DISCRETE_SAMPLE) && !defined(CATNETD_SEARCH_H)) || (defined(PROB_SAMPLE) && !defined(CATNETP_SEARCH_H))

#ifdef DISCRETE_SAMPLE
#define CATNETD_SEARCH_H
#else
#ifdef PROB_SAMPLE
#define CATNETP_SEARCH_H
#else
#error "No proper CATNET class"
#endif
#endif

#ifdef DISCRETE_SAMPLE
template<class t_prob>
class CATNETD_SEARCH:  public c_thread, public c_cache {
protected:
	int m_nCatnets;
	CATNETD<t_prob> **m_pCatnets;
	int m_numNodes, *m_pNodeNumCats, **m_pNodeCats, m_numSamples;
public:
	CATNETD_SEARCH() {
		m_nCatnets     = 0;
		m_pCatnets     = 0;
		m_numNodes     = 0;
		m_pNodeNumCats = 0;
		m_pNodeCats    = 0;
	}

	~CATNETD_SEARCH() {
		_release();
	}

	CATNETD<t_prob> **catnets() {
		return m_pCatnets;
	}

#else
template<class t_prob>
class CATNETP_SEARCH:  public c_thread, public c_cache {

protected:
	int m_nCatnets;
	CATNETP<t_prob> **m_pCatnets;
	int m_numNodes, *m_pNodeNumCats, **m_pNodeCats, m_numSamples;
public:
	CATNETP_SEARCH() {
		m_nCatnets     = 0;
		m_pCatnets     = 0;
		m_numNodes     = 0;
		m_pNodeNumCats = 0;
		m_pNodeCats    = 0;
	}

	~CATNETP_SEARCH() {
		_release();
	}

	CATNETP<t_prob> **catnets() {
		return m_pCatnets;
	}
#endif

	int numCatnets() {
		return m_nCatnets;
	}

protected:
	void _release() {
		int i;
		if(m_pCatnets) {
			for(i = 0; i < m_nCatnets; i++)
				if(m_pCatnets[i]) {
					delete m_pCatnets[i];
					m_pCatnets[i] = 0;
				}
			CATNET_FREE(m_pCatnets);
		}
		m_pCatnets = 0;
		m_nCatnets = 0;
		if(m_pNodeCats) {
			for(i = 0; i < m_numNodes; i++) 
				if(m_pNodeCats[i])
					CATNET_FREE(m_pNodeCats[i]);
			CATNET_FREE(m_pNodeCats);
		}
		if(m_pNodeNumCats) 
			CATNET_FREE(m_pNodeNumCats);
	}

public:
	/* pSamples and perturbations are sample=columns and node=rows. */
	/* Each parentsPool[i] is numNodes long ! */
	int search(SEARCH_PARAMETERS *pestim) {

		if(!pestim)
			return 0;
		int numNodes = pestim->m_numNodes;
		int numSamples = pestim->m_numSamples;
		int *perturbations = pestim->m_pPerturbations;
		int maxParentSet = pestim->m_maxParentSet;
		int *parSizes = pestim->m_pParentSizes;
		int maxComplexity = pestim->m_maxComplexity;
		int **parentsPool = pestim->m_parentsPool;
		int **fixedParentsPool = pestim->m_fixedParentsPool;
		double *matEdgeLiks = pestim->m_matEdgeLiks;
		int becho = pestim->m_echo;

		int i, j, k, d, ncomb, ncombMaxLogLik, nnode, nodecomplx;
		int nocache;
		int maxCategories, numSubSamples, complx, bEqualCategories;
		int *parset, parsetsize, *idparset, *fixparset, fixparsetsize, *bernoulibuff;
		int *paux, **pcomblist, ncomblist, maxpars, ballow, bfixallow;

		t_prob fLogLik, fMaxLogLik, tempLogLik, priorLogLik;
		PROB_LIST<t_prob> probMaxNode, *pProbNode;

#ifdef DISCRETE_SAMPLE
		int *pSamples = (int*)pestim->m_pSamples;
		int *pSubSamples, mincat, maxcat;
		CATNETD<t_prob> baseCatnet, *pNewNet, **pCurCatnetList;
#else		
		t_prob *pSamples = (t_prob*)pestim->m_pSamples;
		t_prob *pSubSamples;
		CATNETP<t_prob> baseCatnet, *pNewNet, **pCurCatnetList;
		int nline;
#endif

		int *pSubClasses = 0;
		int *pClasses = pestim->m_pClasses;

		_release();

		if(numNodes < 1 || numSamples < 1 || !pSamples)
			return CATNET_ERR_PARAM;
		if(maxComplexity < numNodes)
			maxComplexity = numNodes;
		
		m_numNodes    = numNodes;
		m_numSamples  = numSamples;
		maxCategories = 0;

		m_pNodeNumCats = (int*)CATNET_MALLOC(numNodes*sizeof(int));
		if (!m_pNodeNumCats) 
			return CATNET_ERR_MEM;
		m_pNodeCats    = (int**)CATNET_MALLOC(numNodes*sizeof(int*));
		if (!m_pNodeCats) { 
			CATNET_FREE(m_pNodeNumCats);
			return CATNET_ERR_MEM;
		}

		memset(m_pNodeCats,    0, numNodes*sizeof(int*));
		memset(m_pNodeNumCats, 0, numNodes*sizeof(int));

		if(pestim->m_pNodeNumCats && pestim->m_pNodeCats) {
			memcpy(m_pNodeNumCats, pestim->m_pNodeNumCats, numNodes*sizeof(int));
			for(i = 0; i < numNodes; i++) {
				m_pNodeCats[i] = (int*)CATNET_MALLOC(m_pNodeNumCats[i]*sizeof(int));
				if (!m_pNodeCats[i])
					return CATNET_ERR_MEM;
				if (pestim->m_pNodeCats && pestim->m_pNodeCats[i])
					memcpy(m_pNodeCats[i], pestim->m_pNodeCats[i], m_pNodeNumCats[i]*sizeof(int));
			}
		}

#ifdef DISCRETE_SAMPLE
		else { 
			for(i = 0; i < numNodes; i++) {
				mincat = INT_MAX;
				maxcat = -INT_MAX;
				for(j = 0; j < numSamples; j++) {
					if(pSamples[j*numNodes + i] == CATNET_NAN)
						continue;
					if(pSamples[j*numNodes + i] < mincat)
						mincat = pSamples[j*numNodes + i];
					if(pSamples[j*numNodes + i] > maxcat)
						maxcat = pSamples[j*numNodes + i];
				}
				m_pNodeNumCats[i] = maxcat - mincat + 1;
				m_pNodeCats[i] = (int*)CATNET_MALLOC(m_pNodeNumCats[i]*sizeof(int));
				if (!m_pNodeCats[i]) {
					return CATNET_ERR_MEM;
				}
				for(j = 0; j < m_pNodeNumCats[i]; j++)
					m_pNodeCats[i][j] = mincat + j;
				/* order m_pNodeNumCats[i] */
				for(j = 0; j < m_pNodeNumCats[i]; j++) {
					for(k = j + 1; k < m_pNodeNumCats[i]; k++) {
						if(m_pNodeCats[i][j] > m_pNodeCats[i][k]) {
							d = m_pNodeCats[i][j]; 
							m_pNodeCats[i][j] = m_pNodeCats[i][k];
							m_pNodeCats[i][k] = d;
						}
					}
				}
			}
		}
		for(i = 0; i < numNodes; i++) { 
			for(j = 0; j < numSamples; j++) {
				if(pSamples[j*numNodes + i] == CATNET_NAN) {
					// CATNET_NANs will be assigned m_pNodeNumCats[i]
					pSamples[j*numNodes + i] = m_pNodeNumCats[i];
					continue;
				}
				for(d = 0; d < m_pNodeNumCats[i]; d++)
					if(m_pNodeCats[i][d] == pSamples[j*numNodes + i])
						break;
				if(d >= m_pNodeNumCats[i]) {
					CATNET_ERR("Incorrect categories");
					return CATNET_ERR_PARAM;
				}
				pSamples[j*numNodes + i] = d;
			}
			if(maxCategories < m_pNodeNumCats[i])
				maxCategories = m_pNodeNumCats[i];
			if(i > 1 && m_pNodeNumCats[i] != m_pNodeNumCats[0])
				bEqualCategories = 0;
		}
#else
		nline = 0;
		maxCategories = 0;
		for(i = 0; i < numNodes; i++) { 
			nline += m_pNodeNumCats[i];
			if(maxCategories < m_pNodeNumCats[i])
				maxCategories = m_pNodeNumCats[i];
			if(i > 1 && m_pNodeNumCats[i] != m_pNodeNumCats[0])
				bEqualCategories = 0;
		}
#endif

		bEqualCategories = 1;
		for(i = 0; i < numNodes; i++) {
			if(i > 1 && m_pNodeNumCats[i] != m_pNodeNumCats[0])
				bEqualCategories = 0;
		}

		parset    = (int*)CATNET_MALLOC(numNodes*sizeof(int));
		idparset  = (int*)CATNET_MALLOC(numNodes*sizeof(int));
		fixparset = (int*)CATNET_MALLOC(numNodes*sizeof(int));
		if (!parset || !idparset || !fixparset) {
			return CATNET_ERR_MEM;
		}
		bernoulibuff = 0;
		if(matEdgeLiks)
			bernoulibuff = (int*)CATNET_MALLOC(numNodes*sizeof(int));

		/* parent pools */
		if(pestim->m_maxParentsPool >= 1 && pestim->m_matNodeCondLiks != 0 && parentsPool != 0) {
			double fs,ff,fsum;
			for(i = 0; i < numNodes; i++) {
				for(j = 0; j < pestim->m_maxParentsPool; j++) 
					parentsPool[i][j] = -1;
				if(i < pestim->m_maxParentsPool) {
					for(j = 0; j < i; j++) 
						parentsPool[i][j] = j;
					continue;
				}
				/* sample from {m_matNodeCondLiks[i*numNodes+.]} */
				fsum = 0;
				for(j = 0; j < i; j++) {
					fsum += pestim->m_matNodeCondLiks[i*numNodes+j];
				}
				k = 0;
				ncomb = 0;
				GetRNGstate();
				while(k < pestim->m_maxParentsPool && ncomb < numNodes*numNodes) {
					ncomb++;
					ff = fsum * (double)unif_rand();
					fs = 0;
					j = 0;
					while(fs <= ff && j < i) {
						fs += pestim->m_matNodeCondLiks[i*numNodes+j];
						j++;
					}
					j--;
					for(d = 0; d < k; d++) {
						if(parentsPool[i][d] == j)
							break;
					}
					if(d != k || j >= numNodes)
						continue;
					parentsPool[i][k] = j;
					k++;
				}
				PutRNGstate();
			}
		}

#ifdef DISCRETE_SAMPLE
		m_nCatnets = maxComplexity + 1;
		m_pCatnets = (CATNETD<t_prob>**)CATNET_MALLOC(m_nCatnets*sizeof(CATNETD<t_prob>));
		if (!m_pCatnets) {
			return CATNET_ERR_MEM;
		}
		memset(m_pCatnets, 0, m_nCatnets*sizeof(CATNETD<t_prob>*));

		pCurCatnetList = (CATNETD<t_prob>**)CATNET_MALLOC(m_nCatnets*sizeof(CATNETD<t_prob>*));
		
pSubSamples = 0;
		if(perturbations)
			pSubSamples = (int*)CATNET_MALLOC(numNodes*numSamples*sizeof(int));

		/* create a network without edges*/
		pNewNet = new CATNETD<t_prob>(numNodes, 0/*maxParentSet*/, maxCategories, 0, 0, 0, m_pNodeNumCats);
#else
		m_nCatnets = maxComplexity + 1;
		m_pCatnets = (CATNETP<t_prob>**)CATNET_MALLOC(m_nCatnets*sizeof(CATNETP<t_prob>));
		if (!m_pCatnets) {
			return CATNET_ERR_MEM;
		}
		memset(m_pCatnets, 0, m_nCatnets*sizeof(CATNETP<t_prob>*));

		pCurCatnetList = (CATNETP<t_prob>**)CATNET_MALLOC(m_nCatnets*sizeof(CATNETP<t_prob>*));

		pSubSamples = 0;
		if(perturbations) 
			pSubSamples = (t_prob*)CATNET_MALLOC(nline*numSamples*sizeof(t_prob));
		/* create a network without edges*/
		pNewNet = new CATNETP<t_prob>(numNodes, 0/*maxParentSet*/, maxCategories, 0, 0, 0, m_pNodeNumCats);
#endif

		if (!pCurCatnetList || !pNewNet) {
			if (m_pCatnets)
				CATNET_FREE(m_pCatnets);
			if (pCurCatnetList)
				CATNET_FREE(pCurCatnetList);
			if (pNewNet)
				delete pNewNet;
			CATNET_FREE(parset);
			CATNET_FREE(fixparset);
			CATNET_FREE(idparset);
			return CATNET_ERR_MEM;
		}

		pSubClasses = 0;
		if(pClasses)
			pSubClasses = (int*)CATNET_MALLOC(numSamples*sizeof(int));

		/* set parents */
		for(nnode = 0; nnode < numNodes; nnode++) {
			fixparsetsize = 0;
			if(fixedParentsPool && fixedParentsPool[nnode]) {	
				for(j = 0; j < nnode; j++) {
					ballow = 1;
					if(parentsPool && parentsPool[nnode]) {
						ballow = 0;
						for(k = 0; k < pestim->m_maxParentsPool; k++) {
							if(parentsPool[nnode][k] < 0)
								break;
							if(j == parentsPool[nnode][k])
								ballow = 1;
						}
					}
					if(parentsPool && !parentsPool[nnode])
						ballow = 0;
					bfixallow = 0;
					if(fixedParentsPool && fixedParentsPool[nnode]) {
						for(k = 0; k < pestim->m_maxParentsPool; k++) {
							if(fixedParentsPool[nnode][k] < 0)
								break;
							if(j == fixedParentsPool[nnode][k])
								bfixallow = 1;
						}
					}
					if(!ballow)
						continue;
					if(bfixallow) {
					  fixparset[fixparsetsize] = j;
					  fixparsetsize++;
					}
				}
				if(fixparsetsize > 0) {
					pNewNet -> setParents(nnode, fixparset, fixparsetsize);
				}
			}
			// at this point m_pProbLists are not initialized
			// set sample probabilities and calculate log-likelihood
			numSubSamples = numSamples;
			if(perturbations && pSubSamples) {
				numSubSamples = 0;
				for(j = 0; j < numSamples; j++) {
					if(!perturbations[j * numNodes + nnode]) {
#ifdef DISCRETE_SAMPLE
						memcpy(pSubSamples + numSubSamples*numNodes, 
							pSamples + j*numNodes, numNodes*sizeof(int));
#else
						memcpy(pSubSamples + numSubSamples*nline, 
							pSamples + j*nline, nline*sizeof(t_prob));
#endif
						if(pClasses && pSubClasses) 
							pSubClasses[numSubSamples] = pClasses[j];
						numSubSamples++;
					}
				}
				if(pClasses)
					pNewNet->setNodeSampleProbKL(nnode, pSubSamples, numSubSamples, pSubClasses, /*bNormalize = */0, /*bUsePearson = */pestim->m_klmode==1);
				else
					pNewNet->setNodeSampleProb(nnode, pSubSamples, numSubSamples);
			}
			else {
				if(pClasses)
					pNewNet->setNodeSampleProbKL(nnode, pSamples, numSamples, pClasses, /*bNormalize = */0, /*bUsePearson = */pestim->m_klmode==1);
				else
					pNewNet->setNodeSampleProb(nnode, pSamples, numSamples);
			}

			priorLogLik = 0;
			if(matEdgeLiks && bernoulibuff) {
				tempLogLik = 1;
				memset(bernoulibuff, 0, numNodes*sizeof(int));
				for(i = 0; i < fixparsetsize; i++) {
					tempLogLik *= (t_prob)matEdgeLiks[fixparset[i]*numNodes+nnode];
					bernoulibuff[fixparset[i]] = 1;
				}
				if(tempLogLik > 0)
					priorLogLik = (t_prob)log((double)tempLogLik);
				else 
					priorLogLik = -FLT_MAX;
				tempLogLik = 1;
				for(i = 0; i < numNodes; i++) 
					if(i != nnode && bernoulibuff[i] == 0) {
						tempLogLik *= (1 - (t_prob)matEdgeLiks[i*numNodes+nnode]);
					}
				if(tempLogLik > 0)
					tempLogLik = (t_prob)log((double)tempLogLik);
				else 
					tempLogLik = -FLT_MAX;
				if(priorLogLik == -FLT_MAX || tempLogLik == -FLT_MAX)
					priorLogLik = -FLT_MAX;
				else {
					priorLogLik += tempLogLik;
				}
				pNewNet -> setNodePriorProb(nnode, priorLogLik);
			}
		}

		baseCatnet.init(numNodes, maxParentSet, maxCategories, 0, 0, 0, m_pNodeNumCats);

		complx = pNewNet->complexity();
		m_pCatnets[complx] = pNewNet;

		/* main loop of consequential non-empty-parenthood-node additions */
		for(nnode = 1; nnode < numNodes; nnode++) {

			if(_wait_stop_event(4/*millisecs*/) == 0) {
				Rprintf("STOP signal detected\n");
				break;
			}

			if(becho) {
				Rprintf("processing node %d\n", nnode+1);
				Rprintf("    [#parents][#combinations] = ");
			}

			fixparsetsize = 0;
			parsetsize = 0;
			for(j = 0; j < nnode; j++) {
				ballow = 1;
				if(parentsPool && parentsPool[nnode]) {
					ballow = 0;
					for(k = 0; k < pestim->m_maxParentsPool; k++) {
						if(parentsPool[nnode][k] < 0)
							break;
						if(j == parentsPool[nnode][k])
							ballow = 1;
					}
				}
				if(parentsPool && !parentsPool[nnode])
					ballow = 0;
				bfixallow = 0;
				if(fixedParentsPool && fixedParentsPool[nnode]) {
					for(k = 0; k < pestim->m_maxParentsPool; k++) {
						if(fixedParentsPool[nnode][k] < 0)
							break;
						if(j == fixedParentsPool[nnode][k])
							bfixallow = 1;
					}
				}
				if(!ballow)
					continue;
				if(bfixallow) {
				  fixparset[fixparsetsize] = j;
				  fixparsetsize++;
				}
				else {
				  parset[parsetsize] = j;
				  parsetsize++;
				}
			}
			/* extend the content before sending to cache; parsetsize + fixparsetsize < numNodes */
			memcpy(parset + parsetsize, fixparset, fixparsetsize*sizeof(int));

			/* check out wheather the parent pool has equal number of categories */
			bEqualCategories = 1;
			for(j = 0; j < parsetsize + fixparsetsize; j++) {
				if(j > 0 && m_pNodeNumCats[parset[j]] != m_pNodeNumCats[parset[0]])
					bEqualCategories = 0;
			}

			maxpars = maxParentSet;
			if(parSizes && parSizes[nnode] < maxParentSet)
				maxpars = parSizes[nnode];

			if(maxpars > parsetsize + fixparsetsize)
				maxpars = parsetsize + fixparsetsize;

			memset(pCurCatnetList, 0, m_nCatnets*sizeof(CATNETD<t_prob>*));

			if(bEqualCategories) {

			for(d = fixparsetsize + 1; d <= maxpars; d++) {

				//if(_wait_stop_event(1/*millisecs*/) == 0)
				//	break;

				nocache = 1;
				if(!pestim->m_pCacheMutex) {
					nocache = !getCachedProb(parset, parsetsize + fixparsetsize, nnode, 
						idparset, d, &probMaxNode, &fMaxLogLik);
				}
				else {
					MUTEX_LOCK(pestim->m_pCacheMutex);
					nocache = !getCachedProb(parset, parsetsize + fixparsetsize, nnode, 
						idparset, d, &probMaxNode, &fMaxLogLik);
					MUTEX_UNLOCK(pestim->m_pCacheMutex);
				}

				if(nocache) { 

					pcomblist = 0;
					ncomblist = 0;
					_combination_sets<int>(pcomblist, ncomblist, 0, parset, parsetsize, 0, d - fixparsetsize);

				        if(fixparsetsize > 0) {
				        	if(!pcomblist || ncomblist < 1) {
				        	    	pcomblist = (int**)CATNET_MALLOC(1*sizeof(int*));
							if (!pcomblist)
								return CATNET_ERR_MEM;
				            		pcomblist[0] = 0;	
				            		ncomblist    = 1;
				          	}
				        	for(k = 0; k < ncomblist; k++) {
				            		paux = (int*)CATNET_MALLOC(d*sizeof(int));
							if (!paux)
								return CATNET_ERR_MEM;
							for(j = 0; j < fixparsetsize; j++)
				            			paux[j] = fixparset[j];
					        	if(pcomblist[k] && d > fixparsetsize) {
				            			memcpy(paux + fixparsetsize, pcomblist[k], (d-fixparsetsize)*sizeof(int));
							}
				            		if(pcomblist[k])
				            			CATNET_FREE(pcomblist[k]); 
				           		pcomblist[k] = paux;
						}
					}

					if(becho)
						Rprintf("[%d]%d  ", d, ncomblist);

					fMaxLogLik = -FLT_MAX;
					ncombMaxLogLik = -1;
					probMaxNode.reset();

					for(ncomb = 0; ncomb < ncomblist; ncomb++) {
						// add pcomplist[j] parent set to nnode
						baseCatnet.setParents(nnode, pcomblist[ncomb], d);
			     
						// add perturbation
						numSubSamples = numSamples;
						if(perturbations && pSubSamples) {
							numSubSamples = 0;
							for(j = 0; j < numSamples; j++) {
								if(!perturbations[j * numNodes + nnode]) {
#ifdef DISCRETE_SAMPLE
									memcpy(pSubSamples + numSubSamples*numNodes, 
										pSamples + j*numNodes, numNodes*sizeof(int));
#else
									memcpy(pSubSamples + numSubSamples*nline, 
										pSamples + j*nline, nline*sizeof(t_prob));
#endif									
									if(pClasses) pSubClasses[numSubSamples] = pClasses[j];
									numSubSamples++;
								}
							}
							if(pClasses)
								fLogLik = baseCatnet.setNodeSampleProbKL(nnode, pSubSamples, numSubSamples, pSubClasses, /*bNormalize = */0, /*bUsePearson = */pestim->m_klmode==1);
							else
								fLogLik = baseCatnet.setNodeSampleProb(nnode, pSubSamples, numSubSamples);
						}
						else {
							if(pClasses)
								fLogLik = baseCatnet.setNodeSampleProbKL(nnode, pSamples, numSamples, pClasses, /*bNormalize = */0, /*bUsePearson = */pestim->m_klmode==1);
							else
								fLogLik = baseCatnet.setNodeSampleProb(nnode, pSamples, numSamples);
						}

						/*prior*/
						priorLogLik = 0;
						if(matEdgeLiks && bernoulibuff) {
							tempLogLik = 1;
							memset(bernoulibuff, 0, numNodes*sizeof(int));
							for(i = 0; i < d; i++) {
								tempLogLik *= (t_prob)matEdgeLiks[pcomblist[ncomb][i]*numNodes+nnode];
								bernoulibuff[pcomblist[ncomb][i]] = 1;
							}
							if(tempLogLik > 0)
								priorLogLik = (t_prob)log((double)tempLogLik);
							else 
								priorLogLik = -FLT_MAX;
							tempLogLik = 1;
							for(i = 0; i < numNodes; i++) 
								if(i != nnode && bernoulibuff[i] == 0) {
									tempLogLik *= (1 - (t_prob)matEdgeLiks[i*numNodes+nnode]);
}
							if(tempLogLik > 0)
								tempLogLik = (t_prob)log((double)tempLogLik);
							else 
								tempLogLik = -FLT_MAX;
							if(priorLogLik == -FLT_MAX || tempLogLik == -FLT_MAX)
								priorLogLik = -FLT_MAX;
							else {
								priorLogLik += tempLogLik;
							}
							fLogLik += priorLogLik;
						}

						if(fMaxLogLik < fLogLik) {
							fMaxLogLik = fLogLik;
							ncombMaxLogLik = ncomb;
							pProbNode = baseCatnet.getNodeProb(nnode);
							if(pProbNode) {
								probMaxNode = *pProbNode;
								probMaxNode.priorlik = priorLogLik;
							}
						}
					} /* for ncomb */

					if(ncombMaxLogLik >= 0)
						memcpy(idparset, pcomblist[ncombMaxLogLik], d*sizeof(int));

					/* release combination set */
        				for(ncomb = 0; ncomb < ncomblist; ncomb++) {
        	  				if(pcomblist[ncomb])
        	    					CATNET_FREE(pcomblist[ncomb]);
        	  				pcomblist[ncomb] = NULL;
					} /* for ncomb */
        				CATNET_FREE(pcomblist);
        				pcomblist = 0;
					ncomblist = 0;

					if(ncombMaxLogLik < 0){
						/* retain the same m_pCatnets list */
						continue;
					}

					if(pestim->m_pCacheMutex) {
						MUTEX_LOCK(pestim->m_pCacheMutex);
						setCachedProb(parset, parsetsize + fixparsetsize, 
							nnode, idparset, d, &probMaxNode, fMaxLogLik);
						MUTEX_UNLOCK(pestim->m_pCacheMutex);
					}
					else
						setCachedProb(parset, parsetsize + fixparsetsize, 
							nnode, idparset, d, &probMaxNode, fMaxLogLik);
				} /* if(!getCachedProb) */

				nodecomplx = m_pNodeNumCats[nnode]-1;
				for(k = 0; k < d; k++)
					nodecomplx *= m_pNodeNumCats[idparset[k]];

				for(k = 0; k < m_nCatnets; k++) {
					if(!m_pCatnets[k])
						continue;
					complx = m_pCatnets[k]->complexity() - 
						m_pCatnets[k]->nodeComplexity(nnode) + nodecomplx;
					pProbNode = m_pCatnets[k]->getNodeProb(nnode);
					if(complx > maxComplexity || !pProbNode) 
						continue;
					tempLogLik = m_pCatnets[k]->loglik() - 
							(pProbNode->loglik+pProbNode->priorlik) + fMaxLogLik;
					if(!pCurCatnetList[complx] && tempLogLik > -FLT_MAX) {
#ifdef DISCRETE_SAMPLE
						pCurCatnetList[complx] = new CATNETD<t_prob>;
#else
						pCurCatnetList[complx] = new CATNETP<t_prob>;
#endif
					}
					if(pCurCatnetList[complx] && 
						pCurCatnetList[complx]->loglik() < tempLogLik) {
							*pCurCatnetList[complx] = *m_pCatnets[k];
							pCurCatnetList[complx]->setParents(nnode, idparset, d);
							pCurCatnetList[complx]->setNodeProb(nnode, &probMaxNode);
						}
					}
				} /* for k */

			} /* for d */

			else /*if(!bEqualCategories)*/ {

			for(d = fixparsetsize + 1; d <= maxpars; d++) {
				
				if(_wait_stop_event(0/*millisecs*/) == 0)
					break;

				pcomblist = 0;
				ncomblist = 0;
				_combination_sets<int>(pcomblist, ncomblist, 0, parset, parsetsize, 0, d - fixparsetsize);

				if(fixparsetsize > 0) {
					if(!pcomblist || ncomblist < 1) {
					    	pcomblist = (int**)CATNET_MALLOC(1*sizeof(int*));
						if (!pcomblist)
							return CATNET_ERR_MEM;
			         		pcomblist[0] = 0;	
						ncomblist = 1;
					}
					for(k = 0; k < ncomblist; k++) {
				        	paux = (int*)CATNET_MALLOC(d*sizeof(int));
						if (!paux)
							return CATNET_ERR_MEM;
						for(j = 0; j < fixparsetsize; j++)
				            		paux[j] = fixparset[j];
					        if(pcomblist[k] && d > fixparsetsize) {
				            		memcpy(paux + fixparsetsize, pcomblist[k], (d-fixparsetsize)*sizeof(int));
						}
				            	if(pcomblist[k])
				            		CATNET_FREE(pcomblist[k]); 
				           	pcomblist[k] = paux;
					}
				}
			
				if(becho)
					Rprintf("[%d]%d  ", d, ncomblist);

				fMaxLogLik = -FLT_MAX;
				ncombMaxLogLik = -1;
				probMaxNode.reset();

				for(ncomb = 0; ncomb < ncomblist; ncomb++) {
					nocache = 1;
					if(!pestim->m_pCacheMutex) {
						nocache = !getCachedProb(pcomblist[ncomb], d, nnode, idparset, d, &probMaxNode, &fMaxLogLik);
					}
					else {
						MUTEX_LOCK(pestim->m_pCacheMutex);
						nocache = !getCachedProb(pcomblist[ncomb], d, nnode, idparset, d, &probMaxNode, &fMaxLogLik);
						MUTEX_UNLOCK(pestim->m_pCacheMutex);
					}

					if(nocache) { 
						memcpy(idparset, pcomblist[ncomb], d*sizeof(int));
						// add pcomplist[j] parent set to nnode
						baseCatnet.setParents(nnode, idparset, d);
			     
						// add perturbation
						numSubSamples = numSamples;
						if(perturbations && pSubSamples) {
							numSubSamples = 0;
							for(j = 0; j < numSamples; j++) {
								if(!perturbations[j * numNodes + nnode]) {
#ifdef DISCRETE_SAMPLE
									memcpy(pSubSamples + numSubSamples*numNodes, 
										pSamples + j*numNodes, numNodes*sizeof(int));
#else
									memcpy(pSubSamples + numSubSamples*nline, 
										pSamples + j*nline, nline*sizeof(t_prob));
#endif
									if(pClasses) pSubClasses[numSubSamples] = pClasses[j];
									numSubSamples++;
								}
							}
							if(pClasses)
								fMaxLogLik = baseCatnet.setNodeSampleProbKL(nnode, pSubSamples, numSubSamples, pSubClasses, /*bNormalize = */0, /*bUsePearson = */pestim->m_klmode==1);
							else
								fMaxLogLik = baseCatnet.setNodeSampleProb(nnode, pSubSamples, numSubSamples);
						}
						else {
							if(pClasses)
								fMaxLogLik = baseCatnet.setNodeSampleProbKL(nnode, pSamples, numSamples, pClasses, /*bNormalize = */0, /*bUsePearson = */pestim->m_klmode==1);
							else
								fMaxLogLik = baseCatnet.setNodeSampleProb(nnode, pSamples, numSamples);
						}

						/*prior*/
						priorLogLik = 0;
						if(matEdgeLiks && bernoulibuff) {
							tempLogLik = 1;
							memset(bernoulibuff, 0, numNodes*sizeof(int));
							for(i = 0; i < d; i++) {
								tempLogLik *= (t_prob)matEdgeLiks[pcomblist[ncomb][i]*numNodes+nnode];
								bernoulibuff[pcomblist[ncomb][i]] = 1;
							}
							if(tempLogLik > 0)
								priorLogLik = (t_prob)log((double)tempLogLik);
							else 
								priorLogLik = -FLT_MAX;
							tempLogLik = 1;
							for(i = 0; i < numNodes; i++) 
								if(i != nnode && bernoulibuff[i] == 0)
									tempLogLik *= (1 - (t_prob)matEdgeLiks[i*numNodes+nnode]);
							if(tempLogLik > 0)
								tempLogLik = (t_prob)log((double)tempLogLik);
							else 
								tempLogLik = -FLT_MAX;
							if(priorLogLik == -FLT_MAX || tempLogLik == -FLT_MAX)
								priorLogLik = -FLT_MAX;
							else {
								priorLogLik += tempLogLik;
							}
							fMaxLogLik += priorLogLik;
						}

						pProbNode = baseCatnet.getNodeProb(nnode);
						if(pProbNode)
							probMaxNode = *pProbNode;
						probMaxNode.priorlik = priorLogLik;
						if(pestim->m_pCacheMutex) {
							MUTEX_LOCK(pestim->m_pCacheMutex);
							setCachedProb(idparset, d, nnode, idparset, d,
									&probMaxNode, fMaxLogLik);
							MUTEX_UNLOCK(pestim->m_pCacheMutex);
						}
						else
							setCachedProb(idparset, d, nnode, idparset, d,
									&probMaxNode, fMaxLogLik);
					} /* if(!getCachedProb) */

					/* find nnode-complexity for pcomblist[ncomb] parent set */
					nodecomplx = m_pNodeNumCats[nnode]-1;
					for(k = 0; k < d; k++)
						nodecomplx *= m_pNodeNumCats[idparset[k]];
					for(k = 0; k < m_nCatnets; k++) {
						if(!m_pCatnets[k])
							continue;
						pProbNode = m_pCatnets[k]->getNodeProb(nnode);
						complx = m_pCatnets[k]->complexity() - 
							m_pCatnets[k]->nodeComplexity(nnode) +
							nodecomplx;
						if(complx > maxComplexity || !pProbNode) 
							continue;
						tempLogLik = m_pCatnets[k]->loglik() - 
							(pProbNode->loglik + pProbNode->priorlik) + fMaxLogLik;
						if(!pCurCatnetList[complx] && tempLogLik > -FLT_MAX) {
#ifdef DISCRETE_SAMPLE
							pCurCatnetList[complx] = new CATNETD<t_prob>;
#else
							pCurCatnetList[complx] = new CATNETP<t_prob>;
#endif
						}
						if(pCurCatnetList[complx] && 
							pCurCatnetList[complx]->loglik() < tempLogLik) {
								*pCurCatnetList[complx] = *m_pCatnets[k];
								pCurCatnetList[complx]->setParents(nnode, idparset, d);
								pCurCatnetList[complx]->setNodeProb(nnode, &probMaxNode);
						}
					} /* for k */

				} /* for ncomb */
				/* release combination set */
        			for(ncomb = 0; ncomb < ncomblist; ncomb++) {
          				if(pcomblist[ncomb])
            					CATNET_FREE(pcomblist[ncomb]);
          				pcomblist[ncomb] = NULL;
				} /* for ncomb */
        			CATNET_FREE(pcomblist);
        			pcomblist = 0;
				ncomblist = 0;
			} /* for d */

			} /* if(!bEqualCategories) */

			if(becho)
				Rprintf("\n");

			for(j = 0; j < m_nCatnets; j++) {
				if(m_pCatnets[j]) {
					if(pCurCatnetList[j]) {
						if(m_pCatnets[j]->loglik() < pCurCatnetList[j]->loglik()) {
							delete m_pCatnets[j];
							m_pCatnets[j] = pCurCatnetList[j];
							pCurCatnetList[j] = 0;
						}
						else {
							delete pCurCatnetList[j];
							pCurCatnetList[j] = 0;
						}
					}
				}
				else {
					m_pCatnets[j] = pCurCatnetList[j];
					pCurCatnetList[j] = 0;
				}
			}
		} // for(nnode = 0; nnode < numNodes; nnode++)

		CATNET_FREE(pCurCatnetList);
		CATNET_FREE(parset);
		CATNET_FREE(fixparset);
		CATNET_FREE(idparset);

		if(bernoulibuff)
			CATNET_FREE(bernoulibuff);

		if(pSubSamples)
			CATNET_FREE(pSubSamples);
		if(pSubClasses)
			CATNET_FREE(pSubClasses);

		for(j = 0; j < m_nCatnets; j++) {
			if(m_pCatnets[j]) {
				m_pCatnets[j] -> normalizeProbabilities();
				m_pCatnets[j] -> setCategoryIndices(m_pNodeNumCats, m_pNodeCats);
			}
		}

		if(m_pNodeCats) {
			for(i = 0; i < m_numNodes; i++) 
				if(m_pNodeCats[i])
					CATNET_FREE(m_pNodeCats[i]);
			CATNET_FREE(m_pNodeCats);
		}
		m_pNodeCats = 0;
		if(m_pNodeNumCats) 
			CATNET_FREE(m_pNodeNumCats);
		m_pNodeNumCats = 0;

		return CATNET_ERR_OK;
	}

};

#endif /* CATNET_SEARCH_H */

