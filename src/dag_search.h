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
 * dag_search.h
 *
 *  Created on: Nov 8, 2011
 *      Author: nbalov
 */

#include "utils.h"
#include "thread.h"
#include "dag_list.h"
#include "search_params.h"

#if (defined(DISCRETE_SAMPLE) && !defined(DAGD_SEARCH_H)) || (defined(PROB_SAMPLE) && !defined(DAGP_SEARCH_H))

#ifdef DISCRETE_SAMPLE
#define DAGD_SEARCH_H
#else
#ifdef PROB_SAMPLE
#define DAGP_SEARCH_H
#else
#error "No proper DAG class"
#endif
#endif

#ifdef DISCRETE_SAMPLE
template<class t_prob, class t_ind, class t_sample, int n_maxpars, int n_cats>
class DAGD_SEARCH: public DAG_LIST<t_prob, t_ind>, public c_thread {
#else
template<class t_prob, class t_ind, int n_maxpars, int n_cats>
class DAGP_SEARCH: public DAG_LIST<t_prob, t_ind>, public c_thread {
private:
	int m_nline;
	t_prob *m_pProbVect;
#endif

private:
	int m_probSize;
	t_prob *m_prob;
public: 
	using DAG_LIST<t_prob, t_ind>::m_numNodes;
	using DAG_LIST<t_prob, t_ind>::m_numParSlots;
	/* each node has (n_maxpars+1)*n_maxpars of t_ind */
	using DAG_LIST<t_prob, t_ind>::m_parSlots;
	/* each node has n_maxpars of t_prob */
	using DAG_LIST<t_prob, t_ind>::m_parLogliks;
	using DAG_LIST<t_prob, t_ind>::m_parComplx;
	using DAG_LIST<t_prob, t_ind>::m_parSampleSize;
	using DAG_LIST<t_prob, t_ind>::m_numDags;
	using DAG_LIST<t_prob, t_ind>::m_dagPars; 

	/*void *operator new(size_t size) {
		return(CATNET_MALLOC(size));
	}

	void operator delete(void *pobj) {
		CATNET_FREE(pobj);
	}*/

#ifdef DISCRETE_SAMPLE
	DAGD_SEARCH() {
		m_prob = 0;
		m_probSize = 0;
	}
	~DAGD_SEARCH() {
		if(m_prob)
			CATNET_FREE(m_prob);
		m_prob = 0;
		m_probSize = 0;
	}
#else
	DAGP_SEARCH() {
		m_prob = 0;
		m_probSize = 0;
		m_nline = 0;
		m_pProbVect = 0;
	}
	~DAGP_SEARCH() {
		if(m_prob)
			CATNET_FREE(m_prob);
		m_prob = 0;
		if(m_pProbVect)
			CATNET_FREE(m_pProbVect);
		m_pProbVect = 0;
	}
#endif

	// sets sample conditional probability and returns its log-likelihood
#ifdef DISCRETE_SAMPLE
	t_prob nodeLoglik(t_ind nnode, t_sample *pSamples, int& nsamples, t_ind *pars, int numPars) {

		int i, k, pad, nProbSize, samp, ncount;
		/* pSamples have categories in the range [1, n_cats] */
		t_prob loglik, psum;

		if (!m_prob || nnode < 0 || nnode >= m_numNodes || !pSamples || nsamples < 1)
			return (t_prob)-FLT_MAX;

		nProbSize = n_cats;
		for(i = 0; i < numPars; i++)
			nProbSize *= n_cats;
		if(nProbSize > m_probSize)
			return (t_prob)-FLT_MAX;
		memset(m_prob, 0, nProbSize*sizeof(t_prob));

		ncount = 0;
		for (k = 0; k < nsamples; k++) {
			pad = 0;
			for (i = 0; i < numPars; i++) {
				samp = pSamples[k * m_numNodes + pars[i]];
				if(samp < 1 || samp > n_cats) {
					pad = -1;
					break;
				}
				pad = n_cats * pad + samp - 1;
			}
			pad *= n_cats;
			if (pad < 0 || pad >= nProbSize)
				continue;
			samp = pSamples[k * m_numNodes + nnode]-1;
			if (samp >= 0 && samp < n_cats) {
				m_prob[pad+samp]++;
				ncount++;
			}
		}
		loglik = 0; k = 0;
		while (k < nProbSize) {
			psum = 0;
			for (i = 0; i < n_cats; i++)
				psum += m_prob[k + i];
			if (psum > 0) {
				psum = 1 / psum;
				for (i = 0; i < n_cats; i++) 
					if(m_prob[k + i] > 0)
						loglik += m_prob[k + i] * log(m_prob[k + i] * psum);
			}
			k += n_cats;
		}
		if(ncount > 1) 
			loglik /= (t_prob)ncount;
		nsamples = ncount;
		return(loglik);
	}

	t_prob nodeLoglikKL(t_ind nnode, t_sample *pSamples, int& nsamples, t_ind *pars, 
			    int numPars, int *pClasses, int bUsePearson = 0) {

		int i, k, pad, nProbSize, samp, ncount;
		/* pSamples have categories in the range [1, n_cats] */
		t_prob loglik, psum, *prob2, psum2;
		if (nnode < 0 || nnode >= m_numNodes || !pSamples || nsamples < 1)
			return (t_prob)-FLT_MAX;

		nProbSize = n_cats;
		for(i = 0; i < numPars; i++)
			nProbSize *= n_cats;
		if(!m_prob || nProbSize > m_probSize)
			return (t_prob)-FLT_MAX;
		memset(m_prob, 0, nProbSize*sizeof(t_prob));
		prob2 = (t_prob*)CATNET_MALLOC(nProbSize*sizeof(t_prob));
		if (!prob2)
			return (t_prob)-FLT_MAX;
		memset(prob2, 0, nProbSize*sizeof(t_prob));

		ncount = 0;
		for (k = 0; k < nsamples; k++) {
			pad = 0;
			for (i = 0; i < numPars; i++) {
				samp = pSamples[k * m_numNodes + pars[i]];
				if(samp < 1 || samp > n_cats) {
					pad = -1;
					break;
				}
				pad = n_cats * pad + samp - 1;
			}
			pad *= n_cats;
			if (pad < 0 || pad >= nProbSize)
				continue;
			samp = pSamples[k * m_numNodes + nnode]-1;
			if (samp >= 0 && samp < n_cats) {
				m_prob[pad+samp]++;
				ncount++;
				if(pClasses[k] == 0)
					prob2[pad+samp]++;
			}
		}
		loglik = 0; k = 0;
		while (k < nProbSize) {
			psum = 0;
			for (i = 0; i < n_cats; i++)
				psum += m_prob[k + i];
			psum2 = 0;
			for (i = 0; i < n_cats; i++)
				psum2 += prob2[k + i];
			if (psum > 0 && psum2 > 0) {
				for (i = 0; i < n_cats; i++) 
					prob2[k + i] /= psum2;
				if(bUsePearson) {
					for (i = 0; i < n_cats; i++) 
						if(prob2[k + i] > 0) {
							psum2 = m_prob[k + i]/psum - prob2[k + i];
							loglik += psum * psum * psum2 * psum2 / m_prob[k + i];
						}
				}
				else {
					for (i = 0; i < n_cats; i++) 
						if(m_prob[k + i] > 0 && prob2[k + i] > 0)
							loglik += m_prob[k + i] * log(m_prob[k + i] / (psum*prob2[k + i]));
				}
			}
			k += n_cats;
		}
		if(ncount > 1) 
			loglik /= (t_prob)ncount;
		nsamples = ncount;
		CATNET_FREE(prob2);
		return(loglik);
	}
#else
	void _unlist_parent_probs(t_ind nnode, t_ind *pnodepars, int nodepars, t_prob *psample, 
				int &listind, t_prob *plist = 0, int parid = -1, int parcat = -1, 
				t_prob ptemp = 1, t_prob fact = 0) {

		int ncat, sind;
		if(!plist || listind < 0)
			return;

		/* if node probabilities are normalized, call with fact=1 */
		if(fact == 0) {
			int j;
			t_prob aux;
			fact = 1;
			sind = nnode*n_cats;
			//for(i = 0; i < nnode; i++)
			//	sind += m_numCategories[i];
			aux = 0;
			for(ncat = 0; ncat < n_cats; ncat++)
				aux += psample[sind + ncat];
			fact *= aux;
			if(nodepars > 0) { 
				for(j = 0; j < nodepars; j++) {
					sind = pnodepars[j]*n_cats;
					//for(i = 0; i < pnodepars[j]; i++)
					//	sind += m_numCategories[i];
					aux = 0;
					for(ncat = 0; ncat < n_cats; ncat++)
						aux += psample[sind + ncat];
					fact *= aux;
				}
			}
			if(fact > 0)
				fact = 1/fact;
		}
		if(nodepars == 0 || (parid == nodepars-1 && parcat >= 0 && parcat < n_cats)) {
			sind = nnode*n_cats;
			//for(i = 0; i < nnode; i++)
			//	sind += m_numCategories[i];
			for(ncat = 0; ncat < n_cats; ncat++)
				plist[listind++] = ptemp * psample[sind + ncat] * fact;
			return;
		}
		parid++;
		sind = pnodepars[parid]*n_cats;
		//for(i = 0; i < pnodepars[parid]; i++)
		//	sind += m_numCategories[i];
		for(parcat = 0; parcat < n_cats; parcat++)
			_unlist_parent_probs(nnode, pnodepars, nodepars, psample, listind, 
					     plist, parid, parcat, ptemp*psample[sind+parcat], fact);
	};

	t_prob nodeLoglik(t_ind nnode, t_prob *pSamples, int& nsamples, t_ind *pars, int numPars) {

		int i, k, nProbSize, svalid, listind, ncount;

		/* pSamples have categories in the range [1, n_cats] */
		t_prob loglik, psum;
		if (!m_prob || nnode < 0 || nnode >= m_numNodes || !pSamples || nsamples < 1)
			return (t_prob)-FLT_MAX;

		nProbSize = n_cats;
		for(i = 0; i < numPars; i++)
			nProbSize *= n_cats;
		memset(m_prob, 0, nProbSize*sizeof(t_prob));

		ncount = 0;
		for (k = 0; k < nsamples; k++) {
			listind = 0;
			_unlist_parent_probs(nnode, pars, numPars, pSamples+m_nline*k, 
					     listind, m_pProbVect, -1, -1, 1, 1);
			svalid = 1;
			for(i = 0; i < nProbSize; i++)
				if(m_pProbVect[i] < 0) /* NA */ {
					svalid = 0;
					break;
				}
			if(!svalid)
				continue;
			for(i = 0; i < nProbSize; i++)
				m_prob[i] += m_pProbVect[i];
			ncount++;
		}
		loglik = 0; k = 0;
		while (k < nProbSize) {
			psum = 0;
			for (i = 0; i < n_cats; i++)
				psum += m_prob[k + i];
			if (psum > 0) {
				psum = 1 / psum;
				for (i = 0; i < n_cats; i++) 
					if(m_prob[k + i] > 0)
						loglik += m_prob[k + i] * log(m_prob[k + i] * psum);
			}
			k += n_cats;
		}
		if(ncount > 1) 
			loglik /= (t_prob)ncount;
		nsamples = ncount;
		return(loglik);
	}
	
	t_prob nodeLoglikKL(t_ind nnode, t_prob *pSamples, int& nsamples, t_ind *pars, 
			    int numPars, int *pClasses, int bUsePearson = 0) {

		int i, k, nProbSize, svalid, listind, ncount;

		/* pSamples have categories in the range [1, n_cats] */
		t_prob loglik, psum, *prob2, psum2;
		if (!m_prob || nnode < 0 || nnode >= m_numNodes || !pSamples || nsamples < 1)
			return (t_prob)-FLT_MAX;

		nProbSize = n_cats;
		for(i = 0; i < numPars; i++)
			nProbSize *= n_cats;
		memset(m_prob, 0, nProbSize*sizeof(t_prob));
		prob2 = (t_prob*)CATNET_MALLOC(nProbSize*sizeof(t_prob));
		if (!prob2)
			return (t_prob)-FLT_MAX;
		memset(prob2, 0, nProbSize*sizeof(t_prob));

		ncount = 0;
		for (k = 0; k < nsamples; k++) {
			listind = 0;
			_unlist_parent_probs(nnode, pars, numPars, pSamples+m_nline*k, 
					     listind, m_pProbVect, -1, -1, 1, 1);
			svalid = 1;
			for(i = 0; i < nProbSize; i++)
				if(m_pProbVect[i] < 0) /* NA */ {
					svalid = 0;
					break;
				}
			if(!svalid)
				continue;
			for(i = 0; i < nProbSize; i++) 
				m_prob[i] += m_pProbVect[i];
			if(pClasses[k] == 0)
				for(i = 0; i < nProbSize; i++) 
					prob2[i] += m_pProbVect[i];
			ncount++;
		}
		loglik = 0; k = 0;
		while (k < nProbSize) {
			psum = 0;
			for (i = 0; i < n_cats; i++)
				psum += m_prob[k + i];
			psum2 = 0;
			for (i = 0; i < n_cats; i++)
				psum2 += prob2[k + i];
			if (psum > 0 && psum2 > 0) {
				for (i = 0; i < n_cats; i++) 
					prob2[k + i] /= psum2;
				if(bUsePearson) {
					for (i = 0; i < n_cats; i++) 
						if(prob2[k + i] > 0) {
							psum2 = m_prob[k + i]/psum - prob2[k + i];
							loglik += psum * psum * psum2 * psum2 / m_prob[k + i];
						}
				}
				else {
					for (i = 0; i < n_cats; i++) 
						if(m_prob[k + i] > 0 && prob2[k + i] > 0)
							loglik += m_prob[k + i] * log(m_prob[k + i] / (psum*prob2[k + i]));
				}
			}
			k += n_cats;
		}
		if(ncount > 1) 
			loglik /= (t_prob)ncount;
		nsamples = ncount;
		CATNET_FREE(prob2);
		return(loglik);
	}
#endif

public:
	/* Each parentsPool[i] is numNodes long ! */
	int search(SEARCH_PARAMETERS *pestim) {

		if(!pestim)
			return 0;
		int numNodes = pestim->m_numNodes;
		if(numNodes < 1)
			return CATNET_ERR_PARAM;
		int numSamples = pestim->m_numSamples;
		int *perturbations = pestim->m_pPerturbations;
		int maxParentSet = pestim->m_maxParentSet;
		int *parSizes = pestim->m_pParentSizes;
		int maxComplexity = pestim->m_maxComplexity;
		int **parentsPool = pestim->m_parentsPool;
		int **fixedParentsPool = pestim->m_fixedParentsPool;
		int becho = pestim->m_echo;

		int i, j, k, d, ncomb, ncombMaxLogLik, nnode, nodecomplx, complx;
		int parsetsize, fixparsetsize, maxpars, numSubSamples, maxSubSamples;
		int *paux, **pcomblist, ncomblist, ballow, bfixallow;
	
		t_ind *parset, *idparset, *fixparset;
		t_prob fLogLik, fMaxLogLik, tempLoglik;
		DAG_PARS<t_prob> *pNewDag, *pDagPars, *pNextDagPars, *pCurDagList, *pCurDag;

#ifdef DISCRETE_SAMPLE
		t_sample *pSamples = (t_sample*)pestim->m_pSamples;
		t_sample *pSubSamples;
#else
		int nProbVect;
		t_prob *pSamples = (t_prob*)pestim->m_pSamples;
		t_prob *pSubSamples;
		m_nline = numNodes*n_cats;

		if(m_pProbVect)
			CATNET_FREE(m_pProbVect);
		nProbVect = n_cats;
		for(i = 0; i < maxParentSet; i++) 
			nProbVect *= n_cats;
		m_pProbVect = (t_prob*) CATNET_MALLOC(nProbVect * sizeof(t_prob));
		if (!m_pProbVect)
			return CATNET_ERR_MEM;			
#endif
		if(numSamples < 1 || !pSamples)
			return CATNET_ERR_PARAM;

		int *pSubClasses = 0;
		int *pClasses = pestim->m_pClasses;
		if(pClasses) {
			pSubClasses = (int*)CATNET_MALLOC(numSamples*sizeof(int));
			if (!pSubClasses)
				return CATNET_ERR_MEM;
		}

		m_probSize = n_cats;
		for(i = 0; i < n_maxpars; i++)
			m_probSize *= n_cats;
		m_prob = (t_prob*)CATNET_MALLOC(m_probSize*sizeof(t_prob));
		if (!m_prob)
			return CATNET_ERR_MEM;

		nodecomplx = (n_cats-1);
		for(k = 0; k < n_maxpars; k++) 
			nodecomplx *= n_cats;
		if(maxComplexity < numNodes*nodecomplx)
			maxComplexity = numNodes*nodecomplx;

		parset    = (t_ind*)CATNET_MALLOC(numNodes*sizeof(t_ind));
		idparset  = (t_ind*)CATNET_MALLOC(numNodes*sizeof(t_ind));
		fixparset = (t_ind*)CATNET_MALLOC(numNodes*sizeof(t_ind));
		if (!parset || !idparset || !fixparset) {
			if (pSubClasses)
				CATNET_FREE(pSubClasses);
			if (parset)
				CATNET_FREE(parset);
			if (fixparset)
				CATNET_FREE(fixparset);
			if (idparset)
				CATNET_FREE(idparset);
			return CATNET_ERR_MEM;
		}

		pSubSamples = 0;
		if(perturbations) 
#ifdef DISCRETE_SAMPLE
			pSubSamples = (t_sample*)CATNET_MALLOC(numNodes*numSamples*sizeof(t_sample));
#else
			pSubSamples = (t_prob*)CATNET_MALLOC(m_nline*numSamples*sizeof(t_prob));
#endif

		if(m_numNodes != numNodes) {
			/* initialize the first DAG */
			DAG_LIST<t_prob, t_ind>::reset();
			m_numNodes = numNodes;
			m_numParSlots   = (t_ind*)CATNET_MALLOC(numNodes*sizeof(t_ind));
			m_parSlots      = (t_ind**)CATNET_MALLOC(numNodes*sizeof(t_ind*));
			m_parLogliks    = (t_prob**)CATNET_MALLOC(numNodes*sizeof(t_prob*));
			m_parComplx     = (t_ind**)CATNET_MALLOC(numNodes*sizeof(t_ind*));
			m_parSampleSize = (t_ind**)CATNET_MALLOC(numNodes*sizeof(t_ind*));
			if (!m_numParSlots || !m_parSlots || !m_parLogliks || 
		            !m_parComplx   || !m_parSampleSize) {
				if (pSubClasses)
					CATNET_FREE(pSubClasses);
				CATNET_FREE(parset);
				CATNET_FREE(fixparset);
				CATNET_FREE(idparset);
				if (pSubSamples)
					CATNET_FREE(pSubSamples);
				return CATNET_ERR_MEM;
			}

			memset(m_parSlots,      0, numNodes*sizeof(t_ind*));
			memset(m_parLogliks,    0, numNodes*sizeof(t_prob*));
			memset(m_parComplx,     0, numNodes*sizeof(t_ind*));
			memset(m_parSampleSize, 0, numNodes*sizeof(t_ind*));

			for(nnode = 0; nnode < numNodes; nnode++) {
				m_numParSlots[nnode] = n_maxpars+1;
				m_parSlots[nnode]      = (t_ind*)CATNET_MALLOC((n_maxpars+1)*n_maxpars*sizeof(t_ind));
				m_parLogliks[nnode]    = (t_prob*)CATNET_MALLOC((n_maxpars+1)*sizeof(t_prob));
				m_parComplx[nnode]     = (t_ind*)CATNET_MALLOC((n_maxpars+1)*sizeof(t_ind));
				m_parSampleSize[nnode] = (t_ind*)CATNET_MALLOC((n_maxpars+1)*sizeof(t_ind));
				if (!m_numParSlots[nnode] || !m_parSlots[nnode] || !m_parLogliks[nnode] || 
				    !m_parComplx[nnode]   || !m_parSampleSize[nnode]) {
					return CATNET_ERR_MEM;
				}
				memset(m_parSlots[nnode],      0, (n_maxpars+1)*n_maxpars*sizeof(t_ind));
				memset(m_parLogliks[nnode],    0, (n_maxpars+1)*sizeof(t_prob));
				memset(m_parComplx[nnode],     0, (n_maxpars+1)*sizeof(t_ind));
				memset(m_parSampleSize[nnode], 0, (n_maxpars+1)*sizeof(t_ind));
			}
			m_dagPars = new DAG_PARS<t_prob>(numNodes, n_maxpars);
			if(!m_dagPars) CATNET_MEM_ERR();
			m_numDags = 1;

			/* set the parents of the base DAG */
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
				}
				if(fixparsetsize > n_maxpars) {
					CATNET_ERR("fixparsetsize > n_maxpars");
				}
				// set 1st parent slot 
				if(fixparsetsize > 0) 
					memcpy(&m_parSlots[nnode][fixparsetsize*n_maxpars], fixparset, fixparsetsize*sizeof(t_ind)); 
				for(k = fixparsetsize; k < n_maxpars; k++)
					m_parSlots[nnode][fixparsetsize*n_maxpars+k] = -1;

				// calculate the log-likelihood
				if(perturbations && pSubSamples) {
					numSubSamples = 0;
					for(j = 0; j < numSamples; j++) {
						if(!perturbations[j * numNodes + nnode]) {
#ifdef DISCRETE_SAMPLE
							memcpy(pSubSamples + numSubSamples*numNodes, pSamples + j*numNodes, numNodes*sizeof(t_sample));
#else
							memcpy(pSubSamples + numSubSamples*m_nline, pSamples + j*m_nline, m_nline*sizeof(t_prob));
#endif
							if(pClasses) pSubClasses[numSubSamples] = pClasses[j];
							numSubSamples++;
						}
					}
					if(pClasses)
						m_parLogliks[nnode][fixparsetsize] = nodeLoglikKL(nnode, pSubSamples, numSubSamples, fixparset, fixparsetsize, pSubClasses, pestim->m_klmode==1);
					else
						m_parLogliks[nnode][fixparsetsize] = nodeLoglik(nnode, pSubSamples, numSubSamples, fixparset, fixparsetsize);
				}
				else {
					numSubSamples = numSamples;
					if(pClasses)
						m_parLogliks[nnode][fixparsetsize] = nodeLoglikKL(nnode, pSamples, numSubSamples, fixparset, fixparsetsize, pClasses, pestim->m_klmode==1);
					else
						m_parLogliks[nnode][fixparsetsize] = nodeLoglik(nnode, pSamples, numSubSamples, fixparset, fixparsetsize);
				}
				m_parSampleSize[nnode][fixparsetsize] = numSubSamples;

				nodecomplx = (n_cats-1);
				for(k = 0; k < fixparsetsize; k++) 
					nodecomplx *= n_cats;
				m_parComplx[nnode][fixparsetsize] = nodecomplx;

				m_dagPars->loglik += m_parLogliks[nnode][fixparsetsize];
				m_dagPars->setNumPar(nnode,fixparsetsize);

				m_dagPars->complx += nodecomplx;
			}
		} /* m_numNodes != numNodes */

		/* main loop of consequential non-empty-parenthood-node additions */
		for(nnode = 0/*1*/; nnode < numNodes; nnode++) {

			if(_wait_stop_event(4/*millisecs*/) == 0) {
				if(becho)
					Rprintf("STOP signal detected\n");
				break;
			}

			fixparsetsize = 0;
			parsetsize = 0;
			for(j = 0; j < numNodes/*nnode*/; j++) {
				if(j == nnode)
					continue;
				if((!parentsPool || !parentsPool[nnode]) && j > nnode)
					continue;
				ballow = 1;
				if(parentsPool && !parentsPool[nnode])
					ballow = 0;
				if(parentsPool && parentsPool[nnode]) {
					ballow = 0;
					for(k = 0; k < pestim->m_maxParentsPool; k++) {
						if(parentsPool[nnode][k] < 0)
							break;
						if(j == parentsPool[nnode][k])
							ballow = 1;
					}
				}
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
			if (fixparsetsize > 0)
				memcpy(parset + parsetsize, fixparset, fixparsetsize*sizeof(t_ind));

			maxpars = maxParentSet;
			if(parSizes && parSizes[nnode] < maxParentSet)
				maxpars = parSizes[nnode];

			if(maxpars > parsetsize + fixparsetsize)
				maxpars = parsetsize + fixparsetsize;
			if(maxpars > n_maxpars) {
				CATNET_WARNING("maxpars > n_maxpars");
				maxpars = n_maxpars;
			}

			if(becho) {
				Rprintf("processing node %d\n", nnode+1);
				Rprintf("    [#parents][#combinations] = ");
			}

			/* make a copy of the current dag list */
			pCurDagList = 0;
			pCurDag = pCurDagList;
			pDagPars = m_dagPars;
			while(pDagPars) {
				pNewDag = new DAG_PARS<t_prob>(pDagPars);
				if(pCurDag)
					pCurDag->next = pNewDag;
				pCurDag = pNewDag;
				if(!pCurDagList)
					pCurDagList = pNewDag;
				pDagPars = pDagPars->next;
			}

			for(d = fixparsetsize + 1; d <= maxpars; d++) {

				pcomblist = 0;
				ncomblist = 0;
				_combination_sets<t_ind>(pcomblist, ncomblist, 0, parset, parsetsize, 0, d - fixparsetsize);
				if (!pcomblist)
					break;

			        if(fixparsetsize > 0) {
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
				maxSubSamples = 0;

				for(ncomb = 0; ncomb < ncomblist; ncomb++) {
					if(perturbations && pSubSamples) {
						numSubSamples = 0;
						for(j = 0; j < numSamples; j++) {
							if(!perturbations[j * numNodes + nnode]) {
#ifdef DISCRETE_SAMPLE
								memcpy(pSubSamples + numSubSamples*numNodes, pSamples + j*numNodes, numNodes*sizeof(t_sample));
#else
								memcpy(pSubSamples + numSubSamples*m_nline, pSamples + j*m_nline, m_nline*sizeof(t_prob));
#endif
								if(pClasses) pSubClasses[numSubSamples] = pClasses[j];
								numSubSamples++;
							}
						}
						if(pClasses)
							fLogLik = nodeLoglikKL(nnode, pSubSamples, numSubSamples, pcomblist[ncomb], d, pSubClasses, pestim->m_klmode==1);
						else
							fLogLik = nodeLoglik(nnode, pSubSamples, numSubSamples, pcomblist[ncomb], d);
					}
					else {
						numSubSamples = numSamples;
						if(pClasses) 
							fLogLik = nodeLoglikKL(nnode, pSamples, numSubSamples, pcomblist[ncomb], d, pClasses, pestim->m_klmode==1);
						else
							fLogLik = nodeLoglik(nnode, pSamples, numSubSamples, pcomblist[ncomb], d);
					}
					if(fMaxLogLik < fLogLik) {
						fMaxLogLik = fLogLik;
						ncombMaxLogLik = ncomb;
						maxSubSamples = numSubSamples;
					}
				} /* for ncomb */
				if(ncombMaxLogLik < 0)
					continue;

				if (pcomblist[ncombMaxLogLik] && d > 0)
					memcpy(idparset, pcomblist[ncombMaxLogLik], d*sizeof(t_ind));
				/* release combination set */
        			for(ncomb = 0; ncomb < ncomblist; ncomb++) {
        	  			if(pcomblist[ncomb])
        	    				CATNET_FREE(pcomblist[ncomb]);
        	  			pcomblist[ncomb] = NULL;
				} /* for ncomb */

        			CATNET_FREE(pcomblist);
        			pcomblist = 0;
				ncomblist = 0;

				nodecomplx = (n_cats-1);
				for(k = 0; k < d; k++) 
					nodecomplx *= (int)n_cats;
				/* save the d-parent set result */
				m_parLogliks[nnode][d] = fMaxLogLik; 
				m_parComplx[nnode][d] = nodecomplx;
				m_parSampleSize[nnode][d] = maxSubSamples; 
				if (d > 0)
					memcpy(&m_parSlots[nnode][d*n_maxpars], idparset, d*sizeof(t_ind)); 
				for(k = d; k < n_maxpars; k++)
					m_parSlots[nnode][d*n_maxpars+k] = -1;

				/* pDagPars are ordered in increasing complexity */
				pDagPars = m_dagPars;
				pCurDag = pCurDagList;
				while(pCurDag) {
					i = pCurDag->getNumPars(nnode);
					j = nodecomplx - m_parComplx[nnode][i];
					if(j <= 0) {
						pCurDag = pCurDag->next;
						pDagPars = pDagPars->next;
						continue;
					}
					complx = pCurDag->complx + j;
					if(complx > maxComplexity) {
						pCurDag = pCurDag->next;
						pDagPars = pDagPars->next;
						continue;
					}
					tempLoglik = pCurDag->loglik - m_parLogliks[nnode][i] + fMaxLogLik;

					pNextDagPars = pDagPars;
					while(pNextDagPars->next) {
						if(pNextDagPars->next->complx >= complx)
							break;
						pNextDagPars = pNextDagPars->next;
					}
					/* new DAG is formed upon pCurDag */
					pNewDag = pNextDagPars->next;
					if(pNewDag && pNewDag->complx == complx) {
						if(tempLoglik > pNewDag->loglik) {
							/* replace */
							pNewDag->copyNumPars(pCurDag);
							pNewDag->setNumPar(nnode, d); 
							pNewDag->loglik = tempLoglik;
						}
					}
					else {
						/* insert a new DAG */
						pNewDag = new DAG_PARS<t_prob>(pCurDag);
						if(!pNewDag) 
							CATNET_MEM_ERR();
						pNewDag->next = pNextDagPars->next;
						pNewDag->setNumPar(nnode, d);
						pNewDag->complx = complx;
						pNewDag->loglik = tempLoglik;
						pNextDagPars->next = pNewDag;
						m_numDags++;
					}

					pCurDag = pCurDag->next;
					pDagPars = pDagPars->next;
				} /* pDagPars */
				
			} /* for d */
		
			/* release the copy */
			delete pCurDagList;
			pCurDagList = 0;

			if(becho)
				Rprintf("\n");

		} // for(nnode = 0; nnode < numNodes; nnode++)


		CATNET_FREE(parset);
		CATNET_FREE(fixparset);
		CATNET_FREE(idparset);

		if(pSubSamples)
			CATNET_FREE(pSubSamples);
		if(pSubClasses)
			CATNET_FREE(pSubClasses);

		CATNET_FREE(m_prob);
		m_prob = 0;
		m_probSize = 0;

#ifndef DISCRETE_SAMPLE
		if(m_pProbVect)
			CATNET_FREE(m_pProbVect);
		m_pProbVect = 0;
#endif

		return CATNET_ERR_OK;
	}

};

#endif /* DAG_SEARCH_H */

