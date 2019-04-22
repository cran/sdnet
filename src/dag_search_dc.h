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
 * dag_search_dc.h
 *
 *  Created on: Nov 8, 2011
 *      Author: nbalov
 */

/* Implements dag_search with different node categories 
 */

#include "utils.h"
#include "thread.h"
#include "dag_list.h"
#include "search_params.h"

#if (defined(DISCRETE_SAMPLE) && !defined(DAGD_SEARCH_DC_H)) || (defined(PROB_SAMPLE) && !defined(DAGP_SEARCH_DC_H))

#ifdef DISCRETE_SAMPLE
#define DAGD_SEARCH_DC_H
#else
#ifdef PROB_SAMPLE
#define DAGP_SEARCH_DC_H
#else
#error "No proper DAG class"
#endif
#endif

#ifdef DISCRETE_SAMPLE
template<class t_prob, class t_ind, class t_sample, int n_maxpars>
class DAGD_SEARCH_DC: public DAG_LIST<t_prob, t_ind>, public c_thread {
#else
template<class t_prob, class t_ind, int n_maxpars>
class DAGP_SEARCH_DC: public DAG_LIST<t_prob, t_ind>, public c_thread {
private:
	int m_nline;
	t_prob *m_pProbVect;
#endif

private:
	int m_probSize;
	t_prob *m_prob;
	int *m_pNodeNumCats;
public: 
	using DAG_LIST<t_prob, t_ind>::m_numNodes;
	using DAG_LIST<t_prob, t_ind>::m_numParSlots;
	/* each node has m_numParSlots*n_maxpars of t_ind */
	using DAG_LIST<t_prob, t_ind>::m_parSlots;
	/* each node has n_numParSlots of t_prob */
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
	DAGD_SEARCH_DC() {
		m_prob = 0;
		m_probSize = 0;
		m_pNodeNumCats = 0;
	}
	~DAGD_SEARCH_DC() {
		if(m_pNodeNumCats)
			CATNET_FREE(m_pNodeNumCats);
		m_pNodeNumCats = 0;
		if(m_prob)
			CATNET_FREE(m_prob);
		m_prob = 0;
		m_probSize = 0;
	}
#else
	DAGP_SEARCH_DC() {
		m_prob = 0;
		m_probSize = 0;
		m_pNodeNumCats = 0;
		m_nline = 0;
		m_pProbVect = 0;
	}
	~DAGP_SEARCH_DC() {
		if(m_pNodeNumCats)
			CATNET_FREE(m_pNodeNumCats);
		m_pNodeNumCats = 0;
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
	t_prob nodeLoglik(t_ind nnode, t_sample *pSamples, int& nsamples, t_ind *pars, 
			  int numPars) {

		int i, k, kpad, pad, nProbSize, samp, ncount, nodecats, *pcats;
		/* pSamples have categories in the range [1, nodecats] */
		t_prob loglik, psum;

		if (!m_pNodeNumCats || !m_prob)
			return (t_prob)-FLT_MAX;
		if (nnode < 0 || nnode >= m_numNodes || !pSamples || nsamples < 1)
			return (t_prob)-FLT_MAX;

		nodecats  = m_pNodeNumCats[nnode];
		nProbSize = nodecats;
		for(i = 0; i < numPars; i++)
			nProbSize *= m_pNodeNumCats[pars[i]];
		if(nProbSize > m_probSize)
			return (t_prob)-FLT_MAX;

		pcats = (int*)CATNET_MALLOC(numPars*sizeof(int));
		if (!pcats)
			return (t_prob)-FLT_MAX;
		for (i = 0; i < numPars; i++) 
			pcats[i] = m_pNodeNumCats[pars[i]];
		memset(m_prob, 0, nProbSize*sizeof(t_prob));

		ncount = 0;
		kpad = 0;		
		for (k = 0; k < nsamples; k++) {
			pad = 0;
			for (i = 0; i < numPars; i++) {
				samp = pSamples[kpad + pars[i]];
				if(samp < 1 || samp > pcats[i]) {
					pad = -1;
					break;
				}
				pad =  pcats[i] * pad + samp - 1;
			}
			if (pad < 0) {
				kpad += m_numNodes;
				continue;
			}
			pad *= nodecats;
			if (pad >= nProbSize) {
				kpad += m_numNodes;
				continue;
			}
			samp = pSamples[kpad + nnode]-1;
			if (samp >= 0 && samp < nodecats) {
				m_prob[pad+samp]++;
				ncount++;
			}
			kpad += m_numNodes;
		}
		CATNET_FREE(pcats);
		loglik = 0; k = 0;
		while (k < nProbSize) {
			psum = 0;
			for (i = 0; i < nodecats; i++)
				psum += m_prob[k + i];
			if (psum > 0) {
				psum = 1 / psum;
				for (i = 0; i < nodecats; i++) 
					if(m_prob[k + i] > 0)
						loglik += m_prob[k + i] * log(m_prob[k + i] * psum);
			}
			k += nodecats;
		}
		if(ncount > 1) 
			loglik /= (t_prob)ncount;
		nsamples = ncount;
		return(loglik);
	}

	t_prob nodeLoglikKL(t_ind nnode, t_sample *pSamples, int& nsamples, t_ind *pars, 
			    int numPars, int *pClasses, int bUsePearson = 0) {

		int i, k, kpad, pad, nProbSize, samp, ncount, nodecats, *pcats;
		/* pSamples have categories in the range [1, nodecats] */
		t_prob loglik, psum,  *prob2, psum2;

		if (!m_pNodeNumCats || !m_prob)
			return (t_prob)-FLT_MAX;
		if (nnode < 0 || nnode >= m_numNodes || !pSamples || nsamples < 1)
			return (t_prob)-FLT_MAX;

		nodecats  = m_pNodeNumCats[nnode];
		nProbSize = nodecats;
		for(i = 0; i < numPars; i++)
			nProbSize *= m_pNodeNumCats[pars[i]];
		if(nProbSize > m_probSize)
			return (t_prob)-FLT_MAX;

		pcats = (int*)CATNET_MALLOC(numPars*sizeof(int));
		if (!pcats)
			return (t_prob)-FLT_MAX;
		for (i = 0; i < numPars; i++) 
			pcats[i] = m_pNodeNumCats[pars[i]];

		memset(m_prob, 0, nProbSize*sizeof(t_prob));
		prob2 = (t_prob*)CATNET_MALLOC(nProbSize*sizeof(t_prob));
		if (!prob2) {
			CATNET_FREE(pcats);
			return (t_prob)-FLT_MAX;
		}
		memset(prob2, 0, nProbSize*sizeof(t_prob));

		ncount = 0;
		kpad = 0;		
		for (k = 0; k < nsamples; k++) {
			pad = 0;
			for (i = 0; i < numPars; i++) {
				samp = pSamples[kpad + pars[i]];
				if(samp < 1 || samp > pcats[i]) {
					pad = -1;
					break;
				}
				pad =  pcats[i] * pad + samp - 1;
			}
			if (pad < 0) {
				kpad += m_numNodes;
				continue;
			}
			pad *= nodecats;
			if (pad >= nProbSize) {
				kpad += m_numNodes;
				continue;
			}
			samp = pSamples[kpad + nnode]-1;
			if (samp >= 0 && samp < nodecats) {
				m_prob[pad+samp]++;
				ncount++;
				if(pClasses[k] == 0)
					prob2[pad+samp]++;
			}
			kpad += m_numNodes;
		}
		CATNET_FREE(pcats);
		loglik = 0; k = 0;
		while (k < nProbSize) {
			psum = 0;
			for (i = 0; i < nodecats; i++)
				psum += m_prob[k + i];
			psum2 = 0;
			for (i = 0; i < nodecats; i++)
				psum2 += prob2[k + i];
			if (psum > 0 && psum2 > 0) {
				for (i = 0; i < nodecats; i++) 
					prob2[k + i] /= psum2;
				if(bUsePearson) {
					for (i = 0; i < nodecats; i++) 
						if(prob2[k + i] > 0) {
							psum2 = m_prob[k + i]/psum - prob2[k + i];
							loglik += psum * psum * psum2 * psum2 / m_prob[k + i];
						}
				}
				else {
					for (i = 0; i < nodecats; i++) 
						if(m_prob[k + i] > 0 && prob2[k + i] > 0)
							loglik += m_prob[k + i] * log(m_prob[k + i] / (psum*prob2[k + i]));
				}
			}
			k += nodecats;
		}
		if(ncount > 1) 
			loglik /= (t_prob)ncount;
		nsamples = ncount;
		CATNET_FREE(prob2);
		return(loglik);
	}
#else
	void _unlist_parent_probs(t_ind nnode, t_ind *pnodepars, int nodepars, 
				  t_prob *psample, int &listind, t_prob *plist = 0, 
				  int parid = -1, int parcat = -1, 
				t_prob ptemp = 1, t_prob fact = 0) {

		int i, ncat, sind;
		if(!plist || listind < 0)
			return;
		/* if node probabilities are normalized, call with fact=1 */
		if(fact == 0) {
			int j;
			t_prob aux;
			fact = 1;
			//sind = nnode*n_cats;
			sind = 0;
			for(i = 0; i < nnode; i++)
				sind += m_pNodeNumCats[i];
			aux = 0;
			for(ncat = 0; ncat < m_pNodeNumCats[nnode]; ncat++)
				aux += psample[sind + ncat];
			fact *= aux;
			if(nodepars > 0) { 
				for(j = 0; j < nodepars; j++) {
					//sind = pnodepars[j]*n_cats;
					sind = 0;
					for(i = 0; i < pnodepars[j]; i++)
						sind += m_pNodeNumCats[i];
					aux = 0;
					for(ncat = 0; ncat < m_pNodeNumCats[pnodepars[j]]; ncat++)
						aux += psample[sind + ncat];
					fact *= aux;
				}
			}
			if(fact > 0)
				fact = 1/fact;
		}
		if(nodepars == 0 || (parid == nodepars-1 && parcat >= 0)) {
			//sind = nnode*n_cats;
			sind = 0;
			for(i = 0; i < nnode; i++)
				sind += m_pNodeNumCats[i];
			for(ncat = 0; ncat < m_pNodeNumCats[nnode]; ncat++)
				plist[listind++] = ptemp * psample[sind + ncat] * fact;
			return;
		}
		parid++;
		//sind = pnodepars[parid]*n_cats;
		sind = 0;
		for(i = 0; i < pnodepars[parid]; i++)
			sind += m_pNodeNumCats[i];
		for(parcat = 0; parcat < m_pNodeNumCats[nnode]; parcat++)
			_unlist_parent_probs(nnode, pnodepars, nodepars, psample, 
			listind, plist, parid, parcat, ptemp*psample[sind+parcat], fact);
	};

	t_prob nodeLoglik(t_ind nnode, t_prob *pSamples, int& nsamples, t_ind *pars, int numPars) {

		int i, k, nProbSize, svalid, listind, nodecats, ncount;
		/* pSamples have categories in the range [1, n_cats] */
		t_prob loglik, psum;

		if (!m_pNodeNumCats || !m_prob)
			return (t_prob)-FLT_MAX;
		if (nnode < 0 || nnode >= m_numNodes || !pSamples || nsamples < 1)
			return (t_prob)-FLT_MAX;

		nodecats  = m_pNodeNumCats[nnode];
		nProbSize = nodecats;
		for(i = 0; i < numPars; i++)
			nProbSize *= m_pNodeNumCats[pars[i]];
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
			for (i = 0; i < nodecats; i++)
				psum += m_prob[k + i];
			if (psum > 0) {
				psum = 1 / psum;
				for (i = 0; i < nodecats; i++) 
					if(m_prob[k + i] > 0)
						loglik += m_prob[k + i] * log(m_prob[k + i] * psum);
			}
			k += nodecats;
		}
		if(ncount > 1) 
			loglik /= (t_prob)ncount;
		nsamples = ncount;
		return(loglik);
	}

	t_prob nodeLoglikKL(t_ind nnode, t_prob *pSamples, int& nsamples, t_ind *pars, 
			    int numPars, int *pClasses, int bUsePearson = 0) {

		int i, k, nProbSize, svalid, listind, nodecats, ncount;
		/* pSamples have categories in the range [1, n_cats] */
		t_prob loglik, psum, *prob2, psum2;

		if (!m_pNodeNumCats || !m_prob)
			return (t_prob)-FLT_MAX;
		if (nnode < 0 || nnode >= m_numNodes || !pSamples || nsamples < 1)
			return (t_prob)-FLT_MAX;

		nodecats  = m_pNodeNumCats[nnode];
		nProbSize = nodecats;
		for(i = 0; i < numPars; i++)
			nProbSize *= m_pNodeNumCats[pars[i]];

		memset(m_prob, 0, nProbSize*sizeof(t_prob));
		prob2 = (t_prob*)CATNET_MALLOC(nProbSize*sizeof(t_prob));
		if (!prob2) {
			return (t_prob)-FLT_MAX;
		}
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
			for (i = 0; i < nodecats; i++)
				psum += m_prob[k + i];
			psum2 = 0;
			for (i = 0; i < nodecats; i++)
				psum2 += prob2[k + i];
			if (psum > 0 && psum2 > 0) {
				for (i = 0; i < nodecats; i++) 
					prob2[k + i] /= psum2;
				if(bUsePearson) {
					for (i = 0; i < nodecats; i++) 
						if(prob2[k + i] > 0) {
							psum2 = m_prob[k + i]/psum - prob2[k + i];
							loglik += psum * psum * psum2 * psum2 / m_prob[k + i];
						}
				}
				else {
					for (i = 0; i < nodecats; i++) 
						if(m_prob[k + i] > 0 && prob2[k + i] > 0)
							loglik += m_prob[k + i] * log(m_prob[k + i] / (psum*prob2[k + i]));
				}
			}
			k += nodecats;
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

		int i, j, k, d,nd, ncomb, nnode, nodecomplx, complx, maxCategories;
		int parsetsize, fixparsetsize, maxpars, numSubSamples, maxSubSamples, bSavePars;
		int *paux, **pcomblist, ncomblist, ballow, bfixallow;
		int *pdiffcats, ndiffcats, ndiffpars, *pCombFlags;

		t_ind *parset, *idparset, *fixparset;
		t_prob fMaxLogLik, tempLoglik, fLogLik;
		DAG_PARS<t_prob> *pNewDag, *pDagPars, *pNextDagPars, *pCurDagList, *pCurDag, *pDagAux;

		/* node categories */
		if(m_pNodeNumCats)
			CATNET_FREE(m_pNodeNumCats);
		m_pNodeNumCats = 0;		
		if(!pestim->m_pNodeNumCats)
			return CATNET_ERR_PARAM;
		m_pNodeNumCats = (int*)CATNET_MALLOC(numNodes*sizeof(int));
		if(!m_pNodeNumCats)
			return CATNET_ERR_MEM;
		memcpy(m_pNodeNumCats, pestim->m_pNodeNumCats, numNodes*sizeof(int));

		pdiffcats = (int*)CATNET_MALLOC(numNodes*sizeof(int));
		if(!pdiffcats)
			return CATNET_ERR_MEM;

		maxCategories = 0;
		for(i = 0; i < numNodes; i++) 
			if(maxCategories < m_pNodeNumCats[i])
				maxCategories = m_pNodeNumCats[i];

#ifdef DISCRETE_SAMPLE
		t_sample *pSamples = (t_sample*)pestim->m_pSamples;
		t_sample *pSubSamples;
#else
		int nProbVect;
		t_prob *pSamples = (t_prob*)pestim->m_pSamples;
		t_prob *pSubSamples;
		m_nline = 0;
		for(i = 0; i < numNodes; i++)
			m_nline += m_pNodeNumCats[i];
		nProbVect = maxCategories;
		for(i = 0; i < maxParentSet; i++) 
			nProbVect *= maxCategories;
		if(m_pProbVect)
			CATNET_FREE(m_pProbVect);
		m_pProbVect = (t_prob*) CATNET_MALLOC(nProbVect * sizeof(t_prob));
		if(!m_pProbVect) {
			CATNET_FREE(pdiffcats);
			return CATNET_ERR_MEM;
		}
#endif
		if(numSamples < 1 || !pSamples)
			return CATNET_ERR_PARAM;

		int *pSubClasses = 0;
		int *pClasses = pestim->m_pClasses;
		if(pClasses) {
			pSubClasses = (int*)CATNET_MALLOC(numSamples*sizeof(int));
		}

		m_probSize = maxCategories;
		for(i = 0; i < n_maxpars; i++) {
			m_probSize *= maxCategories;
		}

		m_prob = (t_prob*)CATNET_MALLOC(m_probSize*sizeof(t_prob));
		if(!m_prob) {
			CATNET_FREE(pdiffcats);
			if (pSubClasses)
				CATNET_FREE(pSubClasses);
			return CATNET_ERR_MEM;
		}

		nodecomplx = (maxCategories-1);
		for(k = 0; k < n_maxpars; k++) 
			nodecomplx *= maxCategories;
		if(maxComplexity > numNodes*nodecomplx)
			maxComplexity = numNodes*nodecomplx;

		parset    = (t_ind*)CATNET_MALLOC(numNodes*sizeof(t_ind));
		idparset  = (t_ind*)CATNET_MALLOC(numNodes*sizeof(t_ind));
		fixparset = (t_ind*)CATNET_MALLOC(numNodes*sizeof(t_ind));

		if(!parset || !idparset || !fixparset) {
			CATNET_FREE(pdiffcats);
			if (pSubClasses)
				CATNET_FREE(pSubClasses);
			if (parset)
				CATNET_FREE(parset);
			if(idparset)
				CATNET_FREE(idparset);
			if(fixparset)
				CATNET_FREE(fixparset);
			return CATNET_ERR_MEM;
		}

		pSubSamples = 0;
		if(perturbations)
#ifdef DISCRETE_SAMPLE
			pSubSamples = (t_sample*)CATNET_MALLOC(numNodes*numSamples*sizeof(t_sample));
#else
			pSubSamples =   (t_prob*)CATNET_MALLOC(m_nline*numSamples*sizeof(t_prob));
#endif

		if(m_numNodes != numNodes) {
			/* initialize the first DAG */
			DAG_LIST<t_prob, t_ind>::reset();
			m_numNodes = numNodes;

			m_numParSlots   =   (t_ind*)CATNET_MALLOC(numNodes*sizeof(t_ind));
			m_parSlots      =  (t_ind**)CATNET_MALLOC(numNodes*sizeof(t_ind*));
			m_parLogliks    = (t_prob**)CATNET_MALLOC(numNodes*sizeof(t_prob*));
			m_parComplx     =  (t_ind**)CATNET_MALLOC(numNodes*sizeof(t_ind*));
			m_parSampleSize =  (t_ind**)CATNET_MALLOC(numNodes*sizeof(t_ind*));

			if(!m_numParSlots || !m_parSlots || !m_parLogliks || 
			   !m_parComplx || !m_parSampleSize) {
				CATNET_FREE(pdiffcats);
				if (pSubClasses)
					CATNET_FREE(pSubClasses);
				CATNET_FREE(parset);
				CATNET_FREE(idparset);
				CATNET_FREE(fixparset);
				return CATNET_ERR_MEM;
			}

			memset(m_numParSlots,   0, numNodes*sizeof(t_ind));
			memset(m_parSlots,      0, numNodes*sizeof(t_ind*));
			memset(m_parLogliks,    0, numNodes*sizeof(t_prob*));
			memset(m_parComplx,     0, numNodes*sizeof(t_ind*));
			memset(m_parSampleSize, 0, numNodes*sizeof(t_ind*));

			memset(pdiffcats, 0, numNodes*sizeof(int));
			ndiffcats = 0;
			for(i = 0; i < numNodes; i++) {
				k = 1; 
				if(ndiffcats > 0) {
					for(j = 0; j < ndiffcats; j++) {
						if(pdiffcats[j] == m_pNodeNumCats[i]) k = 0;
					}
				}
				if(k) {
					pdiffcats[ndiffcats] = m_pNodeNumCats[i];
					ndiffcats++;
				}
			}
			if(ndiffcats < 1)
				ndiffcats = 1;
			ndiffpars = 0;
			for(ncomb = 0; ncomb <= n_maxpars; ncomb++) {
				nd = 1; for(j = 0; j < ncomb; j++) nd *= ndiffcats;
				ndiffpars += nd;
			}

			m_dagPars = new DAG_PARS<t_prob>(numNodes, ndiffpars);
			if(!m_dagPars) CATNET_MEM_ERR();
			m_numDags = 1;

			/* set the parents of the base DAG */
			for(nnode = 0; nnode < numNodes; nnode++) {

				memset(pdiffcats, 0, numNodes*sizeof(int));
				ndiffcats = 0;
				for(i = 0; i < nnode; i++) {
					k = 1; 
					if(ndiffcats > 0)
						for(j = 0; j < ndiffcats; j++) if(pdiffcats[j] == m_pNodeNumCats[i]) k = 0;
					if(k) {
						pdiffcats[ndiffcats] = m_pNodeNumCats[i];
						ndiffcats++;
					}
				}
				if(ndiffcats < 1)
					ndiffpars = 1;
				else {
					ndiffpars = 0;
					for(ncomb = 0; ncomb <= n_maxpars && ncomb <= nnode; ncomb++) {
						nd = 1; for(j = 0; j < ncomb; j++) nd *= ndiffcats;
						ndiffpars += nd;
					}
				}

				m_numParSlots[nnode]   = ndiffpars;
				m_parSlots[nnode]      =  (t_ind*)CATNET_MALLOC(ndiffpars*n_maxpars*sizeof(t_ind));
				m_parLogliks[nnode]    = (t_prob*)CATNET_MALLOC(ndiffpars*sizeof(t_prob));
				m_parComplx[nnode]     =  (t_ind*)CATNET_MALLOC(ndiffpars*sizeof(t_ind));
				m_parSampleSize[nnode] =  (t_ind*)CATNET_MALLOC(ndiffpars*sizeof(t_ind));

				if(!m_numParSlots[nnode] || !m_parSlots[nnode] || !m_parLogliks[nnode] || 
			  	   !m_parComplx[nnode] || !m_parSampleSize[nnode]) {
					CATNET_FREE(pdiffcats);
					if (pSubClasses)
						CATNET_FREE(pSubClasses);
					CATNET_FREE(parset);
					CATNET_FREE(idparset);
					CATNET_FREE(fixparset);
					return CATNET_ERR_MEM;
				}

				memset(m_parSlots[nnode],      0, ndiffpars*n_maxpars*sizeof(t_ind));
				memset(m_parLogliks[nnode],    0, ndiffpars*sizeof(t_prob));
				memset(m_parComplx[nnode],     0, ndiffpars*sizeof(t_ind));
				memset(m_parSampleSize[nnode], 0, ndiffpars*sizeof(t_ind));

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

				nodecomplx = (m_pNodeNumCats[nnode]-1);
				for(k = 0; k < fixparsetsize; k++) 
					nodecomplx *= m_pNodeNumCats[fixparset[k]];

				/* find a slot for this node complexity */
				for(nd = 0; nd < ndiffpars; nd++) 
					if(!m_parComplx[nnode][nd] || m_parComplx[nnode][nd] == nodecomplx)
						break;
				m_parComplx[nnode][nd] = nodecomplx;

				if(fixparsetsize > 0)
					// set 1st parent slot 
					memcpy(&m_parSlots[nnode][nd*n_maxpars], fixparset, fixparsetsize*sizeof(t_ind)); 
				for(k = fixparsetsize; k < n_maxpars; k++)
					m_parSlots[nnode][nd*n_maxpars+k] = -1;

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
						m_parLogliks[nnode][nd] = nodeLoglikKL(nnode, pSubSamples, numSubSamples, fixparset, fixparsetsize, pSubClasses, pestim->m_klmode==1);
					else
						m_parLogliks[nnode][nd] = nodeLoglik(nnode, pSubSamples, numSubSamples, fixparset, fixparsetsize);
				}
				else {
					numSubSamples = numSamples;
					if(pClasses)
						m_parLogliks[nnode][nd] = nodeLoglikKL(nnode, pSamples, numSubSamples, fixparset, fixparsetsize, pClasses, pestim->m_klmode==1);
					else
						m_parLogliks[nnode][nd] = nodeLoglik(nnode, pSamples, numSubSamples, fixparset, fixparsetsize);
				}
				m_parSampleSize[nnode][nd] = numSubSamples;

				m_dagPars->loglik += m_parLogliks[nnode][nd];
				m_dagPars->setNumPar(nnode, nd);

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

			if (parset && fixparsetsize > 0)
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
						if (paux) {
							for(j = 0; j < fixparsetsize; j++)
				            			paux[j] = fixparset[j];
					        	if(pcomblist[k] && d > fixparsetsize) {
				            			memcpy(paux + fixparsetsize, pcomblist[k], (d-fixparsetsize)*sizeof(int));
							}
						}
						if(pcomblist[k])
			            			CATNET_FREE(pcomblist[k]); 				           							pcomblist[k] = paux;
					}
				}

				if(becho)
					Rprintf("[%d]%d  ", d, ncomblist);

				pCombFlags = (int*)CATNET_MALLOC(ncomblist*sizeof(int));
				if (!pCombFlags)
					break;
				memset(pCombFlags, 0, ncomblist*sizeof(int));

				while(1) {
					/* find unprocessed par combination */
					ncomb = ncomblist;
					for(i = 0; i < ncomblist; i++) {
						if(pCombFlags[i] == 2) {
							pCombFlags[i] = 1;
							j++;
						}
						if(ncomb >= ncomblist && !pCombFlags[i])
							ncomb = i;
					}
					if(ncomb >= ncomblist)
						break;

					if (pcomblist[ncomb] && d > 0)
						memcpy(idparset, pcomblist[ncomb], d*sizeof(t_ind));
					nodecomplx = (m_pNodeNumCats[nnode]-1);
					for(k = 0; k < d; k++) 
						nodecomplx *= m_pNodeNumCats[idparset[k]];

					/* identify all par combinations with complexity `nodecomplx' */
					pCombFlags[ncomb] = 2;
					if(ncomb < ncomblist-1) {
						for(i = ncomb+1; i < ncomblist; i++) {
							j = (m_pNodeNumCats[nnode]-1);
							for(k = 0; k < d; k++) 
								j *= m_pNodeNumCats[pcomblist[i][k]];
							if(j == nodecomplx)
								pCombFlags[i] = 2;
						}
					}

					fMaxLogLik = -FLT_MAX;
					maxSubSamples = 0;
					i = ncomb;
					for(; i < ncomblist; i++) {
	
						if(pCombFlags[i] != 2)
							continue;
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
								fLogLik = nodeLoglikKL(nnode, pSubSamples, numSubSamples, pcomblist[i], d, pSubClasses, pestim->m_klmode==1);
							else
								fLogLik = nodeLoglik(nnode, pSubSamples, numSubSamples, pcomblist[i], d);
						}
						else {
							numSubSamples = numSamples;
							if(pClasses)
								fLogLik = nodeLoglikKL(nnode, pSamples, numSubSamples, pcomblist[i], d, pClasses, pestim->m_klmode==1);
							else
								fLogLik = nodeLoglik(nnode, pSamples, numSubSamples, pcomblist[i], d);
						}
						if(fLogLik > fMaxLogLik) {
							fMaxLogLik = fLogLik;
							maxSubSamples = numSubSamples;
							ncomb = i;
							if (pcomblist[i] && d > 0)
								memcpy(idparset, pcomblist[i], d*sizeof(t_ind));
						}
					} /* for(; i < ncomblist; i++) */

					/* find a slot for this node complexity */
					for(nd = 0; nd < m_numParSlots[nnode]; nd++) 
						if(!m_parComplx[nnode][nd] || m_parComplx[nnode][nd] == nodecomplx) 
							break;
 
					if(m_parComplx[nnode][nd] == nodecomplx && m_parLogliks[nnode][nd] >= fMaxLogLik) {
						continue; /* next ncomb */
					}

					/* pDagPars are ordered in increasing complexity */
					pDagPars = m_dagPars;
					pCurDag  = pCurDagList;
					bSavePars = 0;
int kk=0;
					while(pCurDag) {
kk++;
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
						pDagAux = pNextDagPars->next;
						while(pDagAux) {
							if(pDagAux->complx >= complx)
								break;
							pNextDagPars = pDagAux;
							pDagAux = pDagAux->next;
						}
						// new DAG is formed upon pCurDag 
						pNewDag = pNextDagPars->next;
						if(pNewDag && pNewDag->complx == complx) {
							if(tempLoglik > pNewDag->loglik) {
								// replace 
								bSavePars = 1;
								pNewDag->copyNumPars(pCurDag);
								pNewDag->setNumPar(nnode, nd); 
								pNewDag->loglik = tempLoglik;
							}
						}
						else {
							// insert a new DAG //
							bSavePars = 1;
							pNewDag = new DAG_PARS<t_prob>(pCurDag);
							if(!pNewDag) CATNET_MEM_ERR();
							pNewDag->next = pNextDagPars->next;
							pNewDag->setNumPar(nnode, nd);
							pNewDag->complx = complx;
							pNewDag->loglik = tempLoglik;
							pNextDagPars->next = pNewDag;
							m_numDags++;
						}
						pCurDag = pCurDag->next;
						pDagPars = pDagPars->next;
					} // pDagPars

					if(bSavePars) {
						// save the d-parent set result
						m_parLogliks[nnode][nd] = fMaxLogLik; 
						m_parComplx[nnode][nd] = nodecomplx; 
						m_parSampleSize[nnode][nd] = maxSubSamples;
						if (idparset && d > 0)
							memcpy(&m_parSlots[nnode][nd*n_maxpars], idparset, d*sizeof(t_ind)); 
						for(k = d; k < n_maxpars; k++)
							m_parSlots[nnode][nd*n_maxpars+k] = -1;
					}

				} /* while(1) */

				CATNET_FREE(pCombFlags);
				pCombFlags = 0;

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
		
			/* resize */
			for(nd = 0; nd < m_numParSlots[nnode]; nd++) 
				if(!m_parComplx[nnode][nd])
					break;
			m_numParSlots[nnode] = nd;

			/* release the copy */
			delete pCurDagList;
			pCurDagList = 0;

			if(becho)
				Rprintf("\n");

		} // for(nnode = 0; nnode < numNodes; nnode++)

		CATNET_FREE(pdiffcats);
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
		if(m_pNodeNumCats)
			CATNET_FREE(m_pNodeNumCats);
		m_pNodeNumCats = 0;	

		return CATNET_ERR_OK;
	}

};

#endif /* DAG_SEARCH_DC_H */

