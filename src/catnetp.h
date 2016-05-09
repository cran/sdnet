/*  catnet : categorical Bayesian network inference
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
 * catnetp.h
 *
 *  Created on: Sep 18, 2011
 *      Author: nbalov
 */

#include "catnet_base.h"

#ifndef CATNETP_H_
#define CATNETP_H_

template<class t_prob>
class CATNETP: public CATNET<t_prob, t_prob> {

protected:
	using CATNET<t_prob,t_prob>::m_numNodes;
	using CATNET<t_prob,t_prob>::m_nodeNames;
	using CATNET<t_prob,t_prob>::m_maxParents;
	using CATNET<t_prob,t_prob>::m_numParents;
	using CATNET<t_prob,t_prob>::m_parents;
	using CATNET<t_prob,t_prob>::m_maxCategories;
	using CATNET<t_prob,t_prob>::m_numCategories;
	using CATNET<t_prob,t_prob>::m_catIndices;
	using CATNET<t_prob,t_prob>::m_complexity;
	using CATNET<t_prob,t_prob>::m_loglik;
	using CATNET<t_prob,t_prob>::m_pProbLists;

protected:
	void _unlist_parent_probs(int nnode, int *pnodepars, int nodepars, t_prob *psample, 
				int &listind, t_prob *plist = 0, int parid = -1, int parcat = -1, 
				t_prob ptemp = 1, t_prob fact = 0) {
		int i, ncat, sind;
		if(!plist || listind < 0)
			return;
		/* if node probabilities are normalized, call with fact=1 */
		if(fact == 0) {
			int j;
			t_prob aux;
			fact = 1;
			sind = 0;
			for(i = 0; i < nnode; i++)
				sind += m_numCategories[i];
			aux = 0;
			for(ncat = 0; ncat < m_numCategories[nnode]; ncat++)
				aux += psample[sind + ncat];
			fact *= aux;
			if(nodepars > 0) { 
				for(j = 0; j < nodepars; j++) {
					sind = 0;
					for(i = 0; i < pnodepars[j]; i++)
						sind += m_numCategories[i];
					aux = 0;
					for(ncat = 0; ncat < m_numCategories[pnodepars[j]]; ncat++)
						aux += psample[sind + ncat];
					fact *= aux;
				}
			}
			if(fact > 0)
				fact = 1/fact;
		}
		if(nodepars == 0 || (parid == nodepars-1 && parcat >= 0 && parcat < m_numCategories[pnodepars[parid]])) {
			sind = 0;
			for(i = 0; i < nnode; i++)
				sind += m_numCategories[i];
			for(ncat = 0; ncat < m_numCategories[nnode]; ncat++)
				plist[listind++] = ptemp * psample[sind + ncat] * fact;
			return;
		}

		parid++;

		sind = 0;
		for(i = 0; i < pnodepars[parid]; i++)
			sind += m_numCategories[i];

		for(parcat = 0; parcat < m_numCategories[pnodepars[parid]]; parcat++)
			_unlist_parent_probs(nnode, pnodepars, nodepars, psample, listind, plist, parid, parcat, ptemp*psample[sind+parcat], fact);
	};

public:
	CATNETP(): CATNET<t_prob,t_prob>() {
	};

	CATNETP(int nnodes, int maxpars, int maxcats = 2, const char **nodes = 0,
			const int * pnumpars = 0, const int **ppars = 0, const int *pcats = 0): CATNET<t_prob,t_prob>(nnodes, maxpars, maxcats, nodes, pnumpars, ppars, pcats) {
	};

	t_prob nodeSampleLoglik(int nnode, int *pnodepars, int nodepars, t_prob *psamples, int nsamples, int klmode = 0) {
		int j, i, ncats, nline, listind;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, floglik, loglik, *plist;

		if(!psamples || nsamples < 1)
			return 0;

		pnodeprob = m_pProbLists[nnode]->pProbs;
		ncats = m_numCategories[nnode];
		for(i = 0; i < nodepars; i++) 
			ncats *= m_numCategories[pnodepars[i]];
		plist = (t_prob*) CATNET_MALLOC(ncats * sizeof(t_prob));
		if (!plist)
			return 0;
		
		nline = 0;
		for(j = 0; j < m_numNodes; j++) 
			nline += m_numCategories[j];

		loglik = 0;
		for (j = 0; j < nsamples; j++) {

			listind = 0;
			_unlist_parent_probs(nnode, pnodepars, nodepars, psamples+nline*j, listind, plist, -1, -1, 1, 1);
			floglik = 0;

			if(klmode) {

			for(i = 0; i < ncats; i++) {
				if(pnodeprob[i] > 0 && plist[i] > 0)
					/* use the KL distance instead */
					floglik += pnodeprob[i]*(t_prob)log(plist[i]/pnodeprob[i]);
				else if(pnodeprob[i] > 0) {
					floglik = (t_prob)-FLT_MAX;
					break;
				}
			}
			loglik += floglik;

			}			
			else { 

			for(i = 0; i < ncats; i++) {
				/* soft mode log-likelihood */
				floglik += plist[i]*pnodeprob[i];
			}
			if(floglik > 0) 
				loglik += (t_prob)log((double)floglik);
			else
				loglik = (t_prob)-FLT_MAX;

			}// klmode
		}
		if(nsamples > 1 && loglik > (t_prob)-FLT_MAX)
			loglik /= (t_prob)nsamples;
		CATNET_FREE(plist);
		return(loglik);
	}

	t_prob sampleNodeLoglik(int nnode, t_prob *psamples, int nsamples, int klmode = 0) {
		int i, j, ncats, nline, listind;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, floglik, loglik, *plist, faux;
		int nodepars, *pnodepars;

		if(!psamples || nsamples < 1 || nnode < 0 || nnode >= m_numNodes)
			return 0;

		if(!m_pProbLists || !m_pProbLists[nnode])
			return 0;

		pnodeprob = m_pProbLists[nnode]->pProbs;
		pnodepars = m_parents[nnode];
		nodepars = m_numParents[nnode];
	
		ncats = m_numCategories[nnode];
		for(i = 0; i < nodepars; i++) 
			ncats *= m_numCategories[pnodepars[i]];
		plist = (t_prob*) CATNET_MALLOC(ncats * sizeof(t_prob));
		if (!plist)
			return 0;

		nline = 0;
		for(j = 0; j < m_numNodes; j++) 
			nline += m_numCategories[j];

		loglik = 0;
		for (j = 0; j < nsamples; j++) {
			listind = 0;
			_unlist_parent_probs(nnode, pnodepars, nodepars, psamples+nline*j, listind, plist, -1, -1, 1, 1);

			floglik = 0;
			if(klmode == 1) {
				/* Pearson */ 
				for(i = 0; i < ncats; i++) {
					if(plist[i] > 0) {
						/* use the KL distance instead */
						faux = (plist[i] - pnodeprob[i]);
						floglik -= faux*faux / (t_prob)plist[i];
					}
					else {
						floglik = (t_prob)-FLT_MAX;
						break;
					}
				}
				loglik += floglik;
			}
			else if(klmode == 2) {
				/* KL */ 
				for(i = 0; i < ncats; i++) {
					if(pnodeprob[i] > 0 && plist[i] > 0) 
						/* use the KL distance instead */
						floglik += pnodeprob[i]*(t_prob)log(plist[i]/pnodeprob[i]);
					else if(pnodeprob[i] > 0) {
						floglik = (t_prob)-FLT_MAX;
						break;
					}
				}
				loglik += floglik;
			}
			else {
				for(i = 0; i < ncats; i++) {
					/* soft mode log-likelihood */
					floglik += plist[i]*pnodeprob[i];
				}
				if(floglik > 0) 
					loglik += (t_prob)log((double)floglik);
				else
					loglik = (t_prob)-FLT_MAX;
			} // klmode
		}
		CATNET_FREE(plist);

		if(nsamples > 1 && loglik > (t_prob)-FLT_MAX)
			loglik /= (t_prob)nsamples;

		return(loglik);
	}

	t_prob* sampleLoglikVector(t_prob *psamples, int nsamples, int *pert=0, int klmode = 0) {
		int i, j, nnode, ncats, nline, ncount, listind; 
		t_prob *plist, *ploglik, *pnodeprob;
		int nodepars, *pnodepars;

		if(!psamples || nsamples < 1)
			return 0;

		/* maximum probability list length */
		ncats = m_maxCategories;
		for(i = 0; i < m_maxParents; i++) 
			ncats *= m_maxCategories;
		plist = (t_prob*) CATNET_MALLOC(ncats * sizeof(t_prob));
		if (!plist)
			return 0;

		nline = 0;
		for(j = 0; j < m_numNodes; j++) 
			nline += m_numCategories[j];

		ploglik = (t_prob*) CATNET_MALLOC(m_numNodes * sizeof(t_prob));
		if (!ploglik) {
			CATNET_FREE(plist);
			return 0;
		}
		memset(ploglik, 0, m_numNodes * sizeof(t_prob));

		for (nnode = 0; nnode < m_numNodes; nnode++) {
			if(!m_pProbLists[nnode])
				continue;

			pnodepars = m_parents[nnode];
			nodepars = m_numParents[nnode];
			pnodeprob = m_pProbLists[nnode]->pProbs;

			/* get the node probability length */
			ncats = m_numCategories[nnode];
			for(i = 0; i < nodepars; i++) 
				ncats *= m_numCategories[pnodepars[i]];

			ncount = 0;
			for (j = 0; j < nsamples; j++) {
				// check for perturbation
				if(pert && pert[j * m_numNodes + nnode])
					continue; 
				listind = 0;
				_unlist_parent_probs(nnode, pnodepars, nodepars, psamples+nline*j, listind, plist, -1, -1, 1, 1);
				
				if(klmode) {
				
				for(i = 0; i < ncats; i++) {
					if(pnodeprob[i] > 0 && plist[i] > 0) 
						/* use the KL distance instead */
						ploglik[nnode] += pnodeprob[i]*(t_prob)log(plist[i]/pnodeprob[i]);
					else if(pnodeprob[i] > 0) {
						ploglik[nnode] = (t_prob)-FLT_MAX;
						break;
					}
				}

				}
				else { 

				t_prob floglik = 0;
				for(i = 0; i < ncats; i++) 
					floglik += plist[i]*(t_prob)pnodeprob[i];
				if(floglik > 0)
					ploglik[nnode] += (t_prob)log(floglik);
				else {
					ploglik[nnode] = (t_prob)-FLT_MAX;
					break;
				}

				}// klmode

				ncount++;
			}
			if(ncount > 1 && ploglik[nnode] > -FLT_MAX)
				ploglik[nnode] /= (t_prob)ncount;
		}
		CATNET_FREE(plist);
		return ploglik;
	}

	t_prob* bySampleLoglikVector(t_prob *psamples, int nsamples, int *pert=0, int klmode = 0) {
		int i, j, nnode, ncats, nline, listind; 
		t_prob *plist, *ploglik, *pnodeprob;
		int nodepars, *pnodepars;

		if(!psamples || nsamples < 1)
			return 0;

		ncats = m_maxCategories;
		for(i = 0; i < m_maxParents; i++) 
			ncats *= m_maxCategories;
		plist = (t_prob*) CATNET_MALLOC(ncats * sizeof(t_prob));
		if (!plist) {
			return 0;
		}

		nline = 0;
		for(j = 0; j < m_numNodes; j++) 
			nline += m_numCategories[j];

		ploglik = (t_prob*) CATNET_MALLOC(nsamples * sizeof(t_prob));
		if (!ploglik) {
			CATNET_FREE(plist);
			return 0;
		}
		memset(ploglik, 0, nsamples * sizeof(t_prob));

		for (nnode = 0; nnode < m_numNodes; nnode++) {
			if(!m_pProbLists[nnode])
				continue;
			pnodepars = m_parents[nnode];
			nodepars = m_numParents[nnode];
			pnodeprob = m_pProbLists[nnode]->pProbs;

			/* get the node probability length */
			ncats = m_numCategories[nnode];
			for(i = 0; i < nodepars; i++) 
				ncats *= m_numCategories[pnodepars[i]];

			for (j = 0; j < nsamples; j++) {
				// check for perturbation
				if(pert && pert[j * m_numNodes + nnode])
					continue; 
				listind = 0;
				_unlist_parent_probs(nnode, pnodepars, nodepars, psamples+nline*j, listind, plist, -1, -1, 1, 1);

				if(klmode) {

				for(i = 0; i < ncats; i++) {
					if(pnodeprob[i] > 0 && plist[i] > 0) 
						/* use the KL distance instead */
						ploglik[j] += pnodeprob[i]*(t_prob)log(plist[i]/pnodeprob[i]);
					else if(pnodeprob[i] > 0) {
						ploglik[j] = (t_prob)-FLT_MAX;
						break;
					}
				}

				}
				else { 

				t_prob floglik = 0;
				for(i = 0; i < ncats; i++) 
					floglik += plist[i]*(t_prob)pnodeprob[i];
				if(floglik > 0)
					ploglik[j] += (t_prob)log(floglik);
				else {
					ploglik[j] = (t_prob)-FLT_MAX;
					break;
				}

				}// klmode


			}
		}

		CATNET_FREE(plist);
		return ploglik;
	}

	// sets sample conditional probability and returns its log-likelihood
	t_prob setNodeSampleProb(int nnode, t_prob *psamples, int nsamples, int bNormalize = 0) {
		int i, j, nline, *pnodepars, nodepars, ncats, listind;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, *plist, floglik;
		
		if (!m_pProbLists || !psamples || nsamples < 1) {
			return (t_prob)-FLT_MAX;
		}

		if(!m_pProbLists[nnode]) {
			// error
			return (t_prob)-FLT_MAX;
		}
		else
			m_pProbLists[nnode]->set_zero();

		pnodeprob = m_pProbLists[nnode]->pProbs;
		pnodepars = m_parents[nnode];
		nodepars  = m_numParents[nnode];

		ncats = m_numCategories[nnode];
		for(i = 0; i < nodepars; i++) 
			ncats *= m_numCategories[pnodepars[i]];
		plist = (t_prob*) CATNET_MALLOC(ncats * sizeof(t_prob));
		if (!plist) {
			return 0;
		}

		nline = 0;
		for(j = 0; j < m_numNodes; j++) 
			nline += m_numCategories[j];
		for (j = 0; j < nsamples; j++) {
			listind = 0;
			_unlist_parent_probs(nnode, pnodepars, nodepars, psamples+nline*j, listind, plist, -1, -1, 1, 1);
			for(i = 0; i < ncats; i++) 
				pnodeprob[i] += plist[i];
		}
		CATNET_FREE(plist);

		/* keep sample sizes */
		m_pProbLists[nnode]->sampleSize = nsamples;

		/* at this point m_pProbLists[nnode] has counts not probabilities 
		  find m_pProbLists[nnode]->loglik */
		m_pProbLists[nnode]->loglikelihood();
		if(nsamples > 1) {
			m_pProbLists[nnode]->loglik /= (t_prob)nsamples;
			m_pProbLists[nnode]->priorlik /= (t_prob)nsamples;
		}
		floglik = m_pProbLists[nnode]->loglik + m_pProbLists[nnode]->priorlik;

		if(bNormalize)
			m_pProbLists[nnode] -> normalize();

		/* need to be calculated next call time */
		m_loglik = 0;

		return floglik;
	}

	// sets sample conditional probability and returns its log-likelihood
	t_prob setNodeSampleProbKL(int nnode, t_prob *psamples, int nsamples, int *pClasses, int bNormalize = 0, int bUsePearson = 0) {
		int i, j, nline, *pnodepars, nodepars, ncats, listind, ncount;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, *plist, floglik;
		PROB_LIST<t_prob> *pPriorProb = 0;

		if (!m_pProbLists || !psamples || nsamples < 1 || !pClasses) 
			return (t_prob)-FLT_MAX;

		pPriorProb = new PROB_LIST<t_prob>;
		if(!m_pProbLists[nnode] || !pPriorProb)
			return (t_prob)-FLT_MAX;
		*pPriorProb = *m_pProbLists[nnode];

		m_pProbLists[nnode]->set_zero();

		pnodepars = m_parents[nnode];
		nodepars = m_numParents[nnode];

		ncats = m_numCategories[nnode];
		for(i = 0; i < nodepars; i++) 
			ncats *= m_numCategories[pnodepars[i]];
		plist = (t_prob*) CATNET_MALLOC(ncats * sizeof(t_prob));
		if (!plist) {
			delete pPriorProb;
			return 0;
		}

		nline = 0;
		for(j = 0; j < m_numNodes; j++) 
			nline += m_numCategories[j];

		pnodeprob = m_pProbLists[nnode]->pProbs;		
		ncount = 0;
		for (j = 0; j < nsamples; j++) {
			listind = 0;
			_unlist_parent_probs(nnode, pnodepars, nodepars, psamples+nline*j, listind, plist, -1, -1, 1, 1);
			for(i = 0; i < ncats; i++) 
				pnodeprob[i] += plist[i];
			ncount++;
		}
		/* keep sample sizes */
		m_pProbLists[nnode]->sampleSize = ncount;

		pPriorProb->set_zero();
		pnodeprob = pPriorProb->pProbs;		
		for (j = 0; j < nsamples; j++) {
			if(pClasses[j] != 0)
				continue;
			listind = 0;
			_unlist_parent_probs(nnode, pnodepars, nodepars, psamples+nline*j, listind, plist, -1, -1, 1, 1);
			for(i = 0; i < ncats; i++) 
				pnodeprob[i] += plist[i];
		}
		/* keep sample sizes */
		pPriorProb->sampleSize = nsamples;

		CATNET_FREE(plist);

		//pPriorProb->normalize();
		m_pProbLists[nnode]->setPrior(pPriorProb, bUsePearson);

		/* at this point m_pProbLists[nnode] has counts not probabilities 
		  find m_pProbLists[nnode]->loglik */
		m_pProbLists[nnode]->loglikelihood();
		if(ncount > 1) {
			m_pProbLists[nnode]->loglik /= (t_prob)ncount;
			m_pProbLists[nnode]->priorlik /= (t_prob)ncount;
		}
		floglik = m_pProbLists[nnode]->loglik + m_pProbLists[nnode]->priorlik;

		if(bNormalize)
			m_pProbLists[nnode] -> normalize();

		/* need to be calculated next call time */
		m_loglik = 0;

		delete pPriorProb;
		return floglik;
	}

};

#endif /* CATNETP_H_ */
