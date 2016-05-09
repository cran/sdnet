/*
 *  catnetd : categorical Bayesian network inference
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
 * catnetd.h
 *
 *  Created on: Sep 18, 2009
 *      Author: nbalov
 */

#include "catnet_base.h"

#ifndef CATNETD_H_
#define CATNETD_H_

template<class t_prob>
class CATNETD: public CATNET<int, t_prob> {

protected:
	using CATNET<int,t_prob>::m_numNodes;
	using CATNET<int,t_prob>::m_nodeNames;
	using CATNET<int,t_prob>::m_maxParents;
	using CATNET<int,t_prob>::m_numParents;
	using CATNET<int,t_prob>::m_parents;
	using CATNET<int,t_prob>::m_maxCategories;
	using CATNET<int,t_prob>::m_numCategories;
	using CATNET<int,t_prob>::m_catIndices;
	using CATNET<int,t_prob>::m_complexity;
	using CATNET<int,t_prob>::m_loglik;
	using CATNET<int,t_prob>::m_pProbLists;

public:
	CATNETD(): CATNET<int,t_prob>() {
	};

	CATNETD(int nnodes, int maxpars, int maxcats = 2, const char **nodes = 0,
			const int * pnumpars = 0, const int **ppars = 0, const int *pcats = 0): CATNET<int,t_prob>(nnodes, maxpars, maxcats, nodes, pnumpars, ppars, pcats) {
	};

	t_prob nodeSampleLoglik(int nnode, int *pnodepars, int nodepars, int *psamples, int nsamples, int klmode = 0) {
		int j, ipar;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, floglik = 0;
		int *pnodesample, samp, ncount;

		if(!psamples || nsamples < 1)
			return 0;

		pnodesample = 0;
		if (m_maxParents > 0) {
			pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		}

		pnodepars = m_parents[nnode];
		ncount = 0;
		for (j = 0; j < nsamples; j++) {
			if (pnodesample) {
				for (ipar = 0; ipar < nodepars; ipar++) {
					pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
				}
			}
			pnodeprob = m_pProbLists[nnode]->find_slot(0, pnodesample, 0);
			if(!pnodeprob)
				continue;
			samp = psamples[j * m_numNodes + nnode];
			if (samp >= 0 && samp < m_numCategories[nnode]) {
				if(pnodeprob[samp] > 0)
					floglik += (t_prob)log((double)pnodeprob[samp]);
				else {
					floglik = (t_prob)-FLT_MAX;
					break;
				}
				ncount++;
			}
		}
		if(ncount > 1 && floglik > (t_prob)-FLT_MAX)
			floglik /= (t_prob)ncount;
		CATNET_FREE(pnodesample);
		return(floglik);
	}

	t_prob sampleLoglik(int *psamples, int nsamples) {
		int i, j, ipar;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, nodeloglik;
		int nodepars, ncount;
		int *pnodepars, *pnodesample, samp;

		if(!psamples || nsamples < 1)
			return 0;

		pnodesample = 0;
		if (m_maxParents > 0) {
			pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		}

		m_loglik = 0;
		for (i = 0; i < m_numNodes; i++) {
			if(!m_pProbLists[i])
				continue;
			pnodepars = m_parents[i];
			nodepars = m_numParents[i];
			ncount = 0;
			nodeloglik = 0;
			for (j = 0; j < nsamples; j++) {
				if (pnodesample) {
					for (ipar = 0; ipar < nodepars; ipar++) {
						pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
					}
				}
				pnodeprob = m_pProbLists[i]->find_slot(0, pnodesample, 0);
				if(!pnodeprob)
					continue;
				samp = psamples[j * m_numNodes + i];
				if (pnodeprob && samp >= 0 && samp < m_numCategories[i]) {
					if(pnodeprob[samp] > 0)
						nodeloglik += (double)log((double)pnodeprob[samp]);
					else {
						nodeloglik = (double)-FLT_MAX;
						break;
					}
					ncount++;
				}
			}
			if(ncount > 1 && nodeloglik > (t_prob)-FLT_MAX)
				nodeloglik /= (t_prob)ncount;
			if(nodeloglik == -FLT_MAX) {
				m_loglik = -FLT_MAX;
				break;
			}
			else
				m_loglik += nodeloglik;
		}
		CATNET_FREE(pnodesample);
		return m_loglik;
	}

	t_prob* sampleLoglikVector(int *psamples, int nsamples, int *pert=0, int klmode = 0) {
		int i, j, ipar;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, *ploglik;
		int nodepars, ncount, *pnodepars, *pnodesample, samp;

		if(!psamples || nsamples < 1)
			return 0;

		pnodesample = 0;
		if (m_maxParents > 0) {
			pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		}

		ploglik = (t_prob*) CATNET_MALLOC(m_numNodes * sizeof(t_prob));
		if (!ploglik) {
			CATNET_FREE(pnodesample);
			return 0;
		}
		memset(ploglik, 0, m_numNodes * sizeof(t_prob));

		for (i = 0; i < m_numNodes; i++) {
			if(!m_pProbLists[i])
				continue;
			pnodepars = m_parents[i];
			nodepars = m_numParents[i];

			ncount = 0;
			for (j = 0; j < nsamples; j++) {
				// check for perturbation
				if(pert && pert[j * m_numNodes + i])
					continue;
				if (pnodesample) {
					for (ipar = 0; ipar < nodepars; ipar++) {
						pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
					}
				}
				pnodeprob = m_pProbLists[i]->find_slot(0, pnodesample, 0);
				if(!pnodeprob)
					continue;
				samp = psamples[j * m_numNodes + i];
				if (pnodeprob && samp >= 0 && samp < m_numCategories[i]) {
					if(pnodeprob[samp] > 0)
						ploglik[i] += (t_prob)log((double)pnodeprob[samp]);
					else {
						ploglik[i] = (t_prob)-FLT_MAX;
						break;
					}
					ncount++;
				}
			}
			if(ncount > 1 && ploglik[i] > -FLT_MAX)
				ploglik[i] /= (t_prob)ncount;
		}
		CATNET_FREE(pnodesample);
		return ploglik;
	}

	t_prob* bySampleLoglikVector(int *psamples, int nsamples, int *pert=0, int klmode = 0) {
		int i, j, ipar;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, *ploglik;
		int nodepars, *pnodepars, *pnodesample, samp;

		if(!psamples || nsamples < 1)
			return 0;

		pnodesample = 0;
		if (m_maxParents > 0) {
			pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		}

		ploglik = (t_prob*) CATNET_MALLOC(nsamples * sizeof(t_prob));
		if (!ploglik) {
			CATNET_FREE(pnodesample);
			return 0;
		}
		memset(ploglik, 0, nsamples * sizeof(t_prob));

		for (i = 0; i < m_numNodes; i++) {
			if(!m_pProbLists[i])
				continue;
			pnodepars = m_parents[i];
			nodepars = m_numParents[i];

			for (j = 0; j < nsamples; j++) {
				// check for perturbation
				if(pert && pert[j * m_numNodes + i])
					continue; 
				if (pnodesample) {
					for (ipar = 0; ipar < nodepars; ipar++) {
						pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
					}
				}
				pnodeprob = m_pProbLists[i]->find_slot(0, pnodesample, 0);
				if(!pnodeprob)
					continue;
				samp = psamples[j * m_numNodes + i];
				if (pnodeprob && samp >= 0 && samp < m_numCategories[i]) {
					if(pnodeprob[samp] > 0)
						ploglik[j] += (t_prob)log((double)pnodeprob[samp]);
					else {
						ploglik[j] = (t_prob)-FLT_MAX;
						break;
					}
				}
			}
		}

		CATNET_FREE(pnodesample);
		return ploglik;
	}

	t_prob sampleNodeLoglik(int nnode, int *psamples, int nsamples, int klmode = 0) {
		int j, ipar;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, loglik;
		int nodepars, ncount;
		int *pnodepars, *pnodesample, samp;

		if(!psamples || nsamples < 1 || nnode < 0 || nnode >= m_numNodes)
			return 0;

		if(!m_pProbLists || !m_pProbLists[nnode])
			return 0;

		loglik = 0;

		pnodesample = 0;
		if (m_maxParents > 0) {
			pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		}

		pnodepars = m_parents[nnode];
		nodepars = m_numParents[nnode];

		ncount = 0;	
		for (j = 0; j < nsamples; j++) {
			if (pnodesample) {
				for (ipar = 0; ipar < nodepars; ipar++) 
					pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
			}
			pnodeprob = m_pProbLists[nnode]->find_slot(0, pnodesample, 0);

			samp = psamples[j * m_numNodes + nnode];
			if(klmode == 1) {
				/* Pearson */ 
			}
			if(klmode == 2) {
				/* KL */ 
			}
			else {
				// if there are NAs, pnodeprob maybe 0
				if (pnodeprob && samp >= 0 && samp < m_numCategories[nnode]) {
					if(pnodeprob[samp] > 0)
						loglik += (t_prob)log((double)pnodeprob[samp]);
					else {
						loglik = (t_prob)-FLT_MAX;
						break;
					}
					ncount++;
				}
			} // klmode
		}
		CATNET_FREE(pnodesample);

		if(ncount > 1 && loglik > (t_prob)-FLT_MAX)
			loglik /= (t_prob)ncount;

		return(loglik);
	}

	// sets sample conditional probability and returns its log-likelihood
	t_prob setNodeSampleProb(int nnode, int *psamples, int nsamples, int bNormalize = 0) {
		int i, j;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, floglik;
		int *pnodesample, *pnodepars, samp, ncount;

		if (!m_pProbLists || !psamples || nsamples < 1) {
			return (t_prob)-FLT_MAX;
		}

		if(!m_pProbLists[nnode]) {
			return (t_prob)-FLT_MAX;
		}
		m_pProbLists[nnode]->set_zero();

		pnodesample = 0;
		if (m_maxParents > 0) {
			pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		}

		pnodepars = m_parents[nnode];
		ncount = 0;
		for (j = 0; j < nsamples; j++) {
			if (pnodesample) {
				for (i = 0; i < m_numParents[nnode]; i++) {
					if (pnodepars[i] < 0 || pnodepars[i] >= m_numNodes)
						CATNET_ERR("pnodepars[i]");
					pnodesample[i] = psamples[j * m_numNodes + pnodepars[i]];
				}
			}
			pnodeprob = m_pProbLists[nnode]->find_slot(0, pnodesample, 0);
			// if there are NAs, pnodeprob maybe 0
			samp = psamples[j * m_numNodes + nnode];
			if (pnodeprob && samp >= 0 && samp < m_numCategories[nnode]) {
				pnodeprob[samp]++;
				ncount++;
			}
		}

		CATNET_FREE(pnodesample);

		/* keep sample sizes */
		m_pProbLists[nnode]->sampleSize = ncount;

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

		return(floglik);
	}

	// sets sample conditional probability and returns its log-likelihood
	t_prob setNodeSampleProbKL(int nnode, int *psamples, int nsamples, int *pClasses, int bNormalize = 0, int bUsePearson = 0) {
		int i, j;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, floglik;
		int *pnodesample, *pnodepars, samp, ncount, nclasscount;
		PROB_LIST<t_prob> *pClassProb = 0;

		if (!m_pProbLists || !pClasses || !psamples || nsamples < 1)
			return (t_prob)-FLT_MAX;
		
		pClassProb = new PROB_LIST<t_prob>;
		if(!m_pProbLists[nnode] || !pClassProb)
			return (t_prob)-FLT_MAX;
		*pClassProb = *m_pProbLists[nnode];

		m_pProbLists[nnode]->set_zero();
		pClassProb->set_zero();

		pnodesample = 0;
		if (m_maxParents > 0) {
			pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		}

		pnodepars = m_parents[nnode];
		ncount = 0;
		nclasscount = 0;
		for (j = 0; j < nsamples; j++) {
			if (pnodesample) {
				for (i = 0; i < m_numParents[nnode]; i++) {
					if (pnodepars[i] < 0 || pnodepars[i] >= m_numNodes)
						CATNET_ERR("pnodepars[i]");
					pnodesample[i] = psamples[j * m_numNodes + pnodepars[i]];
				}
			}
			pnodeprob = m_pProbLists[nnode]->find_slot(0, pnodesample, 0);
			// if there are NAs, pnodeprob maybe 0
			samp = psamples[j * m_numNodes + nnode];
			if (pnodeprob && samp >= 0 && samp < m_numCategories[nnode]) {
				pnodeprob[samp]++;
				ncount++;
			}
			if(pClasses[j] != 0)
				continue;
			pnodeprob = pClassProb->find_slot(0, pnodesample, 0);
			// if there are NAs, pnodeprob maybe 0
			samp = psamples[j * m_numNodes + nnode];
			if (pnodeprob && samp >= 0 && samp < m_numCategories[nnode]) {
				pnodeprob[samp]++;
				nclasscount++;
			}
		}
		CATNET_FREE(pnodesample);

		/* keep sample sizes */
		m_pProbLists[nnode]->sampleSize = ncount;
		pClassProb->sampleSize = nclasscount;

		//pClassProb->normalize();
		m_pProbLists[nnode]->setPrior(pClassProb, bUsePearson);

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

		delete pClassProb;
			
		return(floglik);
	}

};

#endif /* CATNETD_H_ */
