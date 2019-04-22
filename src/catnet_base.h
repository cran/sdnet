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
 * catnet_base.h
 *
 *  Created on: Sep 18, 2009
 *      Author: nbalov
 */

#include "utils.h"
#include "problist.h"

#ifndef CATNET_BASE_H_
#define CATNET_BASE_H_

//#define DEBUG_ON	1 

template<class t_sample, class t_prob>
class CATNET {
protected:
	/* nodes are assumed ordered */
	int m_numNodes;
	char **m_nodeNames;
	int m_maxParents;
	int *m_numParents;
	int **m_parents;
	int m_maxCategories;
	int *m_numCategories;
	int **m_catIndices;
	int m_complexity;
	t_prob m_loglik;
	PROB_LIST<t_prob> **m_pProbLists;

public:
	CATNET() {
		_reset();
	}

	CATNET(int nnodes, int maxpars, int maxcats = 2, const char **nodes = 0,
			const int * pnumpars = 0, const int **ppars = 0, const int *pcats = 0) {
		_reset();
		init(nnodes, maxpars, maxcats, nodes, pnumpars, ppars, pcats);
	}

	virtual ~CATNET() {
		_release();
	}

	CATNET<t_sample, t_prob>& operator =(const CATNET<t_sample, t_prob> &cnet) {
		init(cnet.m_numNodes, cnet.m_maxParents, cnet.m_maxCategories,
				(const char**)cnet.m_nodeNames, (const int*)cnet.m_numParents, (const int**)cnet.m_parents,
				(const int*)cnet.m_numCategories, (const PROB_LIST<t_prob> **)cnet.m_pProbLists);
		m_complexity = cnet.m_complexity;
		m_loglik = cnet.m_loglik;
		return *this;
	}

	virtual void _release() {
		int i;
		for (i = 0; i < m_numNodes; i++) {
			if (m_pProbLists && m_pProbLists[i]) {
				delete m_pProbLists[i];
				m_pProbLists[i] = 0;
			}
			if (m_parents && m_parents[i]) {
				CATNET_FREE(m_parents[i]);
				m_parents[i] = 0;
			}
			if (m_nodeNames && m_nodeNames[i]) {
				CATNET_FREE(m_nodeNames[i]);
				m_nodeNames[i] = 0;
			}
			if(m_catIndices && m_catIndices[i]) {
				CATNET_FREE(m_catIndices[i]);
				m_catIndices[i] = 0;
			}
		}
		if (m_numParents)
			CATNET_FREE(m_numParents);
		if (m_parents)
			CATNET_FREE(m_parents);
		if (m_numCategories)
			CATNET_FREE(m_numCategories);
		if (m_nodeNames)
			CATNET_FREE(m_nodeNames);
		if(m_catIndices)
			CATNET_FREE(m_catIndices);
		if (m_pProbLists)
			CATNET_FREE(m_pProbLists);

		_reset();
	}

	virtual void _reset() {
		m_numNodes = 0;
		m_maxParents = 0;
		m_maxCategories = 0;
		m_nodeNames = 0;
		m_numParents = 0;
		m_parents = 0;
		m_numCategories = 0;
		m_catIndices = 0;
		m_pProbLists = 0;
		m_complexity = 0;
		m_loglik = 0;
	}

	void init(int nnodes, int maxpars, int maxcats = 2, const char **nodes = 0,
			const int *pnumpars = 0, const int **ppars = 0, const int *pcats = 0, 
			const PROB_LIST<t_prob> **pprobs = 0) {

		_release();

		int i, j, *nodeparcats, nodenamelen;

		if (nnodes < 1 || maxpars < 0)
			return;
		m_numNodes = nnodes;
		m_maxParents = maxpars;
		m_maxCategories = maxcats;

		m_numParents    = (int*)   CATNET_MALLOC(m_numNodes * sizeof(int));
		
		m_parents       = (int**)  CATNET_MALLOC(m_numNodes * sizeof(int*));
		
		m_numCategories = (int*)   CATNET_MALLOC(m_numNodes * sizeof(int));

		m_catIndices    = (int**)  CATNET_MALLOC(m_numNodes * sizeof(int*));		

 		m_nodeNames     = (char**) CATNET_MALLOC(m_numNodes * sizeof(char*));

		m_pProbLists    = (PROB_LIST<t_prob>**) CATNET_MALLOC(m_numNodes * sizeof(PROB_LIST<t_prob>*));

		if (!m_numParents || !m_parents || !m_numCategories || !m_catIndices 
                 || !m_nodeNames || !m_pProbLists) {
			if (m_numParents)
				CATNET_FREE(m_numParents);
			if (m_parents)
				CATNET_FREE(m_parents);
			if (m_numCategories)
				CATNET_FREE(m_numCategories);
			if (m_nodeNames)
				CATNET_FREE(m_nodeNames);
			if(m_catIndices)
				CATNET_FREE(m_catIndices);
			if (m_pProbLists)
				CATNET_FREE(m_pProbLists);
			m_numNodes      = 0;
			m_maxParents    = 0;
			m_maxCategories = 0;
			m_nodeNames     = 0;
			m_numParents    = 0;
			m_parents       = 0;
			m_numCategories = 0;
			m_catIndices    = 0;
			m_pProbLists    = 0;
			return;
		}

		memset(m_parents, 0, m_numNodes * sizeof(int*));

		if (pcats != 0)
			memcpy(m_numCategories, pcats, m_numNodes * sizeof(int));
		else {
			for (i = 0; i < m_numNodes; i++)
				m_numCategories[i] = m_maxCategories;
		}

		memset(m_catIndices, 0, m_numNodes * sizeof(int*));

		if (nodes) {	
			for (i = 0; i < m_numNodes; i++) {
				m_nodeNames[i] = 0;
				if (!nodes[i])
					continue;
				nodenamelen = strlen(nodes[i]);
				m_nodeNames[i] = (char*) CATNET_MALLOC((nodenamelen+1) * sizeof(char));
				if (m_nodeNames[i])
					strcpy(m_nodeNames[i], nodes[i]);
			}
		}
		else {
			memset(m_nodeNames, 0, m_numNodes * sizeof(char*));
		}

		memset(m_numParents, 0, m_numNodes * sizeof(int));
		if (pnumpars != 0 && ppars != 0) {
			memcpy(m_numParents, pnumpars, m_numNodes * sizeof(int));
			for (i = 0; i < m_numNodes; i++) {
				if (m_numParents[i] <= 0) 
					continue;
				m_parents[i] = (int*) CATNET_MALLOC(m_numParents[i] * sizeof(int));
				if (m_parents[i]) {
					memset(m_parents[i], 0, m_numParents[i] * sizeof(int));
					if(ppars[i])
						memcpy(m_parents[i], ppars[i], m_numParents[i] * sizeof(int));
				}
				if (m_numParents[i] > m_maxParents)
					m_maxParents = m_numParents[i];
			}
		}
		
		memset(m_pProbLists, 0, m_numNodes * sizeof(PROB_LIST<t_prob>*));
		nodeparcats = (int*)CATNET_MALLOC((m_maxParents+1)*sizeof(int));
		if (nodeparcats) {	
			for (i = 0; i < m_numNodes; i++) {
				if (pprobs && pprobs[i]) {
					m_pProbLists[i] = new PROB_LIST<t_prob>;
					*m_pProbLists[i] = *pprobs[i];
				}
				else {
					for(j = 0; j < m_numParents[i]; j++)
						nodeparcats[j] = m_numCategories[m_parents[i][j]]; 
					m_pProbLists[i] = new PROB_LIST<t_prob>(m_numCategories[i], m_maxCategories, m_numParents[i], nodeparcats);
				}
			}
			CATNET_FREE(nodeparcats);
		}
	}

	void setNodeNames(char **pnames, const int *porder) {
		int i;
		const char *str;
		if (!porder || !pnames) 
			return;
		if(!m_nodeNames) {
			m_nodeNames = (char**) CATNET_MALLOC(m_numNodes * sizeof(char*));
			if (m_nodeNames)
				memset(m_nodeNames, 0, m_numNodes * sizeof(char*));
		}
		for (i = 0; i < m_numNodes; i++) {
			m_nodeNames[i] = 0;
			if (porder[i] < 1 || porder[i] > m_numNodes)
				break;
			str = pnames[porder[i]-1];
			if(m_nodeNames[i])
				CATNET_FREE(m_nodeNames[i]);
			m_nodeNames[i] = (char*) CATNET_MALLOC((strlen(str)+1) * sizeof(char));
			if (m_nodeNames[i] && str)
				strcpy((char*)m_nodeNames[i], str);
		}
	}

	void setNodesOrder(const int *porder) {
		int i;
		char str[256];
		if (!porder) 
			return;
		if(!m_nodeNames) {
			m_nodeNames = (char**) CATNET_MALLOC(m_numNodes * sizeof(char*));
			memset(m_nodeNames, 0, m_numNodes * sizeof(char*));
		}
		for (i = 0; i < m_numNodes; i++) {
			m_nodeNames[i] = 0;
			if (porder[i] < 1 || porder[i] > m_numNodes)
				break;
			sprintf(str, "N%d", (int)porder[i]);
			if(m_nodeNames[i])
				CATNET_FREE(m_nodeNames[i]);
			m_nodeNames[i] = (char*) CATNET_MALLOC((strlen(str)+1) * sizeof(char));
			if (m_nodeNames[i])
				strcpy((char*)m_nodeNames[i], str);
		}
	}

	void setCategoryIndices(int *pNumCats, int **pCatIndices) {
		int i;
		// reset node names only
		if (!pNumCats || !pCatIndices || !m_numCategories) 
			return;
		if(!m_catIndices) {
			m_catIndices = (int**) CATNET_MALLOC(m_numNodes * sizeof(int*));
			memset(m_catIndices, 0, m_numNodes * sizeof(int*));
		}
		if(!m_catIndices)
			return;
		for (i = 0; i < m_numNodes; i++) {
			if(pNumCats[i] != m_numCategories[i])
				break;
			if(!m_catIndices[i])
				m_catIndices[i] = (int*) CATNET_MALLOC(m_numCategories[i] * sizeof(int));
			memcpy(m_catIndices[i], pCatIndices[i], m_numCategories[i] * sizeof(int));
		}
	}

	void setCondProb(int node, t_prob *pcondprob, int condprobsize) {
		if (m_numNodes < 1 || node < 0 || node >= m_numNodes)
			return;
		if (m_pProbLists && m_pProbLists[node])
			delete m_pProbLists[node];
		if (m_numParents[node] > m_maxParents)
			m_maxParents = m_numParents[node];
		if (!m_pProbLists) {
			m_pProbLists = (PROB_LIST<t_prob>**) CATNET_MALLOC(m_numNodes
					* sizeof(PROB_LIST<t_prob>*));
			memset(m_pProbLists, 0, m_numNodes * sizeof(PROB_LIST<t_prob>*));
		}
		if(!m_pProbLists)
			return;
		m_pProbLists[node] = 0;
		if (m_numParents[node] < 0 || m_numParents[node] > m_maxParents)
			return;
		int *parcats = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		for (int i = 0; i < m_numParents[node]; i++) {
			parcats[i] = m_numCategories[m_parents[node][i]];
		}
		m_pProbLists[node] = new PROB_LIST<t_prob> (m_numCategories[node],
				m_maxCategories, m_numParents[node], parcats, pcondprob,
				condprobsize);
		CATNET_FREE(parcats);
	}

	int numNodes() {
		return m_numNodes;
	}

	const char** nodeNames() {
		return (const char**)m_nodeNames;
	}

	int maxCategories() {
		return m_maxCategories;
	}

	const int* numCategories() {
		return (const int*)m_numCategories;
	}

	int numCategories(int node) {
		if (node < 0 || node >= m_numNodes)
			return 0;
		return m_numCategories[node];
	}

	int maxParents() {
		return m_maxParents;
	}

	const int** parents() {
		return (const int**)m_parents;
	}

	const int* numParents() {
		return (const int*)m_numParents;
	}

	int numParents(int node) {
		if (node < 0 || node >= m_numNodes)
			return 0;
		return m_numParents[node];
	}

	const int* getParents(int node) {
		if (node < 0 || node >= m_numNodes)
			return 0;
		return (const int*)m_parents[node];
	}

	int setParents(int node, int* parents, int numparents) {
		if (node < 0 || node >= m_numNodes)
			return 0;
		if (!m_numParents || !m_parents)
			return 0;
		if(m_numParents[node] != numparents) {
			m_numParents[node] = numparents;
			if(m_parents[node])
				CATNET_FREE(m_parents[node]);
			m_parents[node] = (int*) CATNET_MALLOC(m_numParents[node] * sizeof(int));
		}
		if (!m_parents[node])
			return 0;
		memcpy(m_parents[node], parents, m_numParents[node] * sizeof(int));

		if (!m_pProbLists)
			m_pProbLists = (PROB_LIST<t_prob>**) CATNET_MALLOC(m_numNodes * sizeof(PROB_LIST<t_prob>*));
		else {
			if(m_pProbLists[node])
				delete m_pProbLists[node];
			m_pProbLists[node] = 0;
		}

		int *parcats = (int*)CATNET_MALLOC(m_numParents[node]*sizeof(int));
		if (!parcats)
			return 0;
		for(int i = 0; i < m_numParents[node]; i++) 
			parcats[i] = m_numCategories[parents[i]];
		m_pProbLists[node] = new PROB_LIST<t_prob>(m_numCategories[node], m_maxCategories, m_numParents[node], parcats);
		CATNET_FREE(parcats);
		             
		if(m_maxParents < numparents)
			m_maxParents = numparents;

		/* need to be calculated next time */
		m_complexity = 0;
		m_loglik = 0;

		return m_numParents[node];
	}

	int complexity() {
		if(m_complexity < m_numNodes)
			return findComplexity();
		return m_complexity;
	}

	int findComplexity() {
		int i, j, c;
		m_complexity = 0;
		for (i = 0; i < m_numNodes; i++) {
			if (!m_parents || !m_parents[i]) {
				m_complexity += (m_numCategories[i]-1);
				continue;
			}
			c = 1;
			for (j = 0; j < m_numParents[i]; j++)
				c *= m_numCategories[m_parents[i][j]];
			m_complexity += c*(m_numCategories[i]-1);
		}
		return m_complexity;
	}

	int nodeComplexity(int nnode) {
		int j, c;
		if(nnode < 0 || nnode >= m_numNodes)
			return(0);
		if (!m_parents || !m_parents[nnode])
			return(m_numCategories[nnode]-1);
		c = (m_numCategories[nnode]-1);
		for (j = 0; j < m_numParents[nnode]; j++)
			c *= m_numCategories[m_parents[nnode][j]];
		return c;
	}

	t_prob loglik() {
		if(m_loglik == 0)
			return findLoglik();
		return m_loglik;
	}

	t_prob findLoglik() {
		int i;
		if(!m_pProbLists)
			return -FLT_MAX;
		m_loglik = 0;
		for (i = 0; i < m_numNodes; i++) {
			if(m_pProbLists[i])
				m_loglik += (m_pProbLists[i]->loglik + m_pProbLists[i]->priorlik);
		}
		return m_loglik;
	}
	
	const PROB_LIST<t_prob>** probLists() {
		return (const PROB_LIST<t_prob>**)m_pProbLists;
	}

	PROB_LIST<t_prob>* getNodeProb(int nnode) {
		if (!m_pProbLists || nnode < 0 || nnode >= m_numNodes)
			return 0;
		return m_pProbLists[nnode];
	}

	const PROB_LIST<t_prob>* setNodeProb(int nnode, const PROB_LIST<t_prob> *pprob) {
		if(!m_pProbLists || nnode < 0 || nnode >= m_numNodes)
			return 0;
		if(!m_pProbLists[nnode])
			m_pProbLists[nnode] = new PROB_LIST<t_prob>;
		if (m_pProbLists[nnode] && pprob) {
			*m_pProbLists[nnode] = *pprob;
		}
		return m_pProbLists[nnode];
	}

	void setNodePriorProb(int nnode, t_prob pprior) {
		if(!m_pProbLists || nnode < 0 || nnode >= m_numNodes || !m_pProbLists[nnode])
			return;
		m_pProbLists[nnode]->priorlik = pprior;
		m_loglik = 0;
	}

	int normalizeProbabilities() {
		int i;
		if(!m_pProbLists)
			return -1;
		for (i = 0; i < m_numNodes; i++) {
			if(m_pProbLists[i]) {
				m_pProbLists[i] -> normalize();
			}
		}
		return 0;
	}

	int *findParentPool(int nnode, int &poolsize) {

		int ipar, par, i, j, bfound, parpoolsize, *parpool, *ppool, *paux;

		poolsize = 0;
		if (nnode < 0 || nnode >= m_numNodes || !m_parents || !m_parents[nnode] || !m_numParents[nnode])
			return 0;

		ppool = 0;
		for (ipar = 0; ipar < m_numParents[nnode]; ipar++) {
			par = m_parents[nnode][ipar];
			parpool = findParentPool(par, parpoolsize);

			i = 0;
			while (i < parpoolsize) {
				bfound = 0;
				for (j = 0; j < poolsize; j++) {
					if (parpool[i] == ppool[j]) {
						bfound = 1;
						break;
					}
				}
				if (!bfound) {
					i++;
					continue;
				}
				for (j = i + 1; j < parpoolsize; j++)
					parpool[j - 1] = parpool[j];
				parpoolsize--;
			}

			paux = (int*) CATNET_MALLOC((poolsize + parpoolsize + 1) * sizeof(int));
			if (!paux) 
				break;

			if (ppool && poolsize > 0)
				memcpy(paux, ppool, poolsize * sizeof(int));
			if (parpoolsize > 0 && parpool)
				memcpy(paux + poolsize, parpool, parpoolsize * sizeof(int));
			// release that
			if(parpool)
				CATNET_FREE(parpool);

			// check par
			bfound = 0;
			for (j = 0; j < poolsize + parpoolsize; j++) {
				if (paux[j] == par) {
					bfound = 1;
					break;
				}
			}
			if (bfound) {
				poolsize += parpoolsize;
			} else {
				paux[poolsize + parpoolsize] = par;
				poolsize += (parpoolsize + 1);
			}

			if (ppool)
				CATNET_FREE(ppool);
			ppool = paux;
		}
		return ppool;
	}

	t_prob *findJointProb(int nnode, int &jointprobsize) {
		t_prob *jointprob, *parprob;
		int i, ii, i0, ic, ipar, j, par, ipool, poolnode, cc;
		int *parpool, *paux, parpoolsize, *blocksize, *paridx, *pcats;
		PROB_LIST<t_prob> *probnode = 0;

		if (nnode < 0 || nnode >= m_numNodes || !m_pProbLists)
			return 0;

		parpool = findParentPool(nnode, parpoolsize);
		// add nnode to the parents list
		paux = (int*) CATNET_MALLOC((parpoolsize + 1) * sizeof(int));
		if (!paux)
			return 0;
		if (parpool && parpoolsize > 0) 
			memcpy(paux, parpool, parpoolsize * sizeof(int));
		if (parpool)
			CATNET_FREE(parpool);
		paux[parpoolsize] = nnode;
		parpoolsize++;
		parpool = paux;

		pcats = 0;
		if (m_maxParents > 0) {
			pcats     = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
			if (!pcats)
				return 0;
		}
		paridx    = 0;
		blocksize = 0;
		if (parpoolsize > 0) {
			paridx    = (int*) CATNET_MALLOC(parpoolsize  * sizeof(int));
			blocksize = (int*) CATNET_MALLOC(parpoolsize  * sizeof(int));
		}

		if (!pcats || !paridx || !blocksize) {
			CATNET_FREE(parpool);
			return 0;
		}

		blocksize[parpoolsize - 1] = 1;
		for (i = parpoolsize - 2; i >= 0; i--) {
			blocksize[i] = blocksize[i + 1] * m_numCategories[parpool[i + 1]];
		}

		jointprobsize = 1;
		for (i = 0; i < parpoolsize; i++) {
			jointprobsize *= m_numCategories[parpool[i]];
		}

		jointprob = (t_prob*) CATNET_MALLOC(jointprobsize * sizeof(t_prob));
		if (!jointprob) {
			CATNET_FREE(pcats);
			CATNET_FREE(paridx);
			CATNET_FREE(blocksize);
			CATNET_FREE(parpool);
			return 0;
		}

		for (i = 0; i < jointprobsize; i++)
			jointprob[i] = 1.0;

		for (ipool = 0; ipool < parpoolsize; ipool++) {
			poolnode = parpool[ipool];
			probnode = m_pProbLists[poolnode];
			if (m_numParents[poolnode] == 0) {
				for (ii = 0; ii < jointprobsize; ii += (blocksize[ipool]
						* m_numCategories[poolnode])) {
					for (ic = 0; ic < m_numCategories[poolnode]; ic++) {
						i0 = ic * blocksize[ipool];
						for (i = 0; i < blocksize[ipool]; i++) {
							jointprob[ii + i0 + i] *= probnode->pProbs[ic];
						}
					}
				}
				continue;
			}

			memset(paridx, 0, parpoolsize * sizeof(int));
			for (ipar = 0; ipar < m_numParents[poolnode]; ipar++) {
				par = m_parents[poolnode][ipar];
				for (i = 0; i < ipool; i++)
					if (par == parpool[i])
						break;
				if (i >= ipool) {
					// bad
					break;
				}
				/* save the node's parent index */
				paridx[i] = ipar + 1;
			}

			if (m_maxParents > 0) {
				memset(pcats, 0, m_maxParents * sizeof(int));
			}
			for (ii = 0; ii < jointprobsize; ii += (blocksize[ipool] * m_numCategories[poolnode])) {

				for (j = 0; j < ipool; j++) {
					if (paridx[j] > 0) {
						ic = (int) (ii / blocksize[j]);
						cc = m_numCategories[m_parents[poolnode][paridx[j] - 1]];
						ic -= cc * (int) (ic / cc);
						pcats[paridx[j] - 1] = ic;
						if(pcats[paridx[j] - 1] >= m_numCategories[parpool[j]])
							pcats[paridx[j] - 1] = m_numCategories[parpool[j]] - 1;
					}
				}

				parprob = probnode->find_slot(0, pcats, 0);
				if (!parprob)
					continue;
				for (ic = 0; ic < m_numCategories[poolnode]; ic++) {
					i0 = ic * blocksize[ipool];
					for (i = 0; i < blocksize[ipool]; i++) {
						jointprob[ii + i0 + i] *= parprob[ic];
					}
				}
			}
		}

		CATNET_FREE(pcats);
		CATNET_FREE(paridx);
		CATNET_FREE(blocksize);
		CATNET_FREE(parpool);

		return jointprob;
	}

	t_prob *marginal_prob(int nnode) {

		t_prob *jointprob, *margprob;
		int jointprobsize = 0, nodecats, ic, k;
		nodecats = m_numCategories[nnode];
		margprob  = (t_prob*) CATNET_MALLOC(nodecats * sizeof(t_prob));
		if (!margprob)
			return 0;
		jointprob = findJointProb(nnode, jointprobsize);
		if (!jointprob)
			return 0;
		for (ic = 0; ic < nodecats; ic++) {
			k = ic;
			margprob[ic] = 0;
			while (k < jointprobsize) {
				margprob[ic] += jointprob[k];
				k += nodecats;
			}
		}
		CATNET_FREE(jointprob);

		return margprob;
	}

	int *getOrder() {
		int *porder = 0, *pfound = 0, bhaspar, i, j, k;
		if(m_numNodes < 1 || !m_numParents || !m_parents)
			return 0;

		porder = (int*)CATNET_MALLOC(m_numNodes * sizeof(int));
		pfound = (int*)CATNET_MALLOC(m_numNodes * sizeof(int));
		if (!porder || !pfound)
			return 0;
		memset(pfound, 0, m_numNodes * sizeof(int));
		for(i = 0; i < m_numNodes; i++) {
			for(j = 0; j < m_numNodes; j++) {
				if(pfound[j])
					continue;
				if(m_numParents[j] <= 0)
					break;
				bhaspar = 0;
				for(k = 0; k < m_numParents[j]; k++) {
					if(pfound[m_parents[j][k]])
						continue;
					bhaspar = 1;
					break;
				}
				if(!bhaspar)
					break;
			}
			if(j >= m_numNodes || pfound[j]) {
				CATNET_FREE(pfound);
				CATNET_FREE(porder);
				return 0;
			}
			porder[i] = j;
			pfound[j] = 1;
		}
		CATNET_FREE(pfound);
		return porder;
	}

	/*
		virtual extensions: sample related functions 
	*/
public:
	virtual t_prob nodeSampleLoglik(int nnode, int *pnodepars, int nodepars, t_sample *psamples, int nsamples, int klmode = 0) {
		return 0;
	};
	virtual t_prob sampleLoglik(t_sample *psamples, int nsamples) {
		return 0;
	};
	virtual t_prob* sampleLoglikVector(t_sample *psamples, int nsamples, int *pert=0, int klmode = 0) {
		return 0;
	};
	virtual t_prob* bySampleLoglikVector(t_sample *psamples, int nsamples, int *pert=0, int klmode = 0) {
		return 0;
	};
	virtual t_prob sampleNodeLoglik(int nnode, t_sample *psamples, int nsamples, int klmode = 0) {
		return 0;
	};
	virtual t_prob setNodeSampleProb(int nnode, t_sample *psamples, int nsamples, int bNormalize = 0) {
		return 0;
	};
	virtual t_prob setNodeSampleProbKL(int nnode, t_sample *psamples, int nsamples, int *pClasses, int bNormalize = 0, int bUsePearson = 0) {
		return 0;
	};

};

#endif /* CATNET_BASE_H_ */

