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
 * dag_list.h
 * includes soft edge constraints in addition to the original CATNET_SEARCH
 *
 *  Created on: Nov 8, 2011
 *      Author: nbalov
 */

#include "utils.h"
#include "search_params.h"

#ifndef DAG_LIST_H
#define DAG_LIST_H

#define BITS_PER_INDEX(nbits, maxpars) { \
		switch(maxpars) {\
		case 0:\
		case 1:	nbits = 1; break;\
		case 2:\
		case 3: nbits = 2; break;\
		case 4:\
		case 5:\
		case 6:\
		case 7: nbits = 4; break;\
		default:if(maxpars < 256) nbits = 8; \
			else nbits = 32; break;\
		}\
}

template<class t_prob>
struct DAG_PARS {
	/* each DAG has numNodes of these */
	int 	numNodes;
	int	nbits;
	int	numParSize;
	char	*numPars;
	t_prob	loglik;
	int	complx;
	DAG_PARS<t_prob> *next;

	/*void *operator new(size_t size) {
		void *pobj = CATNET_MALLOC(size);
		return(pobj);
	}
	void operator delete(void *pobj) {
		CATNET_FREE(pobj);
	}*/

	~DAG_PARS() {
		if(numPars)
			CATNET_FREE(numPars);
		if(next)
			delete next;
		next = 0;
	}
	DAG_PARS() {
		numNodes = 0;
		nbits = 0;
		numPars = 0;
		numParSize = 0;
		loglik = 0;
		complx = 0;
		next = 0;
	}
	DAG_PARS(int nNumNodes, int maxpars) {
		numNodes = nNumNodes;
		BITS_PER_INDEX(nbits, maxpars);
		numParSize = (int)(1+(numNodes*nbits/8));
		numPars = (char*)CATNET_MALLOC(numParSize*sizeof(char));
		if (numPars)
			memset(numPars, 0, numParSize*sizeof(char));
		loglik = 0;
		complx = 0;
		next = 0;
	}
	DAG_PARS(const DAG_PARS<t_prob> *pdag) {
		if (!pdag)
			return;
		numNodes = pdag->numNodes;
		nbits = pdag->nbits;
		numParSize = pdag->numParSize;
		numPars = (char*)CATNET_MALLOC(numParSize*sizeof(char));
		if (numPars && pdag->numPars)
			memcpy(numPars, pdag->numPars, numParSize*sizeof(char));
		loglik = pdag->loglik;
		complx = pdag->complx;
		next = 0;
	}
	void copyNumPars(const DAG_PARS<t_prob> *pdag) {
		if (!pdag)
			return;
		if (numPars && pdag->numPars)
			memcpy(numPars, pdag->numPars, numParSize*sizeof(char));
	}
	int getNumPars(int ind) {
		if (!numPars)
			return 0;
		int nbyte = (int)(ind*nbits/8);
		int noff = ind*nbits-8*nbyte;
		if(nbits <= 8)
			return((numPars[nbyte] >> noff)&((1<<nbits)-1));
		else 
			return(*((int*)(numPars+nbyte)));
	}
	void setNumPar(int ind, int npars) {
		if (!numPars)
			return;
		int nbyte = (int)(ind*nbits/8);
		int noff  = ind*nbits-8*nbyte;
		if(nbits <= 8) 
			numPars[nbyte] = (numPars[nbyte] & ~(((1<<nbits)-1) << noff)) 
					| (npars&((1<<nbits)-1)) << noff;
		else {
			*((int*)(numPars+nbyte)) = npars;
		}
	}
	int compressNumPars(int *pIntBuff, char *pByteBuff, int &nBuffSize, int *pOrder/*[1:numNodes]*/) {
		int nbyte, noff, i, j, k;

		if(!numPars || nBuffSize < numNodes || !pIntBuff || !pByteBuff || !pOrder)
			return 0;

		nbyte = 0;
		noff  = 0;

		if(nbits <= 8) {
			for(i = 0; i < numNodes; i++) {
				pIntBuff[pOrder[i]-1] = (int)((numPars[nbyte] >> noff)&((1<<nbits)-1));
				noff += nbits;
				if(noff >= 8) {
					noff = 0;
					nbyte++;
				}
			}
		}
		else 
			for(i = 0; i < numNodes; i++) 
				pIntBuff[pOrder[i]-1] = *((int*)numPars+i);
		*((int*)pByteBuff) = nbits;
		if(nbits <= 8) {
			nBuffSize = (int)sizeof(int);
			i = 0;
			while(i < numNodes) {
				k = pIntBuff[i];
				j = 1;
				while(i+j < numNodes && pIntBuff[i+j] == k && j < 127) j++;
				if(j <= 2) {
					while(j-- > 0)
						pByteBuff[nBuffSize++] = pIntBuff[i++];
					continue;
				}
				pByteBuff[nBuffSize++] = (char)-j;
				pByteBuff[nBuffSize++] = (char)k;
				i += j;
			}
		}
		else {
			int *pout = (int*)pByteBuff;
			nBuffSize = 1;
			i = 0;
			while(i < numNodes) {
				k = pIntBuff[i];
				j = 1;
				while(i+j < numNodes && pIntBuff[i+j] == k && j < (1<<(8*sizeof(int)-2))) j++;
				if(j <= 2) {
					while(j-- > 0)
						pout[nBuffSize++] = pIntBuff[i++];
					continue;
				}
				pout[nBuffSize++] = -j;
				pout[nBuffSize++] = k;
				i += j;
			}
			nBuffSize *= (int)sizeof(int);
		}

		return nBuffSize;
	}
};

template<class t_prob, class t_ind>
struct DAG_LIST {
	int	m_numNodes;
	int	*m_numParSlots;
	/* each node has m_numParSlots*n_maxpars of t_ind */
	t_ind	**m_parSlots;
	/* each node has m_numParSlots of t_prob */
	t_prob	**m_parLogliks;
	t_ind	**m_parComplx;
	t_ind	**m_parSampleSize;
	int	m_numDags;
	DAG_PARS<t_prob> *m_dagPars; 

	/* use a virtual distructor */
	virtual ~DAG_LIST() {
		reset();
	}

	DAG_LIST() {
		m_numNodes    = 0;
		m_numParSlots = 0;
		m_parSlots    = 0;
		m_parLogliks  = 0;
		m_parComplx   = 0;
		m_parSampleSize = 0;
		m_numDags = 0;
		m_dagPars = 0;
	}

	void reset() {
		for(int i = 0; i < m_numNodes; i++) {
			if(m_parSlots && m_parSlots[i])
				CATNET_FREE(m_parSlots[i]);		
			if(m_parLogliks && m_parLogliks[i])
				CATNET_FREE(m_parLogliks[i]);
			if(m_parComplx && m_parComplx[i])
				CATNET_FREE(m_parComplx[i]);
			if(m_parSampleSize && m_parSampleSize[i])
				CATNET_FREE(m_parSampleSize[i]);
		}
		if(m_numParSlots)
			CATNET_FREE(m_numParSlots);
		if(m_parSlots)
			CATNET_FREE(m_parSlots);
		if(m_parLogliks)
			CATNET_FREE(m_parLogliks);
		if(m_parComplx)
			CATNET_FREE(m_parComplx);
		if(m_parSampleSize)
			CATNET_FREE(m_parSampleSize);
		if(m_dagPars)
			delete m_dagPars;
		m_numNodes    = 0;
		m_numParSlots = 0;
		m_parSlots    = 0;
		m_parLogliks  = 0;
		m_parComplx   = 0;
		m_parSampleSize = 0;
		m_numDags = 0;
		m_dagPars = 0;
	}

	virtual int search(SEARCH_PARAMETERS *pestim) = 0;
};

#endif /* DAG_LIST_H */
