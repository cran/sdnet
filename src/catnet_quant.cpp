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
 * catnet_quant.cpp
 *
 *  Created on: July 12, 2012
 *      Author: nbalov
 */

#include "utils.h"

extern "C" {

extern int g_setseed;
extern size_t g_memcounter;

/* general gaussian mixture */
inline void _gmm(double *pdata, double *pSampleWeights, int ndata, int numCats, double *pmu, double *psig, double *pw, double cover=0.95, int maxiter=8, double feps=CATNET_EPS) {

	int c, i, iter;
	double *qvals, *poldw, *pwc, ff, fsum;

	if (!pdata || !pSampleWeights || !pmu || !psig || !pw) {
		return;
	}

	if(numCats < 1)
		numCats = 1;
	if(maxiter < 2)
		maxiter = 2;
	if(cover < 0 || cover > 1)
		cover = 1;

	poldw = (double*)CATNET_MALLOC((numCats)      *sizeof(double));
	qvals = (double*)CATNET_MALLOC((numCats+1)    *sizeof(double));
	pwc   = (double*)CATNET_MALLOC((ndata*numCats)*sizeof(double));
	if (!poldw || !qvals || !pwc) {
		if (qvals)
			CATNET_FREE(qvals);
		if (poldw)
			CATNET_FREE(poldw);
		if (pwc)
			CATNET_FREE(pwc);
		return;
	}
	
	double fmin, fmax, rang;
	fmin = FLT_MAX;
	fmax = -FLT_MAX;
	for(i =0; i < ndata; i++) {
		if(fmin > pdata[i])
			fmin = pdata[i];
		if(fmax < pdata[i])
			fmax = pdata[i];
	}
	rang = fmax - fmin;
	for(i = 0; i < numCats; i++) {
		c = (int)(ndata*i/numCats);
		qvals[i] = fmin + rang*((1-cover)/2 + (cover*i)/numCats);
	}
	qvals[numCats] = fmax;

	for(i = 0; i < numCats; i++) {
		pmu[i]  = (qvals[i+1] + qvals[i])/2;
		psig[i] = (qvals[i+1] - qvals[i])/2;
		pw[i]   = 1/(double)numCats;
	}

	for(iter=0; iter < maxiter; iter++) {
		memcpy(poldw, pw, numCats*sizeof(double));
		for(c = 0; c < numCats; c++) {
			fsum = 1/sqrt(psig[c]);
			for(i = 0; i < ndata; i++) {
				ff = fsum*(pdata[i]-pmu[c]);
				pwc[i*numCats+c] = pw[c]*exp(-0.5*ff*ff)*fsum;
			}
		}
		for(i = 0; i < ndata; i++) {
			fsum = 0;
			for(c = 0; c < numCats; c++)
				fsum += pwc[i*numCats+c];
			if(fsum < CATNET_EPS)
				fsum = CATNET_EPS;
			fsum = 1/fsum;
			for(c = 0; c < numCats; c++)
				pwc[i*numCats+c] *= fsum;
		}
		for(c = 0; c < numCats; c++) {
			pmu[c] = 0;
			fsum = 0;
			for(i = 0; i < ndata; i++) {
				pmu[c] += pwc[i*numCats+c] * pdata[i];
				fsum += pwc[i*numCats+c]; 
			}
			if(fsum < CATNET_EPS)
				fsum = CATNET_EPS;
			pmu[c] /= fsum;
			psig[c] = 0;
			for(i = 0; i < ndata; i++) {
				ff = pdata[i]-pmu[c];
				psig[c] += pwc[i*numCats+c]*ff*ff;
			}
			psig[c] /= fsum;
			if(psig[c] < CATNET_EPS)
				psig[c] = CATNET_EPS;
			pw[c] = 0;
			for(i = 0; i < ndata; i++)
				pw[c] += pwc[i*numCats+c];
			pw[c] /= ndata;
		}

		fsum = 0;
		for(c = 0; c < numCats; c++) 
			fsum += psig[c];
		fsum /= numCats;
		for(c = 0; c < numCats; c++) 
			psig[c] = fsum;

		fsum = 0;
		for(c = 0; c < numCats; c++) {
			ff = poldw[c]-pw[c];
			fsum += ff*ff;
		}
		if(fsum < feps)
			break;
	}

	for(c = 0; c < numCats; c++)
		psig[c] = sqrt(psig[c]);

	CATNET_FREE(qvals);
	CATNET_FREE(poldw);
	CATNET_FREE(pwc);
}

inline void _qmm(double *pdata, double *pSampleWeights, int ndata, int numCats, double *pmu, double *psig, double *pw, double cover=0.95, int maxiter=8, double feps=CATNET_EPS) {

	int c, i, j, *pranks, *pind, iter;
	double *pwc, *poldw, ff, fsum;

	if (!pdata || !pSampleWeights || !pmu || !psig || !pw) {
		return;
	}

	if(numCats < 1)
		numCats = 1;
	if(maxiter < 2)
		maxiter = 2;
	if(cover < 0 || cover > 1)
		cover = 1;

	pwc    = (double*)CATNET_MALLOC((ndata*numCats)*sizeof(double));
	poldw  = (double*)CATNET_MALLOC((numCats)      *sizeof(double));
	pranks =    (int*)CATNET_MALLOC(ndata          *sizeof(int));
	pind   =    (int*)CATNET_MALLOC(ndata          *sizeof(int));
	if (!poldw || !pwc || !pranks || !pind) {
		if (poldw)
			CATNET_FREE(poldw);
		if (pwc)
			CATNET_FREE(pwc);
		if (pranks)
			CATNET_FREE(pranks);
		if (pind)
			CATNET_FREE(pind);
		return;
	}

	for(i =0; i < ndata; i++) {
		pranks[i] = 0;
		for(j = 0; j < ndata; j++)
			if(pdata[j] < pdata[i])
				pranks[i]++;		
	}
	for(i = 0; i < ndata; i++) {
		for(j = 0; j < ndata; j++) {
			if(pranks[j] == i) {
				pind[i] = j++;
				for(; j < ndata; j++)
					if(pranks[j] == i)
						pranks[j]++;
				break;
			}
		}
	}
	for(i = 0; i < numCats; i++) {
		//c = (int)(ndata*(i+1)/(numCats+1));
		fsum = 0;
		for(c = 0; c < ndata; c++) {
			fsum += pSampleWeights[pind[c]];
			if(fsum >= (double)(i+1)/(double)(numCats+1))
				break;
		}
		if(c < 0) c = 0;
		if(c >= ndata) c = ndata-1;
		pmu[i] = pdata[pind[c]];
	}
	CATNET_FREE(pranks);
	CATNET_FREE(pind);

	for(i = 0; i < numCats; i++) {
		psig[i] = (pmu[numCats-1] - pmu[0])/2;
		psig[i] = psig[i]*psig[i];
		pw[i] = 1/(double)numCats;
	}

	for(iter=0; iter < maxiter; iter++) {

		memcpy(poldw, pw, numCats*sizeof(double));

		for(c = 0; c < numCats; c++) {
			fsum = 1/sqrt(psig[c]);
			for(i = 0; i < ndata; i++) {
				ff = fsum*(pdata[i]-pmu[c]);
				pwc[i*numCats+c] = pw[c]*exp(-0.5*ff*ff)*fsum;
			}
		}
		for(i = 0; i < ndata; i++) {
			fsum = 0;
			for(c = 0; c < numCats; c++)
				fsum += pwc[i*numCats+c];
			if(fsum < CATNET_EPS)
				fsum = CATNET_EPS;
			fsum = 1/fsum;
			for(c = 0; c < numCats; c++)
				pwc[i*numCats+c] *= fsum;
		}
		for(c = 0; c < numCats; c++) {
			fsum = 0;
			for(i = 0; i < ndata; i++)
				fsum += pwc[i*numCats+c]; 
			psig[c] = 0;
			for(i = 0; i < ndata; i++) {
				ff = pdata[i]-pmu[c];
				psig[c] += pwc[i*numCats+c]*ff*ff;
			}
			if(fsum < CATNET_EPS)
				fsum = CATNET_EPS;
			psig[c] /= fsum;
			if(psig[c] < CATNET_EPS)
				psig[c] = CATNET_EPS;
			pw[c] = 0;
			for(i = 0; i < ndata; i++)
				pw[c] += pwc[i*numCats+c];
			pw[c] /= ndata;
		}
		fsum = 0;
		for(c = 0; c < numCats; c++) 
			fsum += psig[c];
		fsum /= numCats;
		for(c = 0; c < numCats; c++) 
			psig[c] = fsum;

		fsum = 0;
		for(c = 0; c < numCats; c++) {
			ff = poldw[c]-pw[c];
			fsum += ff*ff;
		}
		if(fsum < feps)
			break;

	} // iter

	for(c = 0; c < numCats; c++)
		psig[c] = sqrt(psig[c]);

	CATNET_FREE(poldw);
	CATNET_FREE(pwc);
}

/* uniform gaussian mixture */
inline void _umm(double *pdata, double *pSampleWeights, int ndata, int numCats, double *pmu, double *psig, double *pw, double cover=0.95, int maxiter=8, double feps=CATNET_EPS) {

	int c, i;
	double *qvals, *poldw, *pwc, ff, fsum;

	if (!pdata || !pSampleWeights || !pmu || !psig || !pw) {
		return;
	}

	if(numCats < 1)
		numCats = 1;
	if(maxiter < 2)
		maxiter = 2;
	if(cover < 0 || cover > 1)
		cover = 1;

	poldw = (double*)CATNET_MALLOC((numCats)      *sizeof(double));
	qvals = (double*)CATNET_MALLOC((numCats+1)    *sizeof(double));
	pwc   = (double*)CATNET_MALLOC((ndata*numCats)*sizeof(double));
	if (!poldw || !qvals || !pwc) {
		if (qvals)
			CATNET_FREE(qvals);
		if (poldw)
			CATNET_FREE(poldw);
		if (pwc)
			CATNET_FREE(pwc);
		return;
	}

	double fmin, fmax, rang;
	fmin = FLT_MAX;
	fmax = -FLT_MAX;
	for(i =0; i < ndata; i++) {
		if(fmin > pdata[i])
			fmin = pdata[i];
		if(fmax < pdata[i])
			fmax = pdata[i];
	}
	rang = fmax - fmin;
	for(i = 0; i <= numCats; i++) 
		qvals[i] = fmin + rang*((1-cover)/2 + (cover*i)/numCats);

	for(i = 0; i < numCats; i++) {
		pmu[i] = (qvals[i+1] + qvals[i])/2;
		psig[i] = (qvals[i+1] - qvals[i])/2;
		pw[i] = 1/(double)numCats;
	}

	for(c = 0; c < numCats; c++) {
		fsum = 1/sqrt(psig[c]);
		for(i = 0; i < ndata; i++) {
			ff = fsum*(pdata[i]-pmu[c]);
			pwc[i*numCats+c] = pw[c]*exp(-0.5*ff*ff);
		}
	}
	for(i = 0; i < ndata; i++) {
		fsum = 0;
		for(c = 0; c < numCats; c++)
			fsum += pwc[i*numCats+c];
		if(fsum < CATNET_EPS)
			fsum = CATNET_EPS;
		fsum = 1/fsum;
		for(c = 0; c < numCats; c++)
			pwc[i*numCats+c] *= fsum;
	}
	for(c = 0; c < numCats; c++) {
		psig[c] = 0;
		for(i = 0; i < ndata; i++) {
			ff = pdata[i]-pmu[c];
			psig[c] += pwc[i*numCats+c]*ff*ff;
		}
		if(fsum < CATNET_EPS)
			fsum = CATNET_EPS;
		psig[c] /= fsum;
		if(psig[c] < CATNET_EPS)
			psig[c] = CATNET_EPS;
		pw[c] = 0;
		for(i = 0; i < ndata; i++)
			pw[c] += pwc[i*numCats+c];
		pw[c] /= ndata;
	}

	fsum = 0;
	for(c = 0; c < numCats; c++) 
		fsum += psig[c];
	fsum /= numCats;
	for(c = 0; c < numCats; c++) 
		psig[c] = fsum;
	for(c = 0; c < numCats; c++)
		psig[c] = sqrt(psig[c]);

	CATNET_FREE(qvals);
	CATNET_FREE(poldw);
	CATNET_FREE(pwc);
}

SEXP catnetSoftQuant(SEXP rSamples, SEXP rSampleWeights, SEXP rNumCats, SEXP rLearnset, SEXP rCover, SEXP rMode, SEXP rMaxiter, SEXP rEps) {

	int gmode, nnode, numSamples, numNodes, *pNumCats, numCats, c, j, maxiter, nLearnset, ndataline, *pLearnset;
	double *pSamples, *pSampleWeights, *prow, *pdata, fCover, fsum, *pallmu, *pmu, *pallsig, *psig, *pallw, *pw, pp, pmax, feps; 
	int *ddata, cmax;
	SEXP dim;
	SEXP plist, rvec = R_NilValue;

	if(!isMatrix(rSamples))
		error("Data is not a matrix");

	dim        = GET_DIM(rSamples);
	numNodes   = INTEGER(dim)[0];
	numSamples = INTEGER(dim)[1];

	nLearnset = LENGTH(rLearnset);	

	pNumCats = (int*)CATNET_MALLOC(numNodes*sizeof(int));
	if (!pNumCats) {
		return rvec;
	}
	PROTECT(rNumCats = AS_INTEGER(rNumCats));
	memcpy(pNumCats, INTEGER_POINTER(rNumCats), numNodes*sizeof(int));
	UNPROTECT(1);
	numCats = 0;
	for(j = 0; j < numNodes; j++) {
		if(pNumCats[j]<0)
			error("Wrong numcats");
		if(numCats < pNumCats[j])
			numCats = pNumCats[j];
	}	

	prow     = (double*)CATNET_MALLOC(nLearnset*sizeof(double));
	pdata    = (double*)CATNET_MALLOC(numSamples*numNodes*numCats*sizeof(double));
	ddata    =    (int*)CATNET_MALLOC(numSamples*numNodes*sizeof(int));
	pLearnset=    (int*)CATNET_MALLOC(nLearnset*sizeof(int));

	pallw   = (double*)CATNET_MALLOC(numCats*numNodes*sizeof(double));
	pallmu  = (double*)CATNET_MALLOC(numCats*numNodes*sizeof(double));
	pallsig = (double*)CATNET_MALLOC(numCats*numNodes*sizeof(double));

	pmu  = (double*)CATNET_MALLOC(numCats*sizeof(double));
	psig = (double*)CATNET_MALLOC(numCats*sizeof(double));
	pw   = (double*)CATNET_MALLOC(numCats*sizeof(double));

	pSampleWeights = (double*)CATNET_MALLOC(nLearnset*sizeof(double)); 

	if (!prow || !pdata || !ddata || !pLearnset || !pallw || !pallmu ||
            !pallsig || !pmu || !psig || !pw) {
		if (pNumCats)
			CATNET_FREE(pNumCats);
		if (pSampleWeights)
			CATNET_FREE(pSampleWeights);
		if (pLearnset)
			CATNET_FREE(pLearnset);
		if (pw)
			CATNET_FREE(pw);
		if (psig)
			CATNET_FREE(psig);
		if (pmu)
			CATNET_FREE(pmu);
		if (prow)
			CATNET_FREE(prow);
		if (pallw)
			CATNET_FREE(pallw);
		if (pallmu)
			CATNET_FREE(pallmu);
		if (pallsig)
			CATNET_FREE(pallsig);
		if (pdata)
			CATNET_FREE(pdata);
		if (ddata)
			CATNET_FREE(ddata);
		return rvec;
	}

	PROTECT(rLearnset = AS_INTEGER(rLearnset));
	memcpy(pLearnset, INTEGER_POINTER(rLearnset), nLearnset*sizeof(int));
	UNPROTECT(1);
	for(j = 0; j < nLearnset; j++) {
		pLearnset[j]--;
		if(pLearnset[j]<0 || pLearnset[j]>=numSamples)
			error("Wrong learnset");
	}

	PROTECT(rCover = AS_NUMERIC(rCover));
	fCover = NUMERIC_POINTER(rCover)[0];
	UNPROTECT(1);

	PROTECT(rEps = AS_NUMERIC(rEps));
	feps = NUMERIC_POINTER(rEps)[0];
	UNPROTECT(1);

	PROTECT(rMode = AS_CHARACTER(rMode));
	gmode = 0;
	if (!strncmp(CHARACTER_VALUE(rMode), "gauss",    5))
		gmode = 1;
	if (!strncmp(CHARACTER_VALUE(rMode), "uniform",  7))
		gmode = 2;
	if (!strncmp(CHARACTER_VALUE(rMode), "quantile", 8))
		gmode = 3;
	UNPROTECT(1);

	PROTECT(rMaxiter = AS_INTEGER(rMaxiter));
	maxiter = INTEGER_POINTER(rMaxiter)[0];
	UNPROTECT(1);
	if(maxiter < 1) maxiter = 1;
	if(maxiter > 1000) maxiter = 1000;


	for(j = 0; j < nLearnset; j++) 
		pSampleWeights[j] = 1/(double)nLearnset; 
	if(nLearnset == LENGTH(rSampleWeights)) { 
		PROTECT(rSampleWeights = AS_NUMERIC(rSampleWeights)); 
		memcpy(pSampleWeights, NUMERIC_POINTER(rSampleWeights), nLearnset*sizeof(double)); 
		UNPROTECT(1); 
	}

	PROTECT(rSamples = AS_NUMERIC(rSamples));
	pSamples = NUMERIC_POINTER(rSamples);

	for(nnode = 0; nnode < numNodes; nnode++) {

		for(j = 0; j < nLearnset; j++)
			prow[j] = pSamples[pLearnset[j]*numNodes + nnode];

		switch(gmode) {
		case 2: _umm(prow, pSampleWeights, nLearnset, pNumCats[nnode], pmu, psig, pw, fCover, maxiter, feps);
			break;
		case 3: _qmm(prow, pSampleWeights, nLearnset, pNumCats[nnode], pmu, psig, pw, fCover, maxiter, feps);
			break;
		default:
			_gmm(prow, pSampleWeights, nLearnset, pNumCats[nnode], pmu, psig, pw, fCover, maxiter, feps);
		}

		fsum = 0;
		for(c = 0; c < numCats; c++)
			fsum += psig[c];
		fsum /= numCats;
		for(c = 0; c < numCats; c++) {
			/* equal sigmas */
			//psig[c] = fsum;
			pallmu[numCats*nnode+c] = pmu[c];
			pallsig[numCats*nnode+c] = psig[c];
			psig[c] = 1 / psig[c];
			pallw[numCats*nnode+c] = pw[c];
		}

		ndataline = nnode*numCats;
		for(j = 0; j < numSamples; j++) {
			fsum = 0;
			for(c = 0; c < numCats; c++) {
				pp = (pSamples[j*numNodes + nnode] - pmu[c])*psig[c];
				pp = psig[c]*exp(-0.5*pp*pp);
				pdata[ndataline + c] = pp;
				fsum += pp;
			}
			fsum = 1/fsum;
			
			pmax = -1;
			cmax = 0;
			for(c = 0; c < numCats; c++) {
				pdata[ndataline + c] *= fsum;
				if(pmax < pdata[ndataline + c]) {
					pmax = pdata[ndataline + c];
					cmax = c;
				}
			}
			ddata[j*numNodes + nnode] = cmax + 1;
			ndataline += numCats*numNodes;
		}

	}

	UNPROTECT(1); // rSamples

	CATNET_FREE(pSampleWeights);
	CATNET_FREE(pLearnset);
	CATNET_FREE(pNumCats);
	CATNET_FREE(pw);
	CATNET_FREE(psig);
	CATNET_FREE(pmu);
	CATNET_FREE(prow);

	PROTECT(plist = allocVector(VECSXP, 5));

	PROTECT(rvec = NEW_NUMERIC(numSamples*numNodes*numCats));
	memcpy(NUMERIC_POINTER(rvec), pdata, numSamples*numNodes*numCats*sizeof(double));
	SET_VECTOR_ELT(plist, 0, rvec);
	UNPROTECT(1);//rvec

	PROTECT(rvec = NEW_INTEGER(numSamples*numNodes));
	memcpy(INTEGER_POINTER(rvec), ddata, numSamples*numNodes*sizeof(int));
	SET_VECTOR_ELT(plist, 1, rvec);
	UNPROTECT(1);//rvec

	PROTECT(rvec = NEW_NUMERIC(numCats*numNodes));
	memcpy(NUMERIC_POINTER(rvec), pallmu, numCats*numNodes*sizeof(double));
	SET_VECTOR_ELT(plist, 2, rvec);
	UNPROTECT(1);//rvec

	PROTECT(rvec = NEW_NUMERIC(numCats*numNodes));
	memcpy(NUMERIC_POINTER(rvec), pallsig, numCats*numNodes*sizeof(double));
	SET_VECTOR_ELT(plist, 3, rvec);
	UNPROTECT(1);//rvec

	PROTECT(rvec = NEW_NUMERIC(numCats*numNodes));
	memcpy(NUMERIC_POINTER(rvec), pallw, numCats*numNodes*sizeof(double));
	SET_VECTOR_ELT(plist, 4, rvec);
	UNPROTECT(1);//rvec

	UNPROTECT(1);//plist

	CATNET_FREE(pallw);
	CATNET_FREE(pallmu);
	CATNET_FREE(pallsig);
	CATNET_FREE(pdata);
	CATNET_FREE(ddata);

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//Rprintf(str);

	return plist;

}

} // extern "C"
