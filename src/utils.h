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
 * utils.h
 *
 *  Created on: Nov 16, 2009
 *      Author: nbalov
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "thread.h"

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <float.h>
#include <time.h>
#include <stdarg.h>

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Print.h>

#ifndef CATNET_PI
#define CATNET_PI	(double)3.14159265358979323846264338327950288
#endif
#ifndef CATNET_PI2
#define CATNET_PI2	(2*(double)CATNET_PI)
#endif

#define CATNET_NAN	INT_MAX

#define CATNET_ERR_OK		0
#define CATNET_ERR_MEM		-10
#define CATNET_ERR_PARAM	-20
#define CATNET_ERR_PROC		-30

#define CATNET_EPS	(double)1E-16

void * CATNET_MALLOC(size_t nsize);
void CATNET_FREE(void *pMem);
void CATNET_FORMAT_ERR();
void CATNET_MEM_ERR();
void CATNET_PARAM_ERR();
void CATNET_NOTSUPP_ERR();
void CATNET_ERR(const char *str);
void CATNET_WARNING(const char *str);

void _combination_sets(int **&plist, int &nlist, int *curset, int *parset, int nparset, int parid, int parsize);

template<class t_elem>
void _quick_sort(t_elem *plist, int nlist){
	t_elem pivot;
	int j, nless, ngreater;
	if(nlist <= 1)
		return;
	t_elem *paux = (t_elem*)malloc(nlist*sizeof(t_elem));
	nless = 0;
	ngreater = nlist-1;
	pivot = plist[0];
	for(j = 1; j < nlist; j++) {
		if(plist[j] <= pivot) {
			paux[nless] = plist[j];
			nless++;
		}
		else {
			paux[ngreater] = plist[j];
			ngreater--;
		}
	}
	_quick_sort<t_elem>(paux, nless);
	_quick_sort<t_elem>(paux + ngreater + 1, nlist - ngreater - 1);
	paux[nless] = pivot;
	memcpy(plist, paux, nlist*sizeof(t_elem));
	free(paux);
	return;
}

template<class t_elem>
void _order(t_elem *plist, int nlist, int *porder, int decreasing = 1){
	if(!plist || !porder || nlist < 2)
		return;
	int j, i;
	t_elem *paux = (t_elem*)malloc(nlist*sizeof(t_elem));
	memcpy(paux, plist, nlist*sizeof(t_elem));

	_quick_sort<t_elem>(plist, nlist);

	for(j = 1; j < nlist; j++) 
		porder[j] = - 1;
	if(decreasing) {
		for(j = nlist - 1; j >= 0; j--) {
			for(i = 0; i < nlist; i++) {
				if(paux[i] == plist[j]) {
					porder[nlist-1-j] = i;
					paux[i] = (t_elem)FLT_MAX;
				}
			}
		}
	}
	else {
		for(j = 0; j < nlist; j++) {
			for(i = 0; i < nlist; i++) {
				if(paux[i] == plist[j]) {
					porder[j] = i;
					paux[i] = (t_elem)FLT_MAX;
				}
			}
		}
	}
	free(paux);
	return;	
}

template<class t_elem>
t_elem _gen_std_normal_var() {
	/* ISO C pseudo random generator */
	/* include stdlib.h and math.h */
	t_elem u, v;
	GetRNGstate();
	u = (t_elem)unif_rand();
	v = (t_elem)unif_rand();
	PutRNGstate();
	return(sqrt(-2*log(u)) * cos(CATNET_PI2*v));
}

template<class t_elem>
int _gen_permutation(t_elem *psample, int nsample) {
	/* psamples takes values in [1,nsample] */
	int i, j, brep, cc;
	if(nsample < 1 || !psample)
		return -1;
	int *paux = (int*)malloc(nsample*sizeof(int));
	cc = 0;
	
	GetRNGstate();
	while(++cc < 1e5) {
		brep = 0;
		for(i = 0; i < nsample; i++)
			paux[i] = (int)(nsample*nsample*unif_rand());
		for(j = 0; j < nsample; j++) {
			psample[j] = 0;
			for(i = 0; i < nsample; i++) {
				if(i!=j && paux[j] == paux[i]) {
					brep = 1;
					break;
				}
				if(paux[j] >= paux[i])
					psample[j]++;
			}
			if(brep)
				break;
		}
		if(!brep)
			break;
	}
	PutRNGstate();
	if(cc >= 1e5-1) {
		for(j = 0; j < nsample; j++)
			psample[j] = j+1;
	}
	free(paux);
	return 0;
}

template<class t_elem>
int _gen_binomial(int size, t_elem prob) {
	int i, r;
	if(size < 1 || prob <= 0)
		return 0;
	r = 0;
	GetRNGstate();
	for(i = 0; i < size; i++) {
		if((t_elem)unif_rand() <= prob)
			r++;
	}
	PutRNGstate();
	return r;
}


template<class t_prob>
t_prob _gamma_upper_bound(t_prob p, int m, t_prob delta) {
	t_prob fr1, fr2, fp, fres;
	if(p <= 0 || p >= 1 || m < 1)
		return 0;
	if(delta <= 0)
		return 1;
	fp = 1 + log((double)p);
	if(fp < 0) fp = -fp;
	fp /= delta;

	fr1 = fp + sqrt((double)(fp*fp + 2/(p*delta)));
	fr2 = 1/p + fp + sqrt((double)((1/p-fp)*(1/p-fp) + 2/(p*delta)));
	fres = fr1*fr1 + fr2*fr2;

	return(p*(1-p)*fres/(8*m));
}

/* returns increasing subsets of 'parset' of size 'parsize' */
template<class t_ind>
void _combination_sets(t_ind **&plist, int &nlist, t_ind *curset, t_ind *parset, int nparset, int parid, int parsize) 
{
	int i, ancestor;
	if(parid < 0 || parid >= parsize)
		return;
	ancestor = -1;
	if(parid > 0)
		ancestor = curset[parid-1];
	if(parid == parsize - 1) {
		for(i = 0; i < nparset; i++) {
			if(parset[i] <= ancestor)
				continue;
			t_ind **pnewlist = (t_ind**)CATNET_MALLOC((nlist+1)*sizeof(t_ind*));
			if(nlist > 0)
				memcpy(pnewlist, plist, nlist*sizeof(t_ind*));
			pnewlist[nlist] = (t_ind*)CATNET_MALLOC(parsize*sizeof(t_ind));
			if(curset)
				memcpy(pnewlist[nlist], curset, (parsize-1)*sizeof(t_ind));
			pnewlist[nlist][parsize-1] = parset[i];
			CATNET_FREE(plist);
			plist = pnewlist;
			nlist++;
		}
		if(curset) {
			CATNET_FREE(curset);
			curset = 0;
		}
		return;
	}
	for(i = 0; i < nparset; i++) {
		if(parset[i] <= ancestor)
			continue;
		t_ind *pnewset = (t_ind*)CATNET_MALLOC((parid+1)*sizeof(t_ind));
		if(curset && parid > 0)
			memcpy(pnewset, curset, parid*sizeof(t_ind));
		pnewset[parid] = parset[i];
		_combination_sets<t_ind>(plist, nlist, pnewset, parset, nparset, parid+1, parsize);
	}
	if(curset) {
		CATNET_FREE(curset);
		curset = 0;
	}
}

#endif /* UTILS_H_ */

