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
 * utils.c
 *
 *  Created on: Nov 16, 2009
 *      Author: Nikolay Balov
 */

#include "utils.h"

extern size_t g_memcounter;

void * CATNET_MALLOC(size_t nsize) { 
	if(nsize <= 0)
		return 0;
	void *pMem;
	pMem = malloc(nsize);
	if(!pMem) {
		error("Insufficient memory");
		return 0;
	}
	return pMem;
	// memmonitor
	g_memcounter += nsize;
	pMem = malloc(sizeof(int) + nsize);
	if(!pMem) {
	char str[128];
	sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	Rprintf(str);
		error("Insufficient memory");
		return 0;
	}
	*(int*)pMem = nsize;
	pMem = (void*)((int*)pMem + 1);
	return pMem;
}

void CATNET_FREE(void *pMem) {
	if(!pMem)	
		return;
	free(pMem);
	return;
	// memmonitor
	pMem = (void*)((int*)pMem-1);
	size_t nsize = *((int*)pMem);
	g_memcounter -= nsize;
	free(pMem);
}

void CATNET_MEM_ERR() {
	// generate R-error
	error("Insufficient memory");
}

void CATNET_FORMAT_ERR() {
	// generate R-error
	error("Wrong data format");
}

void CATNET_PARAM_ERR() {
	// generate R-error
	error("Wrong parameters");
}

void CATNET_NOTSUPP_ERR() {
	// generate R-error
	error("Unsupported parameters");
}

void CATNET_ERR(const char *str) {
	error(str);
}

void CATNET_WARNING(const char *str) {
	warning(str);
}

