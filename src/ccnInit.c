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


#include "catnet_rexport.h"

int g_setseed = -1;
size_t g_memcounter = 0;

static const R_CallMethodDef R_CallDef[] = {
	{"ccnSearchOrder", (DL_FUNC)&searchOrder, 15},
	{"ccnSearchSA", (DL_FUNC)&searchSA, 23},
	{"ccnParHistogram", (DL_FUNC)&catnetParHistogram, 14},
	{"ccnSetProb", (DL_FUNC)&catnetSetProb, 3},
	{"ccnLoglik", (DL_FUNC)&catnetLoglik, 4},
	{"ccnNodeLoglik", (DL_FUNC)&catnetNodeLoglik, 5},
	{"ccnEntropyPairwise", (DL_FUNC)&catnetEntropyPairwise, 2},
	{"ccnEntropyOrder", (DL_FUNC)&catnetEntropyOrder, 2},
	{"ccnKLPairwise", (DL_FUNC)&catnetKLpairwise, 2},
	{"ccnPearsonPairwise", (DL_FUNC)&catnetPearsonPairwise, 2},
	{"ccnMarginalProb", (DL_FUNC)&catnetMarginalProb, 2},
	{"ccnReleaseCache", (DL_FUNC)&catnetReleaseCache, 0},
	{"ccnSamples", (DL_FUNC)&catnetSamples, 4},
	{"ccnSetSeed", (DL_FUNC)&catnetSetSeed, 1},
	{"ccnCatnetFromDagEvaluate", (DL_FUNC)&createCatnetFromDagEvaluate, 2},
	{"ccnSoftQuant", (DL_FUNC)&catnetSoftQuant, 8},
	{NULL, NULL, 0},
};

void R_init_sdnet(DllInfo *info)
{
	R_registerRoutines(info,NULL,R_CallDef,NULL,NULL);
	R_useDynamicSymbols(info, TRUE);
	g_memcounter = 0;
	// to disable C stack limit
	//R_CStackLimit = (uintptr_t)-1;
}

void R_unload_sdnet(DllInfo *info)
{
	catnetReleaseCache();
}
