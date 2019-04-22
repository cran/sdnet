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

#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

/* these are called by R-functions directly */

SEXP createRCatnet(SEXP cnet);
SEXP createCatnetFromDagEvaluate(SEXP rDagEval, SEXP rDagIndex);
SEXP catnetMarginalProb(SEXP cnet, SEXP rnode);
SEXP catnetJointProb(SEXP cnet, SEXP rnode);
SEXP catnetFindParentPool(SEXP cnet, SEXP rnode);
SEXP showCatnet(SEXP cnet);
SEXP searchOrder(SEXP rSamples, SEXP rPerturbations, 
				SEXP rMaxParents, SEXP rParentSizes, 
				SEXP rMaxComplexity, SEXP rOrder, SEXP rNodeCats, 
				SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rMatEdgeLiks, 
				SEXP rUseCache, SEXP rEcho, SEXP rDagOrCatnet, SEXP rClasses);

SEXP searchSA(SEXP rNodeNames, SEXP rSamples, SEXP rPerturbations, 
			SEXP rMaxParents, SEXP rParentSizes, SEXP rMaxComplexity, SEXP rNodeCats, 
			SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rMaxParentsPool, 
			SEXP rMatEdgeLiks, SEXP rDirProbs, 
			SEXP rModel, SEXP rStartOrder,
			SEXP rTempStart, SEXP rTempCoolFact, SEXP rTempCheckOrders, 
			SEXP rMaxIter, SEXP rOrderShuffles, SEXP rStopDiff, SEXP rStopTime,  
			SEXP rThreads, SEXP rUseCache, SEXP rEcho);
SEXP catnetParHistogram(SEXP rSamples, SEXP rPerturbations, 
			SEXP rMaxParents, SEXP rParentSizes, SEXP rMaxComplexity, SEXP rNodeCats, 
			SEXP rParentsPool, SEXP rFixedParentsPool, 
			SEXP rScore, SEXP rWeight, SEXP rMaxIter,
			SEXP rThreads, SEXP rUseCache, SEXP rEcho);
SEXP catnetLoglik(SEXP cnet, SEXP rSamples, SEXP rPerturbations, SEXP rBySample);
SEXP catnetNodeLoglik(SEXP cnet, SEXP rNode, SEXP rSamples, SEXP rPerturbations);
SEXP catnetSetProb(SEXP cnet, SEXP rSamples, SEXP rPerturbations);
SEXP catnetEntropyPairwise(SEXP rSamples, SEXP rPerturbations);
SEXP catnetEntropyOrder(SEXP rSamples, SEXP rPerturbations);
SEXP catnetKLpairwise(SEXP rSamples, SEXP rPerturbations);
SEXP catnetPearsonPairwise(SEXP rSamples, SEXP rPerturbations);
SEXP catnetSamples(SEXP cnet, SEXP rNumSamples, SEXP rPerturbations, SEXP rNaRate);
SEXP catnetReleaseCache();
SEXP catnetSetSeed(SEXP rSeed);
SEXP catnetSoftQuant(SEXP rSamples, SEXP rSampleWeights, SEXP rNumCats, SEXP rLearnset, SEXP rCover, SEXP rMode, SEXP rMaxiter);


