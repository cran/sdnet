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
 * rcatnet.h
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

#ifndef RCATNET_H_
#define RCATNET_H_

#include "catnetd.h"
#include "catnetp.h"

#define MAX_NODE_NAME	16

class RCatnet : public CATNETD<double> {
public:
	RCatnet(){
	}

	RCatnet(SEXP cnet);

	SEXP genRcatnet(const char * objectName);
	SEXP genRcatnetFromDagEvaluate(SEXP rDagEval, SEXP rDagIndex);

	SEXP genProbList(int node, int paridx, int *pcats);

	RCatnet& operator =(CATNETD<double> &cnet) {
		CATNET<int, double>::init(
 			cnet.numNodes(), cnet.maxParents(), cnet.maxCategories(),
			cnet.nodeNames(), cnet.numParents(), cnet.parents(),
			cnet.numCategories(), cnet.probLists());
		return *this;
	}

	SEXP genSamples(SEXP rNumSamples, SEXP rPerturbations, SEXP rNaRate);

};

class RCatnetP : public CATNETP<double> {
public:
	RCatnetP() {
	}

	RCatnetP(SEXP cnet);

	SEXP genRcatnet(const char * objectName);
	SEXP genRcatnetFromDagEvaluate(SEXP rDagEval, SEXP rDagIndex);

	SEXP genProbList(int node, int paridx, int *pcats);

	RCatnetP& operator =(CATNETP<double> &cnet) {
		CATNET<double, double>::init(
 			cnet.numNodes(), cnet.maxParents(), cnet.maxCategories(),
			cnet.nodeNames(), cnet.numParents(), cnet.parents(),
			cnet.numCategories(), cnet.probLists());
		return *this;
	}

	SEXP genSamples(SEXP rNumSamples, SEXP rPerturbations, SEXP rNaRate);
};

extern "C" {

char *gen_prob_string(int node, SEXP parlist, int paridx, SEXP catlist, SEXP problist, char *str);
SEXP prob_string(SEXP rnodes, SEXP rparents, SEXP rcatlist, SEXP rproblist);
void gen_prob_vector(int node, SEXP parlist, int paridx, SEXP catlist, SEXP problist, double *&pvec, int &nvec);
SEXP prob_vector(SEXP rnodes, SEXP rparents, SEXP rcatlist, SEXP rproblist); 

} // extern "C"


#endif /* RCATNET_H_ */
