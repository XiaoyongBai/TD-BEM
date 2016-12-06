/*
 * TDElementLinearT.h
 *
 *  Created on: Jun 27, 2015
 *      Author: xiaoyong
 */

#ifndef TDELEMENTRAMPHEAVILINEART_H_
#define TDELEMENTRAMPHEAVILINEART_H_

#include "TDElementT.h"

namespace TD_BEM {

class TDElementRampHeaviLinearT: public TDElementT {

public:

	TDElementRampHeaviLinearT(int type);

	virtual ~TDElementRampHeaviLinearT();

	virtual void FormGH_Regular();

	virtual void FormGH_Singular(int singular);

	/* form G and H matrix for bi-linear element */
	void FormGH_Singular_1(int singular, ShapeFunctionT* Shape, double* GE_1, double* GE_2, double* HE_1, double* HE_2);

	void SubElement(int singular);

	void AssembleSubElement(const int* ids, const double* GE_1, const double* GE_2, const double* HE_1, const double* HE_2);

};

} /* namespace TD_BEM */

#endif /* TDELEMENTLINEART_H_ */
