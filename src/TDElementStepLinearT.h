/*
 * TDElementStepLinearT.h
 *
 *  Created on: Sep 21, 2015
 *      Author: xiaoyong
 */

#ifndef TDELEMENTSTEPLINEART_H_
#define TDELEMENTSTEPLINEART_H_

#include "TDElementT.h"

namespace TD_BEM {

class TDElementStepLinearT: public TDElementT{
public:
	TDElementStepLinearT(int type);
	virtual ~TDElementStepLinearT();

	virtual void FormGH_Regular();

	virtual void FormGH_Singular(int singular);

	/* form G and H matrix for bi-linear element */
	void FormGH_Singular_1(int singular, ShapeFunctionT* Shape, double* GE_1, double* GE_2, double* HE_1, double* HE_2, double* C_1, double* C_2);

	void SubElement(int singular);

	void AssembleSubElement(const int* ids, const double* GE_1, const double* GE_2, const double* HE_1, const double* HE_2, const double* C_1, const double* C_2);




};

} /* namespace TD_BEM */

#endif /* TDELEMENTSTEPLINEART_H_ */
