/*
 * SolverLinearT.h
 *
 *  Created on: Jun 27, 2015
 *      Author: xiaoyong
 */

#ifndef SOLVERLINEART_H_
#define SOLVERLINEART_H_

#include "SolverT.h"

namespace TD_BEM {

class SolverHeavisideLinearT: public SolverT
{

public:

	SolverHeavisideLinearT(ModelManagerT* model);

	~SolverHeavisideLinearT();

	virtual void Initialize(void);

	virtual void UpdateGHLinear(const double* G_0, const double* G_1, const double* H_0, const double* H_1);

	virtual void Solve();

	/**
	 * @Function=solve with linear time interpolation
	 */
	void Solve_Averaging();

	int ActualIndex(int );

	virtual void RecordResult(void);
};

} /* namespace TD_BEM */

#endif /* SOLVERLINEART_H_ */
