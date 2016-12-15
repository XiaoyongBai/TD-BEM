/*
 * SolverStepLinearT.h
 *
 *  Created on: Sep 22, 2015
 *      Author: xiaoyong
 */

#ifndef SOLVERSTEPLINEART_H_
#define SOLVERSTEPLINEART_H_

#include "SolverT.h"

namespace TD_BEM {

class SolverStepLinearT : public SolverT{

public:

	SolverStepLinearT(ModelManagerT* model);

	virtual ~SolverStepLinearT();

	virtual void Initialize(void);

	virtual void UpdateGHLinear(const double* G_0, const double* G_1, const double* H_0, const double* H_1);

	virtual void Solve(double a1, double a2);

	/**
	 * @Function=solve with linear time interpolation
	 */
	void Solve_Direct(void);
	void Solve_Averaging(double a1, double a2);


	virtual void RecordResult(void);

};

} /* namespace TD_BEM */

#endif /* SOLVERSTEPLINEART_H_ */
