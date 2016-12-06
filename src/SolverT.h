/*
 * TD_Solver.h
 *
 *  Created on: Jun 2, 2015
 *      Author: xiaoyong
 */

#ifndef SOLVERT_H_
#define SOLVERT_H_

#include "petscksp.h"
#include "petscmat.h"
#include "petscsys.h"
#include "ModelManagerT.h"




namespace TD_BEM{


class SolverT
{
public:
	SolverT(ModelManagerT* model);
	virtual ~SolverT();

	void SetTimeControl(double dt, int Numsteps);

	void LocalIndexing();

	virtual void Initialize();

	void SetCurrStep(int curr_step);

	void SetMaxStep(int MaxStep);

	/**
	 * @Function=update the matrices for linear time interpolation
	 */
	virtual void UpdateGHLinear(const double* G_0, const double* G_1, const double* H_0, const double* H_1);



	virtual void Solve()=0;

	/**
	 * @Function=change columns of fG_Compute and fH_Compute (columns corresponding to fixed displacement dofs)
	 */
	void ExchangeColumns();

	/**
	 * @Function=set fLoad
	 */
	void SetLoad();

	void GetNodeLowHigh(int* low, int* high);
	void GetDofLowHigh(int* low, int* high);

	virtual void RecordResult(void)=0;

	virtual void OutPutTecPlot(void);

	virtual void OutPutVTK(void);
    
    virtual void OutPutHMatrix(void);

protected:

	ModelManagerT* fModel;
	PetscErrorCode fIerr;

	/**
	 * time contrls
	 */
	double fBeta;
	double fDt;
	int fNumStep;
	int fCurrStep;

	int fMaxStep;

	/**
	 * local to global indexing
	 */
	MPI_Comm fComm;
	int fSize;
	int fRank;

	int fNodeLow;
	int fNodeHigh;
	int fDofLow;
	int fDofHigh;

	int fLocNumDof;
	int fGblNumDof;

	int* fIndexRow;
	int* fIndexCol;


	/**
	 * global matrices and field variables
	 */
	Vec* fDis_DB;
	Vec* fTrac_DB;
	Mat* fG_DB;
	Mat* fH_DB;
	Mat* fG_DB_bgn;
	Mat* fH_DB_bgn;

	/**
	 * Matrices for computation
	 */
	Mat fG_Compute;
	Mat fH_Compute;
	Vec fLoad;
	Vec fRHS;


	/**
	 * LinearSolver
	 */
	KSP fKsp;

};


} /* namespace */

#endif /* SOLVERT_H_ */
