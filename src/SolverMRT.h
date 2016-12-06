/*
 * SolverMRT.h
 *
 *  Created on: Oct 25, 2015
 *      Author: xiaoyong
 *
 *
 * Solver for multiregion problems
 */

#ifndef SOLVERMRT_H_
#define SOLVERMRT_H_


#include "petscksp.h"
#include "petscmat.h"
#include "petscsys.h"
#include <string>

#include "ModelManagerT.h"


namespace TD_BEM {

class SolverMRT {
public:
	SolverMRT(ModelManagerT* model_1, ModelManagerT* model_2);

	virtual ~SolverMRT();

	void SetTimeControl(double dt, int Numsteps);

	void LocalIndexing();
	void InterfaceIndexing();

	//virtual void Initialize();

	void SetCurrStep(int curr_step);

	void SetMaxStep(int MaxStep);

	/**
	 * @Function=update the matrices for linear time interpolation
	 */
	virtual void UpdateGHLinear(const double* G_0, const double* G_1, const double* H_0, const double* H_1,
			                    const double* g_0, const double* g_1, const double* h_0, const double* h_1);


	/**
	 * @Function=Assemble G and H of region 1 to the global matrix
	 */
	void Assemble_GH_R1(Mat& G_global, Mat& H_global, const double* G_local, const double* H_local, int insert_or_add);
	void Assemble_GH_R2(Mat& G_global, Mat& H_global, const double* G_local, const double* H_local, int insert_or_add);

	virtual void Solve();

	void Solve_Direct();
	void Solve_Averaging();

	/**
	 * @Function=set fLoad
	 */
	void SetLoad();

	/*****************************************************************
	 * @Function=extract displacements from Global vectors to
	 *           form local displacement vector for region 1 and 2
	 *****************************************************************/
	void ExtractField();

	void GetNodeLowHigh(int* low_1, int* high_1, int* low_2, int* high_2);
	void GetDofLowHigh(int* low_1, int* high_1, int* low_2, int* high_2);

	/*****************************************************************
	 * @Function=Write the results for post processing
	 *****************************************************************/
	virtual void RecordResult(void);

	/*****************************************************************
	 * @Function=Write one model into VTK format
	 *****************************************************************/
	void OutPutVTK(std::string filename, Vec* Dis_DB, ModelManagerT* Model);


protected:

	ModelManagerT* fModel_1;
	ModelManagerT* fModel_2;

	PetscErrorCode fIerr;

	/**
	 * time controls
	 */
	double fBeta;
	double fDt;
	int fNumStep;
	int fCurrStep;

	int fMaxStep;

	/**
	 * MPI indexing
	 */
	MPI_Comm fComm;
	int fSize;
	int fRank;

	/**
	 * dimensions
	 */
	int fNodeLow_1;  //The lowest number of node belonging to this processor, in region 1
	int fNodeHigh_1; //... highest ...

	int fNodeLow_2;
	int fNodeHigh_2;

	int fDofLow_1;
	int fDofLow_2;

	int fDofHigh_1;
	int fDofHigh_2;

	int fNumNodes_1; //total number of nodes in region 1
	int fNumNodes_2;

	int fNumDofs_1; //total number of DOFs in region 1
	int fNumDofs_2;

	int fNumInterfaceNodes; //number of node on interface
	int fNumInterfaceDofs;

	int fNumInteriorNodes_1; //number of interior nodes in region 1
	int fNumInteriorNodes_2;

	int fNumInteriorDofs_1;
	int fNumInteriorDofs_2;

	int fNumRow_Local_1; //Local number of row of local H and G belonging to current processor
	int fNumRow_Local_2;

	int fNumRow_H; //Global number of row of H
	int fNumCol_H; //Global number of column of H

	int fNumRow_G;
	int fNumCol_G;

	/**
	 * Matrices Extraction indexing
	 */
	int* fExtraction_row_1;
	int* fExtraction_row_2;

	int* fExtraction_Hr_col_1;
	int* fExtraction_Hr_col_2;

	int* fExtraction_Hz_col_1;
	int* fExtraction_Hz_col_2;

	int* fExtraction_Gr_col_1;
	int* fExtraction_Gr_col_2;

	int* fExtraction_Gz_col_1;
	int* fExtraction_Gz_col_2;

	/**
	 * Matrices Assemble indexing
	 */
	int* fAssemble_row_1;
	int* fAssemble_row_2;

	int* fAssemble_Hr_col_1;
	int* fAssemble_Hr_col_2;

	int* fAssemble_Hz_col_1;
	int* fAssemble_Hz_col_2;

	int* fAssemble_Gr_col_1;
	int* fAssemble_Gr_col_2;

	int* fAssemble_Gz_col_1;
	int* fAssemble_Gz_col_2;


	/*****************************************
	 * Indices of extracting Field variables
	 ******************************************/
	int* fExtraction_Dis_1; //Global IDs of interface nodes in region 1
	int* fExtraction_Dis_2;

	int* fSet_Dis_1;
	int* fSet_Dis_2;


	/**
	 * global matrices and field variables
	 */
	Vec* fDis_DB_1; //displacements for region 1
	Vec* fTrac_DB_1; //Traction for region 1

	Vec* fDis_DB_2;
	Vec* fTrac_DB_2;

	Vec* fDis_DB_Global; //Displacements combined for computation, composed of Ur1, Ur2, Uz, Tz
	Vec* fTrac_DB_Global; //Traction combined for computation, compoised of Tr1, Tr2

	Mat* fG_DB; //Combined G matrices
	Mat* fH_DB; //Combined H matrices
	Mat* fG_DB_bgn;
	Mat* fH_DB_bgn;

	/**
	 * Matrices for computation
	 */
	Mat fG_Compute;
	Mat fH_Compute;
	Vec fLoad; //load vector for current step
	Vec fRHS; //right hand side of linear equation system


	/**
	 * LinearSolver
	 */
	KSP fKsp;
};

} /*namespace TD_BEM */

#endif /* SOLVERMRT_H_ */
