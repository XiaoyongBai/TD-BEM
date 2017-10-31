/*
 * SolverLinearT.cpp
 *
 *  Created on: Jun 27, 2015
 *      Author: xiaoyong
 */

#include "SolverHeavisideLinearT.h"
#include "MathOperationT.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

using namespace TD_BEM;
using namespace GreenFunction;

SolverHeavisideLinearT::SolverHeavisideLinearT(ModelManagerT* model):
	SolverT(model)
{
	// TODO Auto-generated constructor stub

}

SolverHeavisideLinearT::~SolverHeavisideLinearT() {
	// TODO Auto-generated destructor stub
}

void SolverHeavisideLinearT::Initialize()
{
	SolverT::Initialize();

	int UBC_num;
	int* UBC_DOFs;
	double* UBC_Values;
	int FBC_num;
	int* FBC_DOFs;
	double* FBC_Values;

	fModel->BoundaryCondition_Local(UBC_num, &UBC_DOFs, &UBC_Values, FBC_num, &FBC_DOFs, &FBC_Values);

	//MathOperationT::PrintVector(UBC_num, UBC_DOFs, "UBC_DOFs");
	//MathOperationT::PrintVector(UBC_num, UBC_Values, "UBC_Values");
	//MathOperationT::PrintVector(FBC_num, FBC_DOFs, "FBC_DOFs");
	//MathOperationT::PrintVector(FBC_num, FBC_Values, "FBC_Values");

	fIerr=VecSetValues(fTrac_DB[0], FBC_num, FBC_DOFs, FBC_Values, INSERT_VALUES);

	fIerr=VecAssemblyBegin(fTrac_DB[0]);
	fIerr=VecAssemblyEnd(fTrac_DB[0]);

	//VecView(fTrac_DB[0], PETSC_VIEWER_STDOUT_WORLD);
}

void SolverHeavisideLinearT::UpdateGHLinear(const double* G_0, const double* G_1, const double* H_0, const double* H_1)
{
	if (fCurrStep>fMaxStep+2)
		return;

	SolverT::UpdateGHLinear(G_0, G_1, H_0, H_1);

	//MathOperationT::PrintMatrix(fLocNumDof, fGblNumDof, G_0, "Global G_0");
	//MathOperationT::PrintMatrix(fLocNumDof, fGblNumDof, G_1, "Global G_1");
	//MathOperationT::PrintMatrix(fLocNumDof, fGblNumDof, H_0, "Global H_0");
	//MathOperationT::PrintMatrix(fLocNumDof, fGblNumDof, H_1, "Global H_1");

	/*************
	 * Arrange the matrices
	 *************/
	if(fCurrStep==1)
	{
		fIerr=MatSetValues(fG_DB[0], fLocNumDof, fIndexRow, fGblNumDof, fIndexCol, G_1, INSERT_VALUES);
		fIerr=MatAssemblyBegin(fG_DB[0], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fG_DB[0], MAT_FINAL_ASSEMBLY);

		fIerr=MatSetValues(fH_DB[0], fLocNumDof, fIndexRow, fGblNumDof, fIndexCol,  H_1, INSERT_VALUES);
		fIerr=MatAssemblyBegin(fH_DB[0], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fH_DB[0], MAT_FINAL_ASSEMBLY);

	}
	else
	{
		fIerr=MatSetValues(fG_DB[fCurrStep-1], fLocNumDof, fIndexRow, fGblNumDof, fIndexCol,  G_1, ADD_VALUES);
		fIerr=MatAssemblyBegin(fG_DB[fCurrStep-1], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fG_DB[fCurrStep-1], MAT_FINAL_ASSEMBLY);

		fIerr=MatSetValues(fH_DB[fCurrStep-1], fLocNumDof, fIndexRow, fGblNumDof, fIndexCol, H_1, ADD_VALUES);
		fIerr=MatAssemblyBegin(fH_DB[fCurrStep-1], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fH_DB[fCurrStep-1], MAT_FINAL_ASSEMBLY);
	}


	fIerr=MatSetValues(fG_DB[fCurrStep], fLocNumDof, fIndexRow, fGblNumDof, fIndexCol, G_0, INSERT_VALUES);
	fIerr=MatAssemblyBegin(fG_DB[fCurrStep], MAT_FINAL_ASSEMBLY);
	fIerr=MatAssemblyEnd(fG_DB[fCurrStep], MAT_FINAL_ASSEMBLY);

	fIerr=MatSetValues(fH_DB[fCurrStep], fLocNumDof, fIndexRow, fGblNumDof, fIndexCol, H_0, INSERT_VALUES);
	fIerr=MatAssemblyBegin(fH_DB[fCurrStep], MAT_FINAL_ASSEMBLY);
	fIerr=MatAssemblyEnd(fH_DB[fCurrStep], MAT_FINAL_ASSEMBLY);


	//set fG_DB_bgn
	if(fCurrStep==1)
	{
		fIerr=MatSetValues(fG_DB_bgn[0], fLocNumDof, fIndexRow, fGblNumDof, fIndexCol, G_0, INSERT_VALUES);
		fIerr=MatAssemblyBegin(fG_DB_bgn[0], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fG_DB_bgn[0], MAT_FINAL_ASSEMBLY);

		/********************************************************************
		 ********************************************************************/

		fIerr=MatSetValues(fH_DB_bgn[0], fLocNumDof, fIndexRow, fGblNumDof, fIndexCol, H_0, INSERT_VALUES);
		fIerr=MatAssemblyBegin(fH_DB_bgn[0], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fH_DB_bgn[0], MAT_FINAL_ASSEMBLY);
	}
	else if(fCurrStep==2)
	{
		fIerr=MatCopy(fG_DB_bgn[0], fG_DB_bgn[1], DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAssemblyBegin(fG_DB_bgn[1], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fG_DB_bgn[1], MAT_FINAL_ASSEMBLY);

		fIerr=MatSetValues(fG_DB_bgn[0], fLocNumDof, fIndexRow, fGblNumDof, fIndexCol, G_0, INSERT_VALUES);
		fIerr=MatAssemblyBegin(fG_DB_bgn[0], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fG_DB_bgn[0], MAT_FINAL_ASSEMBLY);

		/********************************************************************
		 ********************************************************************/
		fIerr=MatCopy(fH_DB_bgn[0], fH_DB_bgn[1], DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAssemblyBegin(fH_DB_bgn[1], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fH_DB_bgn[1], MAT_FINAL_ASSEMBLY);

		fIerr=MatSetValues(fH_DB_bgn[0], fLocNumDof, fIndexRow, fGblNumDof, fIndexCol, H_0, INSERT_VALUES);
		fIerr=MatAssemblyBegin(fH_DB_bgn[0], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fH_DB_bgn[0], MAT_FINAL_ASSEMBLY);
	}
	else
	{
		fIerr=MatCopy(fG_DB_bgn[1], fG_DB_bgn[2], DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAssemblyBegin(fG_DB_bgn[2], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fG_DB_bgn[2], MAT_FINAL_ASSEMBLY);

		fIerr=MatCopy(fG_DB_bgn[0], fG_DB_bgn[1], DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAssemblyBegin(fG_DB_bgn[1], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fG_DB_bgn[1], MAT_FINAL_ASSEMBLY);

		fIerr=MatSetValues(fG_DB_bgn[0], fLocNumDof, fIndexRow, fGblNumDof, fIndexCol, G_0, INSERT_VALUES);
		fIerr=MatAssemblyBegin(fG_DB_bgn[0], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fG_DB_bgn[0], MAT_FINAL_ASSEMBLY);

		/********************************************************************
		 ********************************************************************/

		fIerr=MatCopy(fH_DB_bgn[1], fH_DB_bgn[2], DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAssemblyBegin(fH_DB_bgn[2], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fH_DB_bgn[2], MAT_FINAL_ASSEMBLY);

		fIerr=MatCopy(fH_DB_bgn[0], fH_DB_bgn[1], DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAssemblyBegin(fH_DB_bgn[1], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fH_DB_bgn[1], MAT_FINAL_ASSEMBLY);

		fIerr=MatSetValues(fH_DB_bgn[0], fLocNumDof, fIndexRow, fGblNumDof, fIndexCol, H_0, INSERT_VALUES);
		fIerr=MatAssemblyBegin(fH_DB_bgn[0], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fH_DB_bgn[0], MAT_FINAL_ASSEMBLY);
	}

}


void SolverHeavisideLinearT::Solve()
{

	Solve_Averaging();
}


void SolverHeavisideLinearT::Solve_Averaging()
{
	int UBC_num;
	int* UBC_DOFs;
	double* UBC_Values;
	int FBC_num;
	int* FBC_DOFs;
	double* FBC_Values;

	fModel->BoundaryCondition_Local(UBC_num, &UBC_DOFs, &UBC_Values, FBC_num, &FBC_DOFs, &FBC_Values);

	//if(fCurrStep>30)
	//{
		//fIerr=VecZeroEntries(fLoad);

		//MathOperationT::VecSet(FBC_num, FBC_Values, 0);
	//}

	//double amp=sin(10*M_PI*(fCurrStep-1)*fDt); // for sinusoidal load
	double amp=1; //for constant load

	double* FBC_Values_temp=new double[FBC_num];
	MathOperationT::MemCopy(FBC_num, FBC_Values, FBC_Values_temp);
	MathOperationT::VecScale(FBC_num, FBC_Values_temp, amp);

	//MathOperationT::PrintVector(FBC_num, FBC_Values, "FBC");
	//MathOperationT::PrintVector(FBC_num, FBC_Values_temp, "FBC_temp");

	//MathOperationT::PrintVector(UBC_num, UBC_Values, "solver UBC_values");
	//MathOperationT::PrintVector(FBC_num, FBC_Values, "solver FBC_values");

	if(fCurrStep==1)
	{
		//do nothing
		return;
	}
	else if(fCurrStep==2)
	{
		/********
		 * Form LHS
		 ********/
		int m,n;
		MatGetSize(fG_Compute, &m, &n);

		fIerr=MatCopy(fG_DB[0], fG_Compute, DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAssemblyBegin(fG_Compute, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fG_Compute, MAT_FINAL_ASSEMBLY);

		fIerr=MatCopy(fH_DB[0], fH_Compute, DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAssemblyBegin(fH_Compute, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fH_Compute, MAT_FINAL_ASSEMBLY);

		ExchangeColumns();

		/********
		 * Form RHS
		 ********/
		fIerr=MatMult(fG_Compute, fLoad, fRHS);
		fIerr=MatMultAdd(fG_DB_bgn[1], fTrac_DB[0], fRHS, fRHS);
	}
	else
	{
		/********
		 * Form LHS
		 ********/
		fIerr=MatZeroEntries(fG_Compute);
		fIerr=MatZeroEntries(fH_Compute);

		fIerr=MatAXPY(fG_Compute, 4, fG_DB[0],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAXPY(fG_Compute, 1, fG_DB[1],DIFFERENT_NONZERO_PATTERN);

		fIerr=MatAXPY(fH_Compute, 4, fH_DB[0],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAXPY(fH_Compute, 1, fH_DB[1],DIFFERENT_NONZERO_PATTERN);

		ExchangeColumns();

		//MatView(fG_Compute, PETSC_VIEWER_STDOUT_WORLD);

		/********
		 * Form RHS
		 ********/
		fIerr=MatMult(fG_Compute, fLoad, fRHS);

		//MatView(fG_Compute, PETSC_VIEWER_STDOUT_WORLD);
		//VecView(fLoad, PETSC_VIEWER_STDOUT_WORLD);


		//**add single terms
		Mat matrix_temp;
		fIerr=MatConvert(fG_DB[2], MATSAME, MAT_INITIAL_MATRIX, &matrix_temp);
		fIerr=MatAXPY(matrix_temp, 2, fG_DB[1],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatMultAdd(matrix_temp, fTrac_DB[fCurrStep-2],fRHS, fRHS);

		//MatView(matrix_temp, PETSC_VIEWER_STDOUT_WORLD);
		//VecView(fRHS, PETSC_VIEWER_STDOUT_WORLD);

		fIerr=MatConvert(fH_DB[2], MATSAME, MAT_INITIAL_MATRIX, &matrix_temp);
		fIerr=MatAXPY(matrix_temp, 2, fH_DB[1],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatScale(matrix_temp, -1);
		fIerr=MatMultAdd(matrix_temp, fDis_DB[fCurrStep-2],fRHS, fRHS);


		//*** add terms by loop
		for(int f_i=1; f_i<=fCurrStep-3; f_i++)
		{
			fIerr=MatZeroEntries(matrix_temp);
			fIerr=MatAXPY(matrix_temp, 1, fG_DB[ActualIndex(fCurrStep-f_i)],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatAXPY(matrix_temp, 2, fG_DB[ActualIndex(fCurrStep-f_i-1)],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatAXPY(matrix_temp, 1, fG_DB[ActualIndex(fCurrStep-f_i-2)],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatMultAdd(matrix_temp, fTrac_DB[f_i],fRHS, fRHS);

			fIerr=MatZeroEntries(matrix_temp);
			fIerr=MatAXPY(matrix_temp, 1, fH_DB[ActualIndex(fCurrStep-f_i)],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatAXPY(matrix_temp, 2, fH_DB[ActualIndex(fCurrStep-f_i-1)],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatAXPY(matrix_temp, 1, fH_DB[ActualIndex(fCurrStep-f_i-2)],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatScale(matrix_temp, -1);
			fIerr=MatMultAdd(matrix_temp, fDis_DB[f_i],fRHS, fRHS);
		}

		//** special treatment to the first step
		fIerr=MatZeroEntries(matrix_temp);
		fIerr=MatAXPY(matrix_temp, 1, fG_DB_bgn[0],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAXPY(matrix_temp, 2, fG_DB_bgn[1],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAXPY(matrix_temp, 1, fG_DB_bgn[2],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatMultAdd(matrix_temp, fTrac_DB[0],fRHS, fRHS);

		fIerr=MatZeroEntries(matrix_temp);
		fIerr=MatAXPY(matrix_temp, 1, fH_DB_bgn[0],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAXPY(matrix_temp, 2, fH_DB_bgn[1],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAXPY(matrix_temp, 1, fH_DB_bgn[2],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatScale(matrix_temp, -1);
		fIerr=MatMultAdd(matrix_temp, fDis_DB[0],fRHS, fRHS);

		fIerr=MatDestroy(&matrix_temp);
	}

	cout << "current step " <<fCurrStep << endl;

	Vec result_rough;
	fIerr=VecCreateMPI(fComm, PETSC_DECIDE, fGblNumDof, &result_rough);

	//fIerr=MatView(fH_Compute, PETSC_VIEWER_STDOUT_WORLD);
	//cout << "fLoad is"<<endl;
    	//fIerr=VecView(fLoad, PETSC_VIEWER_STDOUT_WORLD);

	//cout<<"fRHS is " << endl;
    	//fIerr=VecView(fRHS, PETSC_VIEWER_STDOUT_WORLD);

	fIerr=KSPSetOperators(fKsp, fH_Compute, fH_Compute);
	fIerr=KSPSolve(fKsp, fRHS, result_rough);
    	VecAssemblyBegin(result_rough);
    	VecAssemblyEnd(result_rough);

	//cout << "result _rough is \n";
    	//VecView(result_rough, PETSC_VIEWER_STDOUT_WORLD);

    	int solve_step=fCurrStep-1;

    	fIerr=VecCopy(result_rough, fDis_DB[solve_step]);
	fIerr=VecSetValues(fDis_DB[solve_step], UBC_num, UBC_DOFs, UBC_Values, INSERT_VALUES);
	fIerr=VecAssemblyBegin(fDis_DB[solve_step]);
	fIerr=VecAssemblyEnd(fDis_DB[solve_step]);

	//if(fCurrStep==10)
	// VecView(fDis_DB[solve_step], PETSC_VIEWER_STDOUT_WORLD);

	fIerr=VecSetValues(fTrac_DB[solve_step], FBC_num, FBC_DOFs, FBC_Values_temp, INSERT_VALUES);
	int low, high;
	fIerr=VecGetOwnershipRange(fTrac_DB[solve_step], &low, &high);
	for(int u_i=0; u_i<UBC_num; u_i++)
	{
		int u_dof=UBC_DOFs[u_i];

		if(u_dof>=fDofLow && u_dof<fDofHigh)
		{
			double t_temp;
			fIerr=VecGetValues(result_rough, 1, &u_dof, &t_temp);
			fIerr=VecSetValues(fTrac_DB[solve_step], 1, &u_dof, &t_temp, INSERT_VALUES);
		}
	}
	fIerr=VecAssemblyBegin(fTrac_DB[solve_step]);
	fIerr=VecAssemblyEnd(fTrac_DB[solve_step]);

	fIerr=VecDestroy(&result_rough);

	delete[] FBC_Values_temp;

	//fIerr=VecView(fTrac_DB[solve_step], PETSC_VIEWER_STDOUT_WORLD);
	//fIerr=MatMult(fH_Compute, result_rough, fRHS);
}


int SolverHeavisideLinearT::ActualIndex(int nominal)
{
	if(nominal<=fMaxStep)
	{
		return nominal;
	}
	else
	{
		return fMaxStep+1;
	}
}


void SolverHeavisideLinearT::RecordResult()
{
	if(fRank==0)
	{
		ofstream myfile;
		myfile.open("displacement.txt");

		double* temp;
		for(int si=0; si<fCurrStep; si++)
		{
			fIerr=VecGetArray(fDis_DB[si], &temp);
			myfile<<setprecision(6);
			//myfile<< setw(10) << fDt*si << setw(15) << temp[261*3-1] << setw(15) << temp[22*3-1] << setw(15) << temp[442*3-1]  << setw(15) << temp[499*3-1]  << endl;
			myfile<<setw(10) << fDt*si << setw(15) << temp[2] << endl;

			fIerr=VecRestoreArray(fDis_DB[si], &temp);
		}

		myfile.close();
	}

}
