/*
 * SolverStepLinearT.cpp
 *
 *  Created on: Sep 22, 2015
 *      Author: xiaoyong
 */

#include "SolverStepLinearT.h"
#include "MathOperationT.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

using namespace TD_BEM;
using namespace GreenFunction;

SolverStepLinearT::SolverStepLinearT(ModelManagerT* model):
	SolverT(model)
{
	// TODO Auto-generated constructor stub
}

SolverStepLinearT::~SolverStepLinearT() {
	// TODO Auto-generated destructor stub
}

void SolverStepLinearT::Initialize()
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

void SolverStepLinearT::UpdateGHLinear(const double* G_0, const double* G_1, const double* H_0, const double* H_1)
{

    SolverT::UpdateGHLinear(G_0, G_1, H_0, H_1);

	/*************
	 * Arrange the matrices
	 *************/
	if (fCurrStep<=fMaxStep)
	{
        //Print the matrices
        string G1_file, G2_file, H1_file, H2_file;
        
        ostringstream convert_G1, convert_G2, convert_H1, convert_H2; // stream used for the conversion
        
        int rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        
        convert_G1 << "G1_step" <<fCurrStep-1<<"_p_"<<rank<<".txt";
        G1_file=convert_G1.str();
        convert_G2 << "G2_step" <<fCurrStep-1<<"_p_"<<rank<<".txt";
        G2_file=convert_G2.str();
        convert_H1 << "H1_step" <<fCurrStep-1<<"_p_"<<rank<<".txt";
        H1_file=convert_H1.str();
        convert_H2 << "H2_step" <<fCurrStep-1<<"_p_"<<rank<<".txt";
        H2_file=convert_H2.str();
        
        ofstream G1_write, G2_write, H1_write, H2_write;
        
        G1_write.open(G1_file.c_str(), ios_base::out);
        for (int i=0; i<fLocNumDof; i++) {
            for (int j=0; j<fGblNumDof; j++) {
                G1_write<< setw(15) << setprecision(5) << std::scientific<<G_1[i*fGblNumDof+j];
            }
            G1_write<<endl;
        }
        G1_write.close();
        
        G2_write.open(G2_file.c_str(), ios_base::out);
        for (int i=0; i<fLocNumDof; i++) {
            for (int j=0; j<fGblNumDof; j++) {
                G2_write<< setw(15) << setprecision(5) << std::scientific<<G_0[i*fGblNumDof+j];
            }
            G2_write<<endl;
        }
        G2_write.close();
        
        H1_write.open(H1_file.c_str(), ios_base::out);
        for (int i=0; i<fLocNumDof; i++) {
            for (int j=0; j<fGblNumDof; j++) {
                H1_write<< setw(15) << setprecision(5) << std::scientific<<H_1[i*fGblNumDof+j];
            }
            G1_write<<endl;
        }
        G1_write.close();
        
        H2_write.open(H2_file.c_str(), ios_base::out);
        for (int i=0; i<fLocNumDof; i++) {
            for (int j=0; j<fGblNumDof; j++) {
                H2_write<< setw(15) << setprecision(5) << std::scientific<<H_0[i*fGblNumDof+j];
            }
            H2_write<<endl;
        }
        H2_write.close();
        

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
	}



	//set fG_DB_bgn
	if (fCurrStep>fMaxStep+3)
		return;

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

		if (fCurrStep>fMaxStep)
		{
			fIerr=MatScale(fH_DB_bgn[0], 0);

			fIerr=MatScale(fG_DB_bgn[0], 0);
		}
	}

}


void SolverStepLinearT::Solve()
{
	Solve_Averaging();
	//Solve_Direct();
}


void SolverStepLinearT::Solve_Direct()
{
	int UBC_num;
	int* UBC_DOFs;
	double* UBC_Values;
	int FBC_num;
	int* FBC_DOFs;
	double* FBC_Values;

	fModel->BoundaryCondition_Local(UBC_num, &UBC_DOFs, &UBC_Values, FBC_num, &FBC_DOFs, &FBC_Values);

	double* FBC_Values_temp=new double[FBC_num];
	MathOperationT::MemCopy(FBC_num, FBC_Values, FBC_Values_temp);

	/********
	 * Form LHS
	 ********/
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

	//*** add terms by loop
	Mat matrix_temp;
	fIerr=MatConvert(fG_DB[0], MATSAME, MAT_INITIAL_MATRIX, &matrix_temp);

	for(int f_i=1; f_i<=min(fCurrStep, fMaxStep); f_i++)
	{
		fIerr=MatZeroEntries(matrix_temp);
		fIerr=MatAXPY(matrix_temp, 1, fG_DB[f_i],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatMultAdd(matrix_temp, fTrac_DB[fCurrStep-f_i],fRHS, fRHS);

		fIerr=MatZeroEntries(matrix_temp);
		fIerr=MatAXPY(matrix_temp, 1, fH_DB[f_i],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatScale(matrix_temp, -1);
		fIerr=MatMultAdd(matrix_temp, fDis_DB[fCurrStep-f_i],fRHS, fRHS);
	}


	cout << "current step " << fCurrStep << endl;


	/********
	 * reset LHS and RHS due to the effect of interface
	 ********/
	/*int num_interface=fModel->GetNumInterface();
	const int* interface=fModel->GetInterface();

	int* Interface_dof=new int[3*num_interface];

	for (int Ii=0; Ii<3*num_interface; Ii++)
	{
		for(int dj=0; dj<3; dj++)
			Interface_dof[3*Ii+dj]=3*interface[Ii]+dj;
	}*/

	//MathOperationT::PrintVector(3*num_interface, Interface_dof, "fInterfaceDofs");



	Vec result_rough;
	fIerr=VecCreateMPI(fComm, PETSC_DECIDE, fGblNumDof, &result_rough);

	//cout <<"fH_Compute is \n";
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


    fIerr=VecCopy(result_rough, fDis_DB[fCurrStep]);
	fIerr=VecSetValues(fDis_DB[fCurrStep], UBC_num, UBC_DOFs, UBC_Values, INSERT_VALUES);
	fIerr=VecAssemblyBegin(fDis_DB[fCurrStep]);
	fIerr=VecAssemblyEnd(fDis_DB[fCurrStep]);



	//if(fCurrStep==10)
	// VecView(fDis_DB[solve_step], PETSC_VIEWER_STDOUT_WORLD);

	fIerr=VecSetValues(fTrac_DB[fCurrStep], FBC_num, FBC_DOFs, FBC_Values_temp, INSERT_VALUES);
	int low, high;
	fIerr=VecGetOwnershipRange(fTrac_DB[fCurrStep], &low, &high);
	for(int u_i=0; u_i<UBC_num; u_i++)
	{
		int u_dof=UBC_DOFs[u_i];

		if(u_dof>=fDofLow && u_dof<fDofHigh)
		{
			double t_temp;
			fIerr=VecGetValues(result_rough, 1, &u_dof, &t_temp);
			fIerr=VecSetValues(fTrac_DB[fCurrStep], 1, &u_dof, &t_temp, INSERT_VALUES);
		}
	}

	fIerr=VecAssemblyBegin(fTrac_DB[fCurrStep]);
	fIerr=VecAssemblyEnd(fTrac_DB[fCurrStep]);

	fIerr=VecDestroy(&result_rough);

	delete[] FBC_Values_temp;

    //fIerr=VecView(fDis_DB[fCurrStep], PETSC_VIEWER_STDOUT_WORLD);
	//fIerr=MatMult(fH_Compute, result_rough, fRHS);
}





void SolverStepLinearT::Solve_Averaging()
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

	double p0=5.46366e9;
	double curr_time=(fCurrStep-1)*fDt;
	double amp=p0*curr_time*exp(-2000*curr_time);

	amp=1;

	double* FBC_Values_temp=new double[FBC_num];
	MathOperationT::MemCopy(FBC_num, FBC_Values, FBC_Values_temp);
	MathOperationT::VecScale(FBC_num, FBC_Values_temp, amp);

	//MathOperationT::PrintVector(FBC_num, FBC_Values, "FBC");
	//MathOperationT::PrintVector(FBC_num, FBC_Values_temp, "FBC_temp");

	//MathOperationT::PrintVector(UBC_num, UBC_Values, "solver UBC_values");
	//MathOperationT::PrintVector(FBC_num, FBC_Values, "solver FBC_values");

    int solve_step=fCurrStep-1;


	if(fCurrStep==1)
	{
		//do nothing
		delete[] FBC_Values_temp;
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
		fIerr=MatMultAdd(matrix_temp, fTrac_DB[solve_step-1],fRHS, fRHS);

		fIerr=MatDestroy(&matrix_temp);

		//MatView(matrix_temp, PETSC_VIEWER_STDOUT_WORLD);
		//VecView(fRHS, PETSC_VIEWER_STDOUT_WORLD);

		fIerr=MatConvert(fH_DB[2], MATSAME, MAT_INITIAL_MATRIX, &matrix_temp);
		fIerr=MatAXPY(matrix_temp, 2, fH_DB[1],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatScale(matrix_temp, -1);
		fIerr=MatMultAdd(matrix_temp, fDis_DB[solve_step-1],fRHS, fRHS);


		//*** add terms by loop
		for(int f_i=3; f_i<=min(solve_step, fMaxStep); f_i++)
		{
			fIerr=MatZeroEntries(matrix_temp);
			fIerr=MatAXPY(matrix_temp, 1, fG_DB[f_i],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatMultAdd(matrix_temp, fTrac_DB[solve_step+1-f_i],fRHS, fRHS);

			fIerr=MatZeroEntries(matrix_temp);
			fIerr=MatAXPY(matrix_temp, 1, fH_DB[f_i],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatScale(matrix_temp, -1);
			fIerr=MatMultAdd(matrix_temp, fDis_DB[solve_step+1-f_i],fRHS, fRHS);
		}

		for(int f_i=2; f_i<=min(solve_step-1, fMaxStep); f_i++)
		{
			fIerr=MatZeroEntries(matrix_temp);
			fIerr=MatAXPY(matrix_temp, 2, fG_DB[f_i],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatMultAdd(matrix_temp, fTrac_DB[solve_step-f_i],fRHS, fRHS);

			fIerr=MatZeroEntries(matrix_temp);
			fIerr=MatAXPY(matrix_temp, 2, fH_DB[f_i],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatScale(matrix_temp, -1);
			fIerr=MatMultAdd(matrix_temp, fDis_DB[solve_step-f_i],fRHS, fRHS);
		}

		for(int f_i=1; f_i<=min(solve_step-2, fMaxStep); f_i++)
		{
			fIerr=MatZeroEntries(matrix_temp);
			fIerr=MatAXPY(matrix_temp, 1, fG_DB[f_i],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatMultAdd(matrix_temp, fTrac_DB[solve_step-1-f_i],fRHS, fRHS);

			fIerr=MatZeroEntries(matrix_temp);
			fIerr=MatAXPY(matrix_temp, 1, fH_DB[f_i],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatScale(matrix_temp, -1);
			fIerr=MatMultAdd(matrix_temp, fDis_DB[solve_step-1-f_i],fRHS, fRHS);
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

	cout << "current step " << solve_step << endl;

	Vec result_rough;
	fIerr=VecCreateMPI(fComm, PETSC_DECIDE, fGblNumDof, &result_rough);

	//cout <<"fH_Compute is \n";
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

	//fIerr=VecView(fDis_DB[solve_step], PETSC_VIEWER_STDOUT_WORLD);
	//fIerr=MatMult(fH_Compute, result_rough, fRHS);
}




void SolverStepLinearT::RecordResult()
{
	ofstream myfile;

	Vec dis_collection;
	Vec trac_collection;
		
	VecScatter dis_ctx, trac_ctx;

	VecScatterCreateToAll(fDis_DB[fCurrStep-1],&dis_ctx,&dis_collection);
	VecScatterBegin(dis_ctx,fDis_DB[fCurrStep-1],dis_collection,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(dis_ctx,fDis_DB[fCurrStep-1],dis_collection,INSERT_VALUES,SCATTER_FORWARD);
	
	VecScatterCreateToAll(fTrac_DB[fCurrStep-1], &trac_ctx, &trac_collection);
	VecScatterBegin(trac_ctx, fTrac_DB[fCurrStep-1], trac_collection, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(trac_ctx, fTrac_DB[fCurrStep-1], trac_collection, INSERT_VALUES, SCATTER_FORWARD);	

	if(fRank==0)
	{
		if (fCurrStep==1)
			myfile.open("displacement.txt", fstream::out|fstream::trunc);
		else
			myfile.open("displacement.txt", fstream::out|fstream::app);


		double* temp_dis,* temp_trac;

		fIerr=VecGetArray(dis_collection, &temp_dis);
		fIerr=VecGetArray(trac_collection, &temp_trac);
	
		myfile<<setprecision(6);
		//myfile<< setw(10) << fDt*si << setw(15) << temp[261*3-1] << setw(15) << temp[22*3-1] << setw(15) << temp[442*3-1]  << setw(15) << temp[499*3-1]  << endl;
		//myfile<<setw(10) << fDt*(fCurrStep-1) << setw(15) << temp[107] << endl;
		myfile<<setw(10) << fDt*(fCurrStep-1) << setw(15) << temp_dis[77] << setw(15) << temp_trac[2]<< endl;

		/**for buried cavity Jiang **/
		//myfile<< setw(10) << fDt*(fCurrStep-1) << setw(15) << temp[261*3-1] << setw(15) << temp[413*3-1] << setw(15) << temp[442*3-1]  << setw(15) << temp[591*3-1]  << endl;

		/***for buried cavity jiang_2**/
		//myfile<< setw(10) << fDt*(fCurrStep-1) << setw(15) << temp[921*3-1] << setw(15) << temp[1682*3-1] << setw(15) << temp[1831*3-1]<< endl;

		/***for buried cavity jiang_4**/
		//myfile<< setw(10) << fDt*(fCurrStep-1) << setw(15) << temp[3441*3-1] << setw(15) << temp[6562*3-1] << setw(15) << temp[6711*3-1]<< endl;


		fIerr=VecRestoreArray(dis_collection, &temp_dis);
		fIerr=VecRestoreArray(trac_collection, &temp_trac);
		myfile.close();
	}

	VecScatterDestroy(&dis_ctx);
	VecScatterDestroy(&trac_ctx);
	VecDestroy(&trac_collection);
	VecDestroy(&dis_collection);

	//OutPutTecPlot();
	OutPutVTK();
}
