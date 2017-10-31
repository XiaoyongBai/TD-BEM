/*
 * SolverMRT.cpp
 *
 *  Created on: Oct 25, 2015
 *      Author: xiaoyong
 */

#include "SolverMRT.h"
#include "MathOperationT.h"

#include <iostream>
#include <fstream>
#include <iomanip>

#include <string>
#include <sstream>

using namespace std;
using namespace TD_BEM;
using namespace GreenFunction;

SolverMRT::SolverMRT(ModelManagerT* model_1, ModelManagerT* model_2) {

	fModel_1=model_1;
	fModel_2=model_2;

	fComm=PETSC_COMM_WORLD;

	MPI_Comm_size(fComm, &fSize);
	MPI_Comm_rank(fComm, &fRank);

	LocalIndexing();
	InterfaceIndexing();

	fMaxStep=10;
	fBeta=0;
	fDt=model_1->GetDT();
	fNumStep=model_1->GetNumStep();

	fCurrStep=0;

	/**
	 * global matrices and field variables
	 */
	fDis_DB_1=NULL;
	fTrac_DB_1=NULL;

	fDis_DB_2=NULL;
	fTrac_DB_2=NULL;

	fDis_DB_Global=NULL;
	fTrac_DB_Global=NULL;

	fG_DB=NULL;
	fH_DB=NULL;
	fG_DB_bgn=NULL;
	fH_DB_bgn=NULL;

	fIerr=KSPCreate(fComm, &fKsp);

	cout<<"Solver is constructed\n";
}




SolverMRT::~SolverMRT() {
	delete[] fExtraction_row_1;
	delete[] fExtraction_row_2;

	delete[] fExtraction_Hr_col_1;
	delete[] fExtraction_Hr_col_2;

	delete[] fExtraction_Hz_col_1;
	delete[] fExtraction_Hz_col_2;

	delete[] fExtraction_Gr_col_1;
	delete[] fExtraction_Gr_col_2;

	delete[] fExtraction_Gz_col_1;
	delete[] fExtraction_Gz_col_2;


	delete[] fAssemble_row_1;
	delete[] fAssemble_row_2;

	delete[] fAssemble_Hr_col_1;
	delete[] fAssemble_Hr_col_2;

	delete[] fAssemble_Hz_col_1;
	delete[] fAssemble_Hz_col_2;

	delete[] fAssemble_Gr_col_1;
	delete[] fAssemble_Gr_col_2;

	delete[] fAssemble_Gz_col_1;
	delete[] fAssemble_Gz_col_2;

	delete[] fExtraction_Dis_1;
	delete[] fExtraction_Dis_2;

	delete[] fSet_Dis_1;
	delete[] fSet_Dis_2;

	/**************************************************************************
	 **
	 ** Release the memory
	 **
	 * ************************************************************************/
	for (int i=0; i<fNumStep+1; i++)
	{
		VecDestroy(&fDis_DB_1[i]);
		VecDestroy(&fTrac_DB_1[i]);

		VecDestroy(&fDis_DB_2[i]);
		VecDestroy(&fTrac_DB_2[i]);

		VecDestroy(&fDis_DB_Global[i]);
		VecDestroy(&fTrac_DB_Global[i]);
	}

	delete[] fDis_DB_1;
	delete[] fTrac_DB_1;
	delete[] fDis_DB_2;
	delete[] fTrac_DB_2;
	delete[] fDis_DB_Global;
	delete[] fTrac_DB_Global;

	for (int i=0; i<fMaxStep+1; i++)
	{
		MatDestroy(&fG_DB[i]);
		MatDestroy(&fH_DB[i]);
	}

	delete[] fG_DB;
	delete[] fH_DB;

	for (int i=0; i<3; i++)
	{
		MatDestroy(&fG_DB_bgn[i]);
		MatDestroy(&fH_DB_bgn[i]);
	}
	delete[] fG_DB_bgn;
	delete[] fH_DB_bgn;

	MatDestroy(&fG_Compute);
	MatDestroy(&fH_Compute);

	VecDestroy(&fLoad);
	VecDestroy(&fRHS);

	KSPDestroy(&fKsp);

}


void SolverMRT::LocalIndexing()
{
	int sd=fModel_1->GetSD();

	/**********************************************************
	 **
	 ** get dimensions for region 1
	 **
	 **********************************************************/
	fNumNodes_1=fModel_1->GetNND();
	fNumDofs_1=sd*fNumNodes_1;

	fNumInterfaceNodes=fModel_1->GetNumInterface();
	fNumInterfaceDofs=sd*fNumInterfaceNodes;

	fNumInteriorNodes_1=fNumNodes_1-fNumInterfaceNodes;
	fNumInteriorDofs_1=fNumInteriorNodes_1*sd;

	Vec NodeIndex;
	VecCreateMPI(fComm, PETSC_DECIDE, fNumNodes_1, &NodeIndex);
	VecSet(NodeIndex, 0);
	VecGetOwnershipRange(NodeIndex, &fNodeLow_1, &fNodeHigh_1);
	VecDestroy(&NodeIndex);

	fNodeHigh_1=fNodeHigh_1-1;
	fNumRow_Local_1=(fNodeHigh_1-fNodeLow_1+1)*sd;

	fDofLow_1=3*fNodeLow_1;
	fDofHigh_1=3*fNodeHigh_1+2;

	/**********************************************************
	 **
	 ** get dimensions for region 2
	 **
	 **********************************************************/
	fNumNodes_2=fModel_2->GetNND();
	fNumDofs_2=sd*fNumNodes_2;

	fNumInteriorNodes_2=fNumNodes_2-fNumInterfaceNodes;
	fNumInteriorDofs_2=fNumInteriorNodes_2*sd;

	VecCreateMPI(fComm, PETSC_DECIDE, fNumNodes_2, &NodeIndex);
	VecSet(NodeIndex, 0);
	VecGetOwnershipRange(NodeIndex, &fNodeLow_2, &fNodeHigh_2);
	VecDestroy(&NodeIndex);

	fNodeHigh_2=fNodeHigh_2-1;

	fNumRow_Local_2=(fNodeHigh_2-fNodeLow_2+1)*sd;

	fDofLow_2=3*fNodeLow_2;
	fDofHigh_2=3*fNodeHigh_2+2;

	/*****************************************************************
	 *
	 * set the number of global DOFs
	 *
	 ******************************************************************/
	fNumRow_H=fNumDofs_1+fNumDofs_2;
	fNumCol_H=fNumRow_H;

	fNumRow_G=fNumRow_H;
	fNumCol_G=fNumInteriorDofs_1+fNumInteriorDofs_2;
}


void SolverMRT::InterfaceIndexing()
{
	int sd=fModel_1->GetSD();

	/******************************************************************
	 *
	 * Form indexing for region 1
	 *
	 ******************************************************************/

	//find interior nodes for region 1
	const int* InterfaceNodes_1=fModel_1->GetInterface();

	int* lable_1=new int[fNumNodes_1];
	MathOperationT::VecSet(fNumNodes_1, lable_1, 0);
	for (int i=0; i<fNumInterfaceNodes; i++)
		lable_1[InterfaceNodes_1[i]]=1;

	int* InteriorNodes_1=new int[fNumInteriorNodes_1];
	int temp_1=0;
	for (int i=0; i<fNumNodes_1; i++)
	{
		if (lable_1[i]==0)
		{
			InteriorNodes_1[temp_1]=i;
			temp_1++;
		}
	}

	//Construct interface dofs for region 1
	int* InterfaceDOFs_1=new int[sd*fNumInterfaceNodes];

	for(int ni=0; ni<fNumInterfaceNodes;ni++)
	{
		for(int d=0; d<sd; d++)
		{
			InterfaceDOFs_1[ni*sd+d]=InterfaceNodes_1[ni]*sd+d;
		}
	}

	//construct interior DOFs for region 1
	int* InteriorDOFs_1=new int[sd*fNumInteriorNodes_1];

	for(int ni=0; ni<fNumInteriorNodes_1;ni++)
	{
		for(int d=0; d<sd; d++)
		{
			InteriorDOFs_1[ni*sd+d]=InteriorNodes_1[ni]*sd+d;
		}
	}

	//MathOperationT::PrintVector(fNumInterfaceNodes*sd, InterfaceDOFs_1, "Interface dofs in region 1");
	//MathOperationT::PrintVector(fNumInteriorNodes_1, InteriorNodes_1, "Interior nodes in region 1");
	//MathOperationT::PrintVector(fNumInteriorDofs_1, InteriorDOFs_1, "Interior dofs in region 1");


	//form extraction indices for region 1
	fExtraction_row_1=new int[fNumRow_Local_1];
	for(int ri=0; ri<fNumRow_Local_1; ri++)
		fExtraction_row_1[ri]=ri;

	//MathOperationT::PrintVector(fNumRow_Local_1, fExtraction_row_1, "fExtraction_row_1");


	fExtraction_Hr_col_1=new int[fNumInteriorDofs_1];
	for(int ci=0; ci<fNumInteriorDofs_1; ci++)
		fExtraction_Hr_col_1[ci]=InteriorDOFs_1[ci];

	//MathOperationT::PrintVector(fNumInteriorDofs_1, fExtraction_Hr_col_1, "fExtraction_Hr_col_1");

	fExtraction_Hz_col_1=new int[fNumInterfaceDofs];
	for(int ci=0; ci<fNumInterfaceDofs; ci++)
		fExtraction_Hz_col_1[ci]=InterfaceDOFs_1[ci];

	//MathOperationT::PrintVector(fNumInterfaceDofs, fExtraction_Hz_col_1, "fExtraction_Hz_col_1");

	fExtraction_Gr_col_1=new int[fNumInteriorDofs_1];
	for(int ci=0; ci<fNumInteriorDofs_1; ci++)
		fExtraction_Gr_col_1[ci]=InteriorDOFs_1[ci];

	//MathOperationT::PrintVector(fNumInteriorDofs_1, fExtraction_Gr_col_1, "fExtraction_Gr_col_1");

	fExtraction_Gz_col_1=new int[fNumInterfaceDofs];
	for(int ci=0; ci<fNumInterfaceDofs; ci++)
		fExtraction_Gz_col_1[ci]=InterfaceDOFs_1[ci];

	//MathOperationT::PrintVector(fNumInterfaceDofs, fExtraction_Gz_col_1, "fExtraction_Gz_col_1");


	//form assembly indices for region 1
	fAssemble_row_1=new int[fNumRow_Local_1];
	for(int mi=fDofLow_1*sd; mi<=fDofHigh_1; mi++)
		fAssemble_row_1[mi-fDofLow_1]=mi;

	//MathOperationT::PrintVector(fNumRow_Local_1, fAssemble_row_1, "fAssemble_row_1");


	fAssemble_Hr_col_1=new int[fNumInteriorDofs_1];
	for(int ri=0; ri<fNumInteriorDofs_1; ri++)
		fAssemble_Hr_col_1[ri]=ri;

	//MathOperationT::PrintVector(fNumInteriorDofs_1, fAssemble_Hr_col_1, "fAssemble_Hr_col_1");

	fAssemble_Hz_col_1=new int[fNumInterfaceDofs];
	for(int ri=0; ri<fNumInterfaceDofs; ri++)
		fAssemble_Hz_col_1[ri]=ri+fNumInteriorDofs_1+fNumInteriorDofs_2;

	//MathOperationT::PrintVector(fNumInterfaceDofs, fAssemble_Hz_col_1, "fAssemble_Hz_col_1");

	fAssemble_Gr_col_1=new int[fNumInteriorDofs_1];
	for(int ri=0; ri<fNumInteriorDofs_1; ri++)
		fAssemble_Gr_col_1[ri]=ri;

	//MathOperationT::PrintVector(fNumInteriorDofs_1, fAssemble_Gr_col_1, "fAssemble_Gr_col_1");

	fAssemble_Gz_col_1=new int[fNumInterfaceDofs];
	for(int ri=0; ri<fNumInterfaceDofs; ri++)
		fAssemble_Gz_col_1[ri]=ri+fNumInteriorDofs_1+fNumInteriorDofs_2+fNumInterfaceDofs;

	//MathOperationT::PrintVector(fNumInterfaceDofs, fAssemble_Gz_col_1, "fAssemble_Gz_col_1");


	//form indices for extracting field vairables
	fExtraction_Dis_1=new int[fNumDofs_1];
	for(int i=0; i<fNumInteriorDofs_1; i++)
	{
		fExtraction_Dis_1[InteriorDOFs_1[i]]=i;
	}
	for(int i=0; i<fNumInterfaceDofs; i++)
	{
		fExtraction_Dis_1[InterfaceDOFs_1[i]]=i+fNumInteriorDofs_1+fNumInteriorDofs_2;
	}

	fSet_Dis_1=new int[fNumDofs_1];
	for(int i=0; i<fNumDofs_1; i++)
		fSet_Dis_1[i]=i;

	//MathOperationT::PrintVector(fNumDofs_1, fExtraction_Dis_1, "fExtraction_Dis_1");
	//MathOperationT::PrintVector(fNumDofs_1, fSet_Dis_1, "fSet_Dis_1");


	/******************************************************************
	 *
	 * Form indices for region 2
	 *
	 ******************************************************************/

	//find interior nodes for region 2
	const int* InterfaceNodes_2=fModel_2->GetInterface();

	int* lable_2=new int[fNumNodes_2];
	MathOperationT::VecSet(fNumNodes_2, lable_2, 0);
	for (int i=0; i<fNumInterfaceNodes; i++)
		lable_2[InterfaceNodes_2[i]]=1;

	int* InteriorNodes_2=new int[fNumInteriorNodes_2];
	int temp_2=0;
	for (int i=0; i<fNumNodes_2; i++)
	{
		if (lable_2[i]==0)
		{
			InteriorNodes_2[temp_2]=i;
			temp_2++;
		}
	}

	//Construct interface DOFs for region 2
	int* InterfaceDOFs_2=new int[sd*fNumInterfaceNodes];

	for(int ni=0; ni<fNumInterfaceNodes;ni++)
	{
		for(int d=0; d<sd; d++)
		{
			InterfaceDOFs_2[ni*sd+d]=InterfaceNodes_2[ni]*sd+d;
		}
	}

	//construct interior DOFs for region 2
	int* InteriorDOFs_2=new int[sd*fNumInteriorNodes_2];

	for(int ni=0; ni<fNumInteriorNodes_2;ni++)
	{
		for(int d=0; d<sd; d++)
		{
			InteriorDOFs_2[ni*sd+d]=InteriorNodes_2[ni]*sd+d;
		}
	}

	//MathOperationT::PrintVector(fNumInterfaceNodes*sd, InterfaceDOFs_2, "Interface dofs in region 2");
	//MathOperationT::PrintVector(fNumInteriorNodes_2, InteriorNodes_2, "Interior nodes in region 2");
	//MathOperationT::PrintVector(fNumInteriorDofs_2, InteriorDOFs_2, "Interior dofs in region 2");


	//form extraction indices for region 2
	fExtraction_row_2=new int[fNumRow_Local_2];
	for(int ri=0; ri<fNumRow_Local_2; ri++)
		fExtraction_row_2[ri]=ri;

	//MathOperationT::PrintVector(fNumRow_Local_2, fExtraction_row_2, "fExtraction_row_2");

	fExtraction_Hr_col_2=new int[fNumInteriorDofs_2];
	for(int ci=0; ci<fNumInteriorDofs_2; ci++)
		fExtraction_Hr_col_2[ci]=InteriorDOFs_2[ci];

	//MathOperationT::PrintVector(fNumInteriorDofs_2, fExtraction_Hr_col_2, "fExtraction_Hr_col_2");

	fExtraction_Hz_col_2=new int[fNumInterfaceDofs];
	for(int ci=0; ci<fNumInterfaceDofs; ci++)
		fExtraction_Hz_col_2[ci]=InterfaceDOFs_2[ci];

	//MathOperationT::PrintVector(fNumInterfaceDofs, fExtraction_Hz_col_2, "fExtraction_Hz_col_2");

	fExtraction_Gr_col_2=new int[fNumInteriorDofs_2];
	for(int ci=0; ci<fNumInteriorDofs_2; ci++)
		fExtraction_Gr_col_2[ci]=InteriorDOFs_2[ci];

	//MathOperationT::PrintVector(fNumInteriorDofs_2, fExtraction_Gr_col_2, "fExtraction_Gr_col_2");

	fExtraction_Gz_col_2=new int[fNumInterfaceDofs];
	for(int ci=0; ci<fNumInterfaceDofs; ci++)
		fExtraction_Gz_col_2[ci]=InterfaceDOFs_2[ci];

	//MathOperationT::PrintVector(fNumInterfaceDofs, fExtraction_Gz_col_2, "fExtraction_Gz_col_2");


	/******************************************************************
	 *
	 * Find identical nodes in region 1 for interface nodes in region 2
	 *
	 ******************************************************************/
	int* Interface_node_2_to_1= new int[fNumInterfaceNodes];
	int* Interface_dof_2_to_1=new int[fNumInterfaceDofs];

	const double* coords_1=fModel_1->GetCoords();
	const double* coords_2=fModel_2->GetCoords();

	for(int ni=0; ni<fNumInterfaceNodes; ni++)
	{
		int ii=InterfaceNodes_2[ni]
;		double x2=coords_2[ii*sd];
		double y2=coords_2[ii*sd+1];
		double z2=coords_2[ii*sd+2];

		int match=-1;

		for(int nj=0; nj<fNumInterfaceNodes; nj++)
		{
			int ij=InterfaceNodes_1[nj];
			double x1=coords_1[ij*sd];
			double y1=coords_1[ij*sd+1];
			double z1=coords_1[ij*sd+2];

			double r=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
			if(r<1e-10)
			{
				match=nj;
				break;
			}
		}

		if(match==-1)
		{
			cout<<"can't find a match point for interface node " << ni << " in region 2 \n";
			throw "Interface mismatch, interface nodes must match one-to-one\n";
		}
		Interface_node_2_to_1[ni]=match;

		Interface_dof_2_to_1[ni*sd]=3*match;
		Interface_dof_2_to_1[ni*sd+1]=3*match+1;
		Interface_dof_2_to_1[ni*sd+2]=3*match+2;

	}

	//MathOperationT::PrintVector(fNumInterfaceNodes, Interface_node_2_to_1, "Interface nodes from 2 to 1");
	//MathOperationT::PrintVector(fNumInterfaceDofs, Interface_dof_2_to_1, "Interface Dofs from 2 to 1");



	//form assembly indices for region 2
	fAssemble_row_2=new int[fNumRow_Local_2];
	for(int mi=fDofLow_2; mi<=fDofHigh_2; mi++)
		fAssemble_row_2[mi-fDofLow_2]=mi+fNumDofs_1;

	//MathOperationT::PrintVector(fNumRow_Local_2, fAssemble_row_2, "fAssemble_row_2");


	fAssemble_Hr_col_2=new int[fNumInteriorDofs_2];
	for(int ri=0; ri<fNumInteriorDofs_2; ri++)
		fAssemble_Hr_col_2[ri]=ri+fNumInteriorDofs_1;

	//MathOperationT::PrintVector(fNumInteriorDofs_2, fAssemble_Hr_col_2, "fAssemble_Hr_col_2");

	fAssemble_Hz_col_2=new int[fNumInterfaceDofs];
	for(int ri=0; ri<fNumInterfaceDofs; ri++)
		fAssemble_Hz_col_2[ri]=Interface_dof_2_to_1[ri]+fNumInteriorDofs_1+fNumInteriorDofs_2;

	//MathOperationT::PrintVector(fNumInterfaceDofs, fAssemble_Hz_col_2, "fAssemble_Hz_col_2");

	fAssemble_Gr_col_2=new int[fNumInteriorDofs_2];
	for(int ri=0; ri<fNumInteriorDofs_2; ri++)
		fAssemble_Gr_col_2[ri]=ri+fNumInteriorDofs_1;

	//MathOperationT::PrintVector(fNumInteriorDofs_2, fAssemble_Gr_col_2, "fAssemble_Gr_col_2");

	fAssemble_Gz_col_2=new int[fNumInterfaceDofs];
	for(int ri=0; ri<fNumInterfaceDofs; ri++)
		fAssemble_Gz_col_2[ri]=Interface_dof_2_to_1[ri]+fNumInteriorDofs_1+fNumInteriorDofs_2+fNumInterfaceDofs;

	//MathOperationT::PrintVector(fNumInterfaceDofs, fAssemble_Gz_col_2, "fAssemble_Gz_col_2");


	//form indices for extracting field variables
	fExtraction_Dis_2=new int[fNumDofs_2];
	for(int i=0; i<fNumInteriorDofs_2; i++)
	{
		fExtraction_Dis_2[InteriorDOFs_2[i]]=i+fNumInteriorDofs_1;
	}
	for(int i=0; i<fNumInterfaceDofs; i++)
	{
		fExtraction_Dis_2[InterfaceDOFs_2[i]]=Interface_dof_2_to_1[i]+fNumInteriorDofs_1+fNumInteriorDofs_2;
	}

	fSet_Dis_2=new int[fNumDofs_2];
	for(int i=0; i<fNumDofs_2; i++)
		fSet_Dis_2[i]=i;

	//MathOperationT::PrintVector(fNumDofs_2, fExtraction_Dis_2, "fExtraction_Dis_2");
	//MathOperationT::PrintVector(fNumDofs_2, fSet_Dis_2, "fSet_Dis_2");


	delete[] lable_1;
	delete[] InterfaceDOFs_1;
	delete[] InteriorNodes_1;
	delete[] InteriorDOFs_1;

	delete[] lable_2;
	delete[] InterfaceDOFs_2;
	delete[] InteriorNodes_2;
	delete[] InteriorDOFs_2;

	delete[] Interface_node_2_to_1;
	delete[] Interface_dof_2_to_1;

}



void SolverMRT::SetTimeControl(double dt, int Numsteps)
{
	fDt=dt;
	fNumStep=Numsteps;

	/**********************************************************
	 **
	 ** Allocate field variables
	 **
	 **********************************************************/
	fDis_DB_1= new Vec[fNumStep+1];
	fTrac_DB_1= new Vec[fNumStep+1];

	fDis_DB_2= new Vec[fNumStep+1];
	fTrac_DB_2= new Vec[fNumStep+1];

	fDis_DB_Global= new Vec[fNumStep+1];
	fTrac_DB_Global= new Vec[fNumStep+1];

	Vec* temp_d=fDis_DB_1;
	Vec* temp_t=fTrac_DB_1;

	for (int i=0; i<fNumStep+1; i++)
	{
		fIerr=VecCreateMPI(PETSC_COMM_WORLD,  PETSC_DECIDE, fNumDofs_1, temp_d);
		fIerr=VecCreateMPI(PETSC_COMM_WORLD,  PETSC_DECIDE, fNumDofs_1, temp_t);
		VecSet(*temp_d, 0);
		VecSet(*temp_t, 0);
		temp_d++;
		temp_t++;
	}

	temp_d=fDis_DB_2;
	temp_t=fTrac_DB_2;

	for (int i=0; i<fNumStep+1; i++)
	{
		fIerr=VecCreateMPI(PETSC_COMM_WORLD,  PETSC_DECIDE, fNumDofs_2, temp_d);
		fIerr=VecCreateMPI(PETSC_COMM_WORLD,  PETSC_DECIDE, fNumDofs_2, temp_t);
		VecSet(*temp_d, 0);
		VecSet(*temp_t, 0);
		temp_d++;
		temp_t++;
	}

	temp_d=fDis_DB_Global;
	temp_t=fTrac_DB_Global;

	for (int i=0; i<fNumStep+1; i++)
	{
		fIerr=VecCreateMPI(PETSC_COMM_WORLD,  PETSC_DECIDE, fNumCol_H, temp_d);
		fIerr=VecCreateMPI(PETSC_COMM_WORLD,  PETSC_DECIDE, fNumCol_G, temp_t);
		VecSet(*temp_d, 0);
		VecSet(*temp_t, 0);
		temp_d++;
		temp_t++;
	}


	/**********************************************************
	 **
	 ** Allocate Matrices
	 **
	 **********************************************************/
	fG_DB=new Mat[fMaxStep+1];
	fH_DB=new Mat[fMaxStep+1];

	Mat* temp_G=fG_DB;
	Mat* temp_H=fH_DB;

	for (int i=0; i<fMaxStep+1; i++)
	{
		fIerr=MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, fNumRow_G, fNumCol_G, PETSC_NULL, temp_G);
		fIerr=MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, fNumRow_H, fNumCol_H, PETSC_NULL, temp_H);
		temp_G++;
		temp_H++;
	}

	fG_DB_bgn=new Mat[3];
	fH_DB_bgn=new Mat[3];

	temp_G=fG_DB_bgn;
	temp_H=fH_DB_bgn;

	for (int i=0; i<3; i++)
	{
		fIerr=MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, fNumRow_G, fNumCol_G, PETSC_NULL, temp_G);
		fIerr=MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, fNumRow_H, fNumCol_H, PETSC_NULL, temp_H);
		temp_G++;
		temp_H++;
	}

	/**********************************************************
	 **
	 ** Allocate computing matrices
	 **
	 **********************************************************/
	fIerr=MatCreateDense(fComm, PETSC_DECIDE, PETSC_DECIDE, fNumRow_G, fNumCol_G, PETSC_NULL, &fG_Compute);
	fIerr=MatCreateDense(fComm, PETSC_DECIDE, PETSC_DECIDE, fNumRow_H, fNumCol_H, PETSC_NULL, &fH_Compute);
	fIerr=VecCreateMPI(fComm,  PETSC_DECIDE, fNumCol_G, &fLoad);
	fIerr=VecCreateMPI(fComm,  PETSC_DECIDE, fNumCol_H, &fRHS);
}



void SolverMRT::GetNodeLowHigh(int* low_1, int* high_1, int* low_2, int* high_2)
{
	*low_1=fNodeLow_1;
	*high_1=fNodeHigh_1;

	*low_2=fNodeLow_2;
	*high_2=fNodeHigh_2;
}


void SolverMRT::Assemble_GH_R1(Mat& G_global, Mat& H_global, const double* G_local, const double* H_local, int insert_or_add)
{
	/**************************************************
	 * decompose G into 2 parts
	 * ************************************************/
	double* G_interior=new double[fNumRow_Local_1*fNumInteriorDofs_1];
	double* G_interface=new double[fNumRow_Local_1*fNumInterfaceDofs];

	MathOperationT::MatExtraction(fNumRow_Local_1, fNumDofs_1, G_local, fNumRow_Local_1, fExtraction_row_1, fNumInteriorDofs_1, fExtraction_Gr_col_1,G_interior);
	MathOperationT::MatExtraction(fNumRow_Local_1, fNumDofs_1, G_local, fNumRow_Local_1, fExtraction_row_1, fNumInterfaceDofs, fExtraction_Gz_col_1,G_interface);

	/**************************************************
	 * decompose H into 2 parts
	 * ************************************************/
	double* H_interior=new double[fNumRow_Local_1*fNumInteriorDofs_1];
	double* H_interface=new double[fNumRow_Local_1*fNumInterfaceDofs];

	MathOperationT::MatExtraction(fNumRow_Local_1, fNumDofs_1, H_local, fNumRow_Local_1, fExtraction_row_1, fNumInteriorDofs_1, fExtraction_Hr_col_1, H_interior);
	MathOperationT::MatExtraction(fNumRow_Local_1, fNumDofs_1, H_local, fNumRow_Local_1, fExtraction_row_1, fNumInterfaceDofs, fExtraction_Hz_col_1, H_interface);


	/**************************************************
	 * Assemble G
	 * ************************************************/
	if(insert_or_add==0) //insert
	{
		fIerr=MatSetValues(G_global, fNumRow_Local_1, fAssemble_row_1, fNumInteriorDofs_1, fAssemble_Gr_col_1, G_interior, INSERT_VALUES);
		fIerr=MatAssemblyBegin(G_global, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(G_global, MAT_FINAL_ASSEMBLY);
	}
	else if(insert_or_add==1) //add
	{
		fIerr=MatSetValues(G_global, fNumRow_Local_1, fAssemble_row_1, fNumInteriorDofs_1, fAssemble_Gr_col_1, G_interior, ADD_VALUES);
		fIerr=MatAssemblyBegin(G_global, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(G_global, MAT_FINAL_ASSEMBLY);
	}
	else
	{
		cout<<"SolverMRT::Assemble_GH_R1, invalid model of mat set values"<<endl;
		throw "SolverMRT::Assemble_GH_R1, invalid";
	}


	/**************************************************
	 * Assemble H
	 * ************************************************/
	MathOperationT::VecScale(fNumRow_Local_1*fNumInterfaceDofs, G_interface,-1);

	if(insert_or_add==0)
	{
		fIerr=MatSetValues(H_global, fNumRow_Local_1, fAssemble_row_1, fNumInteriorDofs_1, fAssemble_Hr_col_1, H_interior, INSERT_VALUES);
		fIerr=MatAssemblyBegin(H_global, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(H_global, MAT_FINAL_ASSEMBLY);

		fIerr=MatSetValues(H_global, fNumRow_Local_1, fAssemble_row_1, fNumInterfaceDofs, fAssemble_Hz_col_1, H_interface, INSERT_VALUES);
		fIerr=MatAssemblyBegin(H_global, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(H_global, MAT_FINAL_ASSEMBLY);

		fIerr=MatSetValues(H_global, fNumRow_Local_1, fAssemble_row_1, fNumInterfaceDofs, fAssemble_Gz_col_1, G_interface, INSERT_VALUES);
		fIerr=MatAssemblyBegin(H_global, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(H_global, MAT_FINAL_ASSEMBLY);
	}
	else
	{
		fIerr=MatSetValues(H_global, fNumRow_Local_1, fAssemble_row_1, fNumInteriorDofs_1, fAssemble_Hr_col_1, H_interior, ADD_VALUES);
		fIerr=MatAssemblyBegin(H_global, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(H_global, MAT_FINAL_ASSEMBLY);

		fIerr=MatSetValues(H_global, fNumRow_Local_1, fAssemble_row_1, fNumInterfaceDofs, fAssemble_Hz_col_1, H_interface, ADD_VALUES);
		fIerr=MatAssemblyBegin(H_global, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(H_global, MAT_FINAL_ASSEMBLY);

		fIerr=MatSetValues(H_global, fNumRow_Local_1, fAssemble_row_1, fNumInterfaceDofs, fAssemble_Gz_col_1, G_interface, ADD_VALUES);
		fIerr=MatAssemblyBegin(H_global, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(H_global, MAT_FINAL_ASSEMBLY);
	}


	delete[] G_interior;
	delete[] G_interface;
	delete[] H_interior;
	delete[] H_interface;
}


void SolverMRT::Assemble_GH_R2(Mat& G_global, Mat& H_global, const double* G_local, const double* H_local, int insert_or_add)
{
	/**************************************************
	 * decompose G into 2 parts
	 * ************************************************/
	double* G_interior=new double[fNumRow_Local_2*fNumInteriorDofs_2];
	double* G_interface=new double[fNumRow_Local_2*fNumInterfaceDofs];

	MathOperationT::MatExtraction(fNumRow_Local_2, fNumDofs_2, G_local, fNumRow_Local_2, fExtraction_row_2, fNumInteriorDofs_2, fExtraction_Gr_col_2,G_interior);
	MathOperationT::MatExtraction(fNumRow_Local_2, fNumDofs_2, G_local, fNumRow_Local_2, fExtraction_row_2, fNumInterfaceDofs, fExtraction_Gz_col_2,G_interface);

	/**************************************************
	 * decompose H into 2 parts
	 * ************************************************/
	double* H_interior=new double[fNumRow_Local_2*fNumInteriorDofs_2];
	double* H_interface=new double[fNumRow_Local_2*fNumInterfaceDofs];

	MathOperationT::MatExtraction(fNumRow_Local_2, fNumDofs_2, H_local, fNumRow_Local_2, fExtraction_row_2, fNumInteriorDofs_2, fExtraction_Hr_col_2, H_interior);
	MathOperationT::MatExtraction(fNumRow_Local_2, fNumDofs_2, H_local, fNumRow_Local_2, fExtraction_row_2, fNumInterfaceDofs, fExtraction_Hz_col_2, H_interface);


	/**************************************************
	 * Assemble G
	 * ************************************************/
	if(insert_or_add==0) //insert
	{
		fIerr=MatSetValues(G_global, fNumRow_Local_2, fAssemble_row_2, fNumInteriorDofs_2, fAssemble_Gr_col_2, G_interior, INSERT_VALUES);
		fIerr=MatAssemblyBegin(G_global, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(G_global, MAT_FINAL_ASSEMBLY);
	}
	else if(insert_or_add==1) //add
	{
		fIerr=MatSetValues(G_global, fNumRow_Local_2, fAssemble_row_2, fNumInteriorDofs_2, fAssemble_Gr_col_2, G_interior, ADD_VALUES);
		fIerr=MatAssemblyBegin(G_global, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(G_global, MAT_FINAL_ASSEMBLY);
	}
	else
	{
		cout<<"SolverMRT::Assemble_GH_R2, invalid model of mat set values"<<endl;
		throw "SolverMRT::Assemble_GH_R2, invalid";
	}


	/**************************************************
	 * Assemble H
	 * ************************************************/
	if(insert_or_add==0)
	{
		fIerr=MatSetValues(H_global, fNumRow_Local_2, fAssemble_row_2, fNumInteriorDofs_2, fAssemble_Hr_col_2, H_interior, INSERT_VALUES);
		fIerr=MatAssemblyBegin(H_global, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(H_global, MAT_FINAL_ASSEMBLY);

		fIerr=MatSetValues(H_global, fNumRow_Local_2, fAssemble_row_2, fNumInterfaceDofs, fAssemble_Hz_col_2, H_interface, INSERT_VALUES);
		fIerr=MatAssemblyBegin(H_global, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(H_global, MAT_FINAL_ASSEMBLY);

		fIerr=MatSetValues(H_global, fNumRow_Local_2, fAssemble_row_2, fNumInterfaceDofs, fAssemble_Gz_col_2, G_interface, INSERT_VALUES);
		fIerr=MatAssemblyBegin(H_global, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(H_global, MAT_FINAL_ASSEMBLY);
	}
	else
	{
		fIerr=MatSetValues(H_global, fNumRow_Local_2, fAssemble_row_2, fNumInteriorDofs_2, fAssemble_Hr_col_2, H_interior, ADD_VALUES);
		fIerr=MatAssemblyBegin(H_global, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(H_global, MAT_FINAL_ASSEMBLY);

		fIerr=MatSetValues(H_global, fNumRow_Local_2, fAssemble_row_2, fNumInterfaceDofs, fAssemble_Hz_col_2, H_interface, ADD_VALUES);
		fIerr=MatAssemblyBegin(H_global, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(H_global, MAT_FINAL_ASSEMBLY);

		fIerr=MatSetValues(H_global, fNumRow_Local_2, fAssemble_row_2, fNumInterfaceDofs, fAssemble_Gz_col_2, G_interface, ADD_VALUES);
		fIerr=MatAssemblyBegin(H_global, MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(H_global, MAT_FINAL_ASSEMBLY);
	}

	delete[] G_interior;
	delete[] G_interface;
	delete[] H_interior;
	delete[] H_interface;
}


void SolverMRT::UpdateGHLinear(const double* G_0, const double* G_1, const double* H_0, const double* H_1,
		                         const double* g_0, const double* g_1, const double* h_0, const double* h_1)
{
	/*************
	 * Arrange the matrices
	 *************/
	if (fCurrStep<=fMaxStep)
	{

		if(fCurrStep==1)
		{
			Assemble_GH_R1(fG_DB[0],fH_DB[0], G_1, H_1, 0);
			Assemble_GH_R2(fG_DB[0],fH_DB[0], g_1, h_1, 0);
		}
		else
		{
			Assemble_GH_R1(fG_DB[fCurrStep-1],fH_DB[fCurrStep-1], G_1, H_1, 1);
			Assemble_GH_R2(fG_DB[fCurrStep-1],fH_DB[fCurrStep-1], g_1, h_1, 1);
		}

		Assemble_GH_R1(fG_DB[fCurrStep],fH_DB[fCurrStep], G_0, H_0, 0);
		Assemble_GH_R2(fG_DB[fCurrStep],fH_DB[fCurrStep], g_0, h_0, 0);
	}


	//set fG_DB_bgn
	if (fCurrStep>fMaxStep+3)
		return;

	if(fCurrStep==1)
	{
		Assemble_GH_R1(fG_DB_bgn[0],fH_DB_bgn[0], G_0, H_0, 0);
		Assemble_GH_R2(fG_DB_bgn[0],fH_DB_bgn[0], g_0, h_0, 0);
	}
	else if(fCurrStep==2)
	{
		fIerr=MatCopy(fG_DB_bgn[0], fG_DB_bgn[1], DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAssemblyBegin(fG_DB_bgn[1], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fG_DB_bgn[1], MAT_FINAL_ASSEMBLY);

		fIerr=MatCopy(fH_DB_bgn[0], fH_DB_bgn[1], DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAssemblyBegin(fH_DB_bgn[1], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fH_DB_bgn[1], MAT_FINAL_ASSEMBLY);

		Assemble_GH_R1(fG_DB_bgn[0],fH_DB_bgn[0], G_0, H_0, 0);
		Assemble_GH_R2(fG_DB_bgn[0],fH_DB_bgn[0], g_0, h_0, 0);
	}
	else
	{
		fIerr=MatCopy(fG_DB_bgn[1], fG_DB_bgn[2], DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAssemblyBegin(fG_DB_bgn[2], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fG_DB_bgn[2], MAT_FINAL_ASSEMBLY);

		fIerr=MatCopy(fG_DB_bgn[0], fG_DB_bgn[1], DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAssemblyBegin(fG_DB_bgn[1], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fG_DB_bgn[1], MAT_FINAL_ASSEMBLY);

		/********************************************************************
		 ********************************************************************/

		fIerr=MatCopy(fH_DB_bgn[1], fH_DB_bgn[2], DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAssemblyBegin(fH_DB_bgn[2], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fH_DB_bgn[2], MAT_FINAL_ASSEMBLY);

		fIerr=MatCopy(fH_DB_bgn[0], fH_DB_bgn[1], DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAssemblyBegin(fH_DB_bgn[1], MAT_FINAL_ASSEMBLY);
		fIerr=MatAssemblyEnd(fH_DB_bgn[1], MAT_FINAL_ASSEMBLY);

		Assemble_GH_R1(fG_DB_bgn[0],fH_DB_bgn[0], G_0, H_0, 0);
		Assemble_GH_R2(fG_DB_bgn[0],fH_DB_bgn[0], g_0, h_0, 0);

		if (fCurrStep>fMaxStep)
		{
			fIerr=MatScale(fH_DB_bgn[0], 0);
			fIerr=MatScale(fG_DB_bgn[0], 0);
		}
	}

}




void SolverMRT::SetLoad()
{
	//****************************************
	double p0=5.46366e9;
	double curr_time=fCurrStep*fDt;
	double amp=p0*curr_time*exp(-2000*curr_time);

	//amp=1e6;


	int UBC_num;
	int* UBC_DOFs;
	double* UBC_Values;
	int FBC_num;
	int* FBC_DOFs;
	double* FBC_Values;

	fModel_1->BoundaryCondition_Global(UBC_num, &UBC_DOFs, &UBC_Values,
			                          FBC_num, &FBC_DOFs, &FBC_Values);

	double* FBC_Values_temp_1=new double[FBC_num];
	MathOperationT::MemCopy(FBC_num, FBC_Values, FBC_Values_temp_1);

	MathOperationT::VecScale(FBC_num, FBC_Values_temp_1, amp);

	//MathOperationT::PrintVector(FBC_num, FBC_Values, "FBC_Values");
	//MathOperationT::PrintVector(FBC_num, FBC_Values_temp_1, "FBC_Values_temp_1");


	double* load_1=new double[fNumDofs_1];
	MathOperationT::VecSet(fNumDofs_1, load_1,0);

	for (int i=0; i<FBC_num; i++)
		load_1[FBC_DOFs[i]]=FBC_Values_temp_1[i];

	double* load_interior_1=new double[fNumInteriorDofs_1];
	MathOperationT::VecExtraction(fNumDofs_1, load_1, fNumInteriorDofs_1, fExtraction_Gr_col_1, load_interior_1);

	//MathOperationT::PrintVector(fNumInteriorDofs_1, load_interior_1, "load_interior_1");


	if(fRank==0)
	{
		fIerr=VecSetValues(fLoad, fNumInteriorDofs_1, fAssemble_Gr_col_1, load_interior_1, INSERT_VALUES);
	}
	fIerr=VecAssemblyBegin(fLoad);
	fIerr=VecAssemblyEnd(fLoad);



	/*********************************************************************
	 * Add region 2
	 *********************************************************************/
	fModel_2->BoundaryCondition_Global(UBC_num, &UBC_DOFs, &UBC_Values,
			                          FBC_num, &FBC_DOFs, &FBC_Values);

	double* FBC_Values_temp_2=new double[FBC_num];
	MathOperationT::MemCopy(FBC_num, FBC_Values, FBC_Values_temp_2);


	MathOperationT::VecScale(FBC_num, FBC_Values_temp_2, amp);

	//MathOperationT::PrintVector(FBC_num, FBC_Values, "FBC_Values");
	//MathOperationT::PrintVector(FBC_num, FBC_Values_temp_2, "FBC_Values_temp_1");

	double* load_2=new double[fNumDofs_2];
	MathOperationT::VecSet(fNumDofs_2, load_2,0);

	for (int i=0; i<FBC_num; i++)
		load_2[FBC_DOFs[i]]=FBC_Values_temp_2[i];

	double* load_interior_2=new double[fNumInteriorDofs_2];
	MathOperationT::VecExtraction(fNumDofs_2, load_2, fNumInteriorDofs_2, fExtraction_Gr_col_2, load_interior_2);
	if(fRank==0)
	{
		fIerr=VecSetValues(fLoad, fNumInteriorDofs_2, fAssemble_Gr_col_2, load_interior_2, INSERT_VALUES);
	}
	fIerr=VecAssemblyBegin(fLoad);
	fIerr=VecAssemblyEnd(fLoad);


	//VecView(fLoad,PETSC_VIEWER_STDOUT_WORLD);

	delete[] FBC_Values_temp_1;
	delete[] FBC_Values_temp_2;
	delete[] load_1;
	delete[] load_2;
	delete[] load_interior_1;
	delete[] load_interior_2;
}



void SolverMRT::Solve()
{
	Solve_Direct();
	//Solve_Averaging();
}


void SolverMRT::Solve_Direct()
{
	//MatView(fG_DB[0], PETSC_VIEWER_STDOUT_WORLD);
	//MatView(fH_DB[0], PETSC_VIEWER_STDOUT_WORLD);

	/********
	 * Form LHS
	 ********/
	fIerr=MatCopy(fG_DB[0], fG_Compute, DIFFERENT_NONZERO_PATTERN);
	fIerr=MatAssemblyBegin(fG_Compute, MAT_FINAL_ASSEMBLY);
	fIerr=MatAssemblyEnd(fG_Compute, MAT_FINAL_ASSEMBLY);

	fIerr=MatCopy(fH_DB[0], fH_Compute, DIFFERENT_NONZERO_PATTERN);
	fIerr=MatAssemblyBegin(fH_Compute, MAT_FINAL_ASSEMBLY);
	fIerr=MatAssemblyEnd(fH_Compute, MAT_FINAL_ASSEMBLY);

	/********
	 * Form RHS
	 ********/
	fIerr=MatMult(fG_Compute, fLoad, fRHS);

    //fIerr=VecView(fRHS, PETSC_VIEWER_STDOUT_WORLD);


	//*** add terms by loop
	Mat G_temp;
	Mat H_temp;
	fIerr=MatConvert(fG_DB[0], MATSAME, MAT_INITIAL_MATRIX, &G_temp);
	fIerr=MatConvert(fH_DB[0], MATSAME, MAT_INITIAL_MATRIX, &H_temp);

	for(int f_i=1; f_i<=min(fCurrStep, fMaxStep); f_i++)
	{
		fIerr=MatZeroEntries(G_temp);
		fIerr=MatAXPY(G_temp, 1, fG_DB[f_i],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatMultAdd(G_temp, fTrac_DB_Global[fCurrStep-f_i],fRHS, fRHS);

		fIerr=MatZeroEntries(H_temp);
		fIerr=MatAXPY(H_temp, 1, fH_DB[f_i],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatScale(H_temp, -1);
		fIerr=MatMultAdd(H_temp, fDis_DB_Global[fCurrStep-f_i],fRHS, fRHS);
	}
	fIerr=MatDestroy(&G_temp);
	fIerr=MatDestroy(&H_temp);


	cout << "current step " << fCurrStep << endl;

    	Vec result_rough;
	fIerr=VecCreateMPI(fComm, PETSC_DECIDE, fNumRow_H, &result_rough);

	fIerr=KSPSetOperators(fKsp, fH_Compute, fH_Compute);
	fIerr=KSPSolve(fKsp, fRHS, result_rough);

	fIerr=VecCopy(result_rough, fDis_DB_Global[fCurrStep]);

	fIerr=VecCopy(fLoad, fTrac_DB_Global[fCurrStep]);

	ExtractField();

	fIerr=VecDestroy(&result_rough);


    //fIerr=VecView(fDis_DB[fCurrStep], PETSC_VIEWER_STDOUT_WORLD);
	//fIerr=MatMult(fH_Compute, result_rough, fRHS);
}

void SolverMRT::Solve_Averaging()
{
    int solve_step=fCurrStep-1;

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

		/********
		 * Form RHS
		 ********/
		fIerr=MatMult(fG_Compute, fLoad, fRHS);
		fIerr=MatMultAdd(fG_DB_bgn[1], fTrac_DB_Global[0], fRHS, fRHS);
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


		//MatView(fG_Compute, PETSC_VIEWER_STDOUT_WORLD);

		/********
		 * Form RHS
		 ********/
		fIerr=MatMult(fG_Compute, fLoad, fRHS);

		//MatView(fG_Compute, PETSC_VIEWER_STDOUT_WORLD);
		//VecView(fLoad, PETSC_VIEWER_STDOUT_WORLD);


		//**add single terms
		Mat H_temp;
		Mat G_temp;
		fIerr=MatConvert(fG_DB[2], MATSAME, MAT_INITIAL_MATRIX, &G_temp);
		fIerr=MatAXPY(G_temp, 2, fG_DB[1],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatMultAdd(G_temp, fTrac_DB_Global[solve_step-1],fRHS, fRHS);

		//MatView(matrix_temp, PETSC_VIEWER_STDOUT_WORLD);
		//VecView(fRHS, PETSC_VIEWER_STDOUT_WORLD);

		fIerr=MatConvert(fH_DB[2], MATSAME, MAT_INITIAL_MATRIX, &H_temp);
		fIerr=MatAXPY(H_temp, 2, fH_DB[1],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatScale(H_temp, -1);
		fIerr=MatMultAdd(H_temp, fDis_DB_Global[solve_step-1],fRHS, fRHS);


		//*** add terms by loop
		for(int f_i=3; f_i<=min(solve_step, fMaxStep); f_i++)
		{
			fIerr=MatZeroEntries(G_temp);
			fIerr=MatAXPY(G_temp, 1, fG_DB[f_i],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatMultAdd(G_temp, fTrac_DB_Global[solve_step+1-f_i],fRHS, fRHS);

			fIerr=MatZeroEntries(H_temp);
			fIerr=MatAXPY(H_temp, 1, fH_DB[f_i],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatScale(H_temp, -1);
			fIerr=MatMultAdd(H_temp, fDis_DB_Global[solve_step+1-f_i],fRHS, fRHS);
		}

		for(int f_i=2; f_i<=min(solve_step-1, fMaxStep); f_i++)
		{
			fIerr=MatZeroEntries(G_temp);
			fIerr=MatAXPY(G_temp, 2, fG_DB[f_i],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatMultAdd(G_temp, fTrac_DB_Global[solve_step-f_i],fRHS, fRHS);

			fIerr=MatZeroEntries(H_temp);
			fIerr=MatAXPY(H_temp, 2, fH_DB[f_i],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatScale(H_temp, -1);
			fIerr=MatMultAdd(H_temp, fDis_DB_Global[solve_step-f_i],fRHS, fRHS);
		}

		for(int f_i=1; f_i<=min(solve_step-2, fMaxStep); f_i++)
		{
			fIerr=MatZeroEntries(G_temp);
			fIerr=MatAXPY(G_temp, 1, fG_DB[f_i],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatMultAdd(G_temp, fTrac_DB_Global[solve_step-1-f_i],fRHS, fRHS);

			fIerr=MatZeroEntries(H_temp);
			fIerr=MatAXPY(H_temp, 1, fH_DB[f_i],DIFFERENT_NONZERO_PATTERN);
			fIerr=MatScale(H_temp, -1);
			fIerr=MatMultAdd(H_temp, fDis_DB_Global[solve_step-1-f_i],fRHS, fRHS);
		}

		//** special treatment to the first step
		fIerr=MatZeroEntries(G_temp);
		fIerr=MatAXPY(G_temp, 1, fG_DB_bgn[0],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAXPY(G_temp, 2, fG_DB_bgn[1],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAXPY(G_temp, 1, fG_DB_bgn[2],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatMultAdd(G_temp, fTrac_DB_Global[0],fRHS, fRHS);

		fIerr=MatZeroEntries(H_temp);
		fIerr=MatAXPY(H_temp, 1, fH_DB_bgn[0],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAXPY(H_temp, 2, fH_DB_bgn[1],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatAXPY(H_temp, 1, fH_DB_bgn[2],DIFFERENT_NONZERO_PATTERN);
		fIerr=MatScale(H_temp, -1);
		fIerr=MatMultAdd(H_temp, fDis_DB_Global[0],fRHS, fRHS);

		fIerr=MatDestroy(&G_temp);
		fIerr=MatDestroy(&H_temp);
	}

	cout << "current step " << solve_step << endl;

	Vec result_rough;
	fIerr=VecCreateMPI(fComm, PETSC_DECIDE, fNumRow_H, &result_rough);

	fIerr=KSPSetOperators(fKsp, fH_Compute, fH_Compute);
	fIerr=KSPSolve(fKsp, fRHS, result_rough);

    fIerr=VecCopy(result_rough, fDis_DB_Global[solve_step]);

	fIerr=VecCopy(fLoad, fTrac_DB_Global[solve_step]);

	ExtractField();

	fIerr=VecDestroy(&result_rough);


	//fIerr=VecView(fDis_DB[solve_step], PETSC_VIEWER_STDOUT_WORLD);
	//fIerr=MatMult(fH_Compute, result_rough, fRHS);
}


void SolverMRT::ExtractField()
{
	Vec dis_collection;
	VecScatter ctx;

	VecScatterCreateToAll(fDis_DB_Global[fCurrStep-1],&ctx,&dis_collection);
	VecScatterBegin(ctx,fDis_DB_Global[fCurrStep-1],dis_collection,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(ctx,fDis_DB_Global[fCurrStep-1],dis_collection,INSERT_VALUES,SCATTER_FORWARD);

	double* temp;
	fIerr=VecGetArray(dis_collection, &temp);

	double* dis_1=new double[fNumDofs_1];
	double* dis_2=new double[fNumDofs_2];
	MathOperationT::VecExtraction(fNumRow_H, temp, fNumDofs_1, fExtraction_Dis_1, dis_1);
	MathOperationT::VecExtraction(fNumRow_H, temp, fNumDofs_2, fExtraction_Dis_2, dis_2);
	fIerr=VecRestoreArray(dis_collection, &temp);


	if(fRank==0)
	{
		//MathOperationT::PrintVector(fNumDofs_1, fExtraction_Dis_1, "fExtraction_Dis_1");
		//MathOperationT::PrintVector(fNumDofs_2, fExtraction_Dis_2, "fExtraction_Dis_2");

		//MathOperationT::PrintVector(fNumDofs_1, fSet_Dis_1, "fSet_Dis_1");
		//MathOperationT::PrintVector(fNumDofs_2, fSet_Dis_2, "fSet_Dis_2");

		VecSetValues(fDis_DB_1[fCurrStep-1], fNumDofs_1, fSet_Dis_1, dis_1, INSERT_VALUES);
		VecSetValues(fDis_DB_2[fCurrStep-1], fNumDofs_2, fSet_Dis_2, dis_2, INSERT_VALUES);
	}


	fIerr=VecAssemblyBegin(fDis_DB_1[fCurrStep-1]);
	fIerr=VecAssemblyEnd(fDis_DB_1[fCurrStep-1]);

	fIerr=VecAssemblyBegin(fDis_DB_2[fCurrStep-1]);
	fIerr=VecAssemblyEnd(fDis_DB_2[fCurrStep-1]);

	//VecView(fDis_DB_Global[fCurrStep-1],PETSC_VIEWER_STDOUT_WORLD);
	//VecView(fDis_DB_1[fCurrStep-1],PETSC_VIEWER_STDOUT_WORLD);
	//VecView(fDis_DB_2[fCurrStep-1],PETSC_VIEWER_STDOUT_WORLD);


	VecDestroy(&dis_collection);
	VecScatterDestroy(&ctx);

	delete[] dis_1;
	delete[] dis_2;
}



void SolverMRT::GetDofLowHigh(int* low_1, int* high_1, int* low_2, int* high_2)
{
	*low_1=fDofLow_1;
	*high_1=fDofHigh_1;

	*low_2=fDofLow_2;
	*high_2=fDofHigh_2;
}


void SolverMRT::SetCurrStep(int curr_step)
{
	fCurrStep=curr_step;
}

void SolverMRT::SetMaxStep(int MaxStep)
{
	fMaxStep=MaxStep;
}


void SolverMRT::RecordResult()
{
	ofstream myfile;

	Vec dis_collection_1;
	VecScatter ctx_1;

	VecScatterCreateToAll(fDis_DB_1[fCurrStep-1],&ctx_1,&dis_collection_1);
	VecScatterBegin(ctx_1,fDis_DB_1[fCurrStep-1],dis_collection_1,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(ctx_1,fDis_DB_1[fCurrStep-1],dis_collection_1,INSERT_VALUES,SCATTER_FORWARD);

	Vec dis_collection_2;
	VecScatter ctx_2;

	VecScatterCreateToAll(fDis_DB_2[fCurrStep-1],&ctx_2,&dis_collection_2);
	VecScatterBegin(ctx_2,fDis_DB_2[fCurrStep-1],dis_collection_2,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(ctx_2,fDis_DB_2[fCurrStep-1],dis_collection_2,INSERT_VALUES,SCATTER_FORWARD);

	if(fRank==0)
	{
		if (fCurrStep==1)
			myfile.open("displacement.txt", fstream::out|fstream::trunc);
		else
			myfile.open("displacement.txt", fstream::out|fstream::app);


		double* temp_1, *temp_2;

		fIerr=VecGetArray(dis_collection_1, &temp_1);
		fIerr=VecGetArray(dis_collection_2, &temp_2);

		myfile<<setprecision(6);

		/**for coarse mesh **/
		myfile<< setw(10) << fDt*(fCurrStep-1) << setw(15) << temp_2[46*3-1] << endl;


		/**for buried cavity Jiang **/
		//myfile<< setw(10) << fDt*(fCurrStep-1) << setw(15) << temp_2[2602*3-1] <<
		//		 setw(15) << temp_1[1682*3-1] << setw(15) << temp_1[1831*3-1] << endl;


		fIerr=VecRestoreArray(dis_collection_1, &temp_1);
		fIerr=VecRestoreArray(dis_collection_2, &temp_2);

		myfile.close();
	}

	VecScatterDestroy(&ctx_1);
	VecDestroy(&dis_collection_1);
	VecScatterDestroy(&ctx_2);
	VecDestroy(&dis_collection_2);

	OutPutVTK("BEM_1st_", fDis_DB_1, fModel_1);
	OutPutVTK("BEM_2nd_", fDis_DB_2, fModel_2);
}




void SolverMRT::OutPutVTK(std::string filename, Vec* Dis_DB, ModelManagerT* Model)
{
	ofstream myfile;

	Vec dis_collection;
	VecScatter ctx;

	VecScatterCreateToAll(Dis_DB[fCurrStep-1],&ctx,&dis_collection);
	VecScatterBegin(ctx,Dis_DB[fCurrStep-1],dis_collection,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(ctx,Dis_DB[fCurrStep-1],dis_collection,INSERT_VALUES,SCATTER_FORWARD);

	int nnd=Model->GetNND();
	int nel=Model->GetNEL();

	const int* IEN=Model->GetIEN();
	const double* Coords=Model->GetCoords();


	if(fRank==0)
	{
		ostringstream oo;
		oo<<fCurrStep;
		filename+=oo.str();

		filename+=".vtk";

		const char* name=filename.c_str();

		myfile.open(name, fstream::out|fstream::trunc);

		/* write the Head */
		myfile<<"# vtk DataFile Version 3.0 " << endl;
		myfile<<"Boundary element result"<<endl;
		myfile<<"ASCII"<<endl;
		myfile<<"DATASET UNSTRUCTURED_GRID " <<endl;
		myfile<<"POINTS " << nnd << " float"<<endl;

		myfile<<setprecision(6);

		for(int ni=0; ni<nnd; ni++)
		{
			myfile<< setw(15)<<Coords[3*ni] <<setw(15) <<Coords[3*ni+1] <<setw(15) <<Coords[3*ni+2];
			myfile<<endl;
		}

		myfile<<"CELLS  " << nel << " " << nel*5 <<endl;

		for(int ei=0; ei<nel; ei++)
		{
			myfile<< setw(5)<<4 << setw(5)<< IEN[4*ei] << setw(5) <<IEN[4*ei+1] << setw(5) <<IEN[4*ei+2]<< setw(5) <<IEN[4*ei+3];
			myfile<<endl;
		}

		myfile<<"CELL_TYPES "<<nel<<endl;
		for(int ei=0; ei<nel; ei++)
		{
			myfile<<9<<endl;
		}


		double* temp;
		fIerr=VecGetArray(dis_collection, &temp);

		myfile<<"POINT_DATA " << nnd << endl;
		myfile<<"VECTORS " << "Displacement " << "double"<<endl;
		for(int ni=0; ni<nnd; ni++)
		{
     		myfile<< setw(15)<< temp[3*ni] <<setw(15) <<temp[3*ni+1] <<setw(15) <<temp[3*ni+2];
			myfile<<endl;
		}

		fIerr=VecRestoreArray(dis_collection, &temp);


		myfile.close();

	}

	VecScatterDestroy(&ctx);
	VecDestroy(&dis_collection);
}



