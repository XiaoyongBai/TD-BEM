/*
 * ModelManagerT.cpp
 *
 *  Created on: May 28, 2015
 *      Author: xiaoyong
 */

#include "ModelManagerT.h"
#include "MathOperationT.h"

#include <iostream>

using namespace std;
using namespace TD_BEM;
using namespace GreenFunction;


ModelManagerT::ModelManagerT()
{
	fFormulation=1;
	fDT=0;
	fNumStep=1;
	fTimeInterpolation=1;
	fNQ=6;

	fSD=3;
	fIfExt=0;
	fNND=0;
	fENND=0;
	fNEL=0;
	fET=1;
	fIEN=NULL;
	fCoords=NULL;



	fNND_Aux=0;
	fENND_Aux=0;
	fNEL_Aux=0;
	fET_Aux=1;
	fIEN_Aux=NULL;
	fCoords_Aux=NULL;


	fE=0;
	fNu=0;
	fDensity=0;

	fUBC_num_global=0;
	fUBC_DOFs_global=NULL;
	fUBC_Values_global=NULL;
	fFBC_num_global=0;
	fFBC_DOFs_global=NULL;
	fFBC_Values_global=NULL;

	fUBC_num_local=0;
	fUBC_DOFs_local=NULL;
	fUBC_Values_local=NULL;
	fFBC_num_local=0;
	fFBC_DOFs_local=NULL;
	fFBC_Values_local=NULL;

	fNumInterface=0;
	fInterface=NULL;
}

ModelManagerT::~ModelManagerT()
{
	if (fIEN) delete[] fIEN;
	if (fCoords) delete[] fCoords;

	if (fIEN_Aux) delete[] fIEN_Aux;
	if (fCoords_Aux) delete[] fCoords_Aux;

	if (fUBC_DOFs_global) delete[] fUBC_DOFs_global;
	if (fUBC_Values_global) delete[] fUBC_Values_global;

	if (fFBC_DOFs_global) delete[] fFBC_DOFs_global;
	if (fFBC_Values_global) delete[] fFBC_Values_global;

	if (fUBC_DOFs_local) delete[] fUBC_DOFs_local;
	if (fUBC_Values_local) delete[] fUBC_Values_local;

	if (fFBC_DOFs_local)delete[] fFBC_DOFs_local;
	if (fFBC_Values_local) delete[] fFBC_Values_local;

	if (fInterface) delete[] fInterface;
}


int ModelManagerT::GetFormulation(void)
{
	return fFormulation;
}

void ModelManagerT::SetFormulation(int formulation)
{
	fFormulation=formulation;
}

int ModelManagerT::GetNumStep(void)
{
	return fNumStep;
}

void ModelManagerT::SetNumStep(int numstep)
{
	fNumStep=numstep;
}

double ModelManagerT::GetDT(void)
{
	return fDT;
}
void ModelManagerT::SetDT(double dt)
{
	fDT=dt;
}

const int ModelManagerT::GetTimeInterpolation(void)
{
	return fTimeInterpolation;
}

void ModelManagerT::SetTimeInterpolation(int interpolation)
{
	fTimeInterpolation=interpolation;
}

const int ModelManagerT::GetNQ(void)
{
	return fNQ;
}

void ModelManagerT::SetNQ(int nq)
{
	fNQ=nq;
}


int ModelManagerT::GetSD(void)
{
	return fSD;
}

void ModelManagerT::SetSD(int sd)
{
	fSD=sd;
}

int ModelManagerT::GetIfExt(void)
{
	return fIfExt;
}

void ModelManagerT::SetIfExt(int if_ext)
{
	fIfExt=if_ext;
}


int ModelManagerT::GetNND(void)
{
	return fNND;
}

void ModelManagerT::SetNND(int nnd)
{
	fNND=nnd;
}

int ModelManagerT::GetElementType()
{
	return fET;
}

void ModelManagerT::SetElementType(int et)
{
	fET=et;

	switch (et){
	case 1:
		fENND=4;
		break;
	case 2:
		fENND=8;
		break;
	default:
		throw "ModelManagerT::SetElementType, unsupported element type";
		break;
	}
}

int ModelManagerT::GetENND()
{
	return fENND;
}

int ModelManagerT::GetNEL(void)
{
	return fNEL;
}

void ModelManagerT::SetNEL(int nel)
{
	fNEL=nel;
}

const int* ModelManagerT::GetIEN(void)
{
	return fIEN;
}

void ModelManagerT::SetIEN(int* IEN)
{
	fIEN=new int[fNEL*fENND];

	MathOperationT::MemCopy(fNEL*fENND, IEN, fIEN);
}

const double* ModelManagerT::GetCoords(void)
{
	return fCoords;
}

void ModelManagerT::SetCoords(double* coords)
{
	fCoords=new double[fNND*fSD];

	MathOperationT::MemCopy(fNND*fSD, coords, fCoords);
}




int ModelManagerT::GetNND_Aux(void)
{
	return fNND_Aux;
}

void ModelManagerT::SetNND_Aux(int nnd)
{
	fNND_Aux=nnd;
}

int ModelManagerT::GetElementType_Aux()
{
	return fET_Aux;
}

void ModelManagerT::SetElementType_Aux(int et)
{
	fET_Aux=et;

	switch (et){
	case 1:
		fENND_Aux=4;
		break;
	case 2:
		fENND_Aux=8;
		break;
	default:
		throw "ModelManagerT::SetElementType_aux, unsupported element type";
		break;
	}
}

int ModelManagerT::GetENND_Aux()
{
	return fENND_Aux;
}

int ModelManagerT::GetNEL_Aux(void)
{
	return fNEL_Aux;
}

void ModelManagerT::SetNEL_Aux(int nel)
{
	fNEL_Aux=nel;
}

const int* ModelManagerT::GetIEN_Aux(void)
{
	return fIEN_Aux;
}

void ModelManagerT::SetIEN_Aux(int* IEN)
{
	fIEN_Aux=new int[fNEL_Aux*fENND_Aux];

	MathOperationT::MemCopy(fNEL_Aux*fENND_Aux, IEN, fIEN_Aux);
}

const double* ModelManagerT::GetCoords_Aux(void)
{
	return fCoords_Aux;
}

void ModelManagerT::SetCoords_Aux(double* coords)
{
	fCoords_Aux=new double[fNND_Aux*fSD];

	MathOperationT::MemCopy(fNND_Aux*fSD, coords, fCoords_Aux);
}






void ModelManagerT::GetMaterialConstants(double& E, double& nu, double& rho)
{
	E=fE;
	nu=fNu;
	rho=fDensity;
}

void ModelManagerT::SetMaterialConstants(double E, double nu, double rho)
{
	fE=E;
	fNu=nu;
	fDensity=rho;
}

void ModelManagerT::BoundaryCondition_Global(int& UBC_num, int** UBC_DOFs, double** UBC_Values,
						   	   	   	  int& FBC_num, int** FBC_DOFs, double** FBC_Values)
{
    UBC_num=fUBC_num_global;
    FBC_num=fFBC_num_global;
    *UBC_DOFs=fUBC_DOFs_global;
    *UBC_Values=fUBC_Values_global;
    *FBC_DOFs=fFBC_DOFs_global;
    *FBC_Values=fFBC_Values_global;
}

void ModelManagerT::SetBCs(int UBC_num, int* UBC_DOFs, double* UBC_values,
						   int FBC_num, int* FBC_DOFs, double* FBC_values)
{
	fUBC_num_global=UBC_num;
	fUBC_DOFs_global=new int[fUBC_num_global];
	fUBC_Values_global=new double[fUBC_num_global];

	MathOperationT::MemCopy(fUBC_num_global, UBC_DOFs, fUBC_DOFs_global);
	MathOperationT::MemCopy(fUBC_num_global, UBC_values, fUBC_Values_global);

	fFBC_num_global=FBC_num;
	fFBC_DOFs_global=new int[fFBC_num_global];
	fFBC_Values_global=new double[fFBC_num_global];

	MathOperationT::MemCopy(fFBC_num_global, FBC_DOFs, fFBC_DOFs_global);
	MathOperationT::MemCopy(fFBC_num_global, FBC_values, fFBC_Values_global);
}


void ModelManagerT::ReadInput(const char* name)
{
	fInput.ReadInput(name);

	SetFormulation(fInput.GetFormulation());
	SetNumStep(fInput.GetNumStep());
	SetDT(fInput.GetDT());
	SetTimeInterpolation(fInput.GetTimeInterpolation());
	SetNQ(fInput.GetNQ());

	SetSD(fInput.GetSD());
	SetIfExt(fInput.GetIfExt());
	SetNND(fInput.GetNND());
	SetElementType(fInput.GetElementType());
	SetNEL(fInput.GetNEL());

	SetIEN(fInput.GetIEN());

	SetCoords(fInput.GetCoords());



	if (fFormulation==3)
	{
		SetNND_Aux(fInput.GetNND_Aux());
		SetElementType_Aux(fInput.GetElementType_Aux());
		SetNEL_Aux(fInput.GetNEL_Aux());

		SetIEN_Aux(fInput.GetIEN_Aux());

		SetCoords_Aux(fInput.GetCoords_Aux());
	}


	double E, nu, rho;
	fInput.GetMaterialConstants(E, nu, rho);
	SetMaterialConstants(E, nu, rho);

	int UBC_num, FBC_num;
	int* UBC_DOFs, *FBC_DOFs;
	double* UBC_Values, *FBC_Values;
	fInput.BoundaryCondition(UBC_num, &UBC_DOFs, &UBC_Values,
						     FBC_num, &FBC_DOFs, &FBC_Values);

	SetBCs(UBC_num, UBC_DOFs, UBC_Values,
		   FBC_num, FBC_DOFs, FBC_Values);

	int NumInterface;
	int* interface;
	fInput.GetInterface(NumInterface, &interface);

	SetNumInterface(NumInterface);
	SetInterface(interface);

	//MathOperationT::PrintVector(UBC_num, UBC_Values, "model UBC_values");
	//MathOperationT::PrintVector(FBC_num, FBC_Values, "model FBC_values");

}


int ModelManagerT::GetNumInterface(void)
{
	return fNumInterface;
}

void ModelManagerT::SetNumInterface(int num_interface)
{
	fNumInterface=num_interface;
}

const int* ModelManagerT::GetInterface(void)
{
	return fInterface;
}


void ModelManagerT::SetInterface(int* Interface)
{
	fInterface=new int[fNumInterface];

	MathOperationT::MemCopy(fNumInterface, Interface, fInterface);
}



void ModelManagerT::Distributed_BC(int low, int high)
{
	int num_ubc_local=0;

	for(int i=0; i<fUBC_num_global; i++)
	{
		int dof=fUBC_DOFs_global[i];

		if (dof>=low && dof<=high) num_ubc_local++;
	}

	int* ubc_dof_temp=new int[num_ubc_local];
	double* ubc_values_temp=new double[num_ubc_local];

	int ui_temp=0;

	for(int i=0; i<fUBC_num_global; i++)
	{
		int dof=fUBC_DOFs_global[i];

		if (dof>=low && dof<=high)
		{

			ubc_dof_temp[ui_temp]=dof;
			ubc_values_temp[ui_temp]=fUBC_Values_global[i];

			ui_temp++;
		}
	}

	fUBC_num_local=num_ubc_local;
	fUBC_DOFs_local=new int[num_ubc_local];
	fUBC_Values_local=new double[num_ubc_local];

	MathOperationT::MemCopy(num_ubc_local, ubc_dof_temp, fUBC_DOFs_local);
	MathOperationT::MemCopy(num_ubc_local, ubc_values_temp, fUBC_Values_local);

	delete[] ubc_dof_temp;
	delete[] ubc_values_temp;


	int num_fbc_local=0;

	for(int i=0; i<fFBC_num_global; i++)
	{
		int dof=fFBC_DOFs_global[i];

		if (dof>=low && dof<=high) num_fbc_local++;
	}

	int* fbc_dof_temp=new int[num_fbc_local];
	double* fbc_values_temp=new double[num_fbc_local];

	int fi_temp=0;

	for(int i=0; i<fFBC_num_global; i++)
	{
		int dof=fFBC_DOFs_global[i];

		if (dof>=low && dof<=high)
		{

			fbc_dof_temp[fi_temp]=dof;
			fbc_values_temp[fi_temp]=fFBC_Values_global[i];

			fi_temp++;
		}
	}


	fFBC_num_local=num_fbc_local;
	fFBC_DOFs_local=new int[num_fbc_local];
	fFBC_Values_local=new double[num_fbc_local];

	MathOperationT::MemCopy(num_fbc_local, fbc_dof_temp, fFBC_DOFs_local);
	MathOperationT::MemCopy(num_fbc_local, fbc_values_temp, fFBC_Values_local);


	delete[] fbc_dof_temp;
	delete[] fbc_values_temp;

	/*cout << "in Model low is " << low << " high is " << high <<endl;

	cout << "UBC number =" <<fUBC_num_local<<endl;
	MathOperationT::PrintVector(fUBC_num_local, fUBC_DOFs_local, "model UBC_DOFs");
	MathOperationT::PrintVector(fUBC_num_local, fUBC_Values_local, "model UBC_Values");


	cout << "FBC number =" <<fFBC_num_local<<endl;
	MathOperationT::PrintVector(fFBC_num_local, fFBC_DOFs_local, "model FBC_DOFs");
	MathOperationT::PrintVector(fFBC_num_local, fFBC_Values_local, "model FBC_Values");*/

}


void ModelManagerT::BoundaryCondition_Local(int& UBC_num, int** UBC_DOFs, double** UBC_Values,
						   	   	   	  int& FBC_num, int** FBC_DOFs, double** FBC_Values)
{
    UBC_num=fUBC_num_local;
    FBC_num=fFBC_num_local;
    *UBC_DOFs=fUBC_DOFs_local;
    *UBC_Values=fUBC_Values_local;
    *FBC_DOFs=fFBC_DOFs_local;
    *FBC_Values=fFBC_Values_local;
}

