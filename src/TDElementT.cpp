/*
 * TDElementT.cpp
 *
 *  Created on: May 29, 2015
 *      Author: xiaoyong
 */

#include "TDElementT.h"
#include "MathOperationT.h"

#include <iostream>

using namespace TD_BEM;
using namespace GreenFunction;
using namespace std;

TDElementT::TDElementT(int type)
{
	fType=type;

	SetShape(fType); //initialize the shape function


	int sd=3;
	fC_1=new double[sd*sd];
	fC_2=new double[sd*sd];


	fGeomB=new double[sd*sd];

	fGreenFunction=new GreenFunctionT;

	fSubCID=new int[4];
	fSubNodes=new double*[4];
	fSubIDs=new int*[4];

	for(int i=0; i<4; i++)
	{
		fSubNodes[i]=new double[3*4];
		fSubIDs[i]=new int[4];
	}

	switch(fType)
	{
	case 1:
		//do nothing, no sub element is needed
		fSubShape=NULL;
		break;
	case 2:
		fSubShape=new ShapeFunctionT(1);
		break;
	default:
		throw "TDElementT::TDElementT, unsupported element type. There is no sub element for this element";
	}

	fCollocation=NULL;
	fCoords=NULL;


}

TDElementT::~TDElementT()
{

	delete[] fGE_1;
	delete[] fGE_2;


	delete[] fHE_1;
	delete[] fHE_2;


	delete fGreenFunction;
	delete fShape;
	if(fSubShape!=NULL) delete fSubShape;

	delete[] fC_1;
	delete[] fC_2;

	delete[] fGeomB;

	delete[] fSubCID;

	for(int i=0; i<4; i++)
	{
		delete[] fSubNodes[i];
		delete[] fSubIDs[i];
	}

	delete[] fSubNodes;

	delete[] fSubIDs;
}

void TDElementT::SetCollocation(double* collocation)
{
	fCollocation=collocation;
}

void TDElementT::SetCoords(double* coords)
{
	fCoords=coords;

	fShape->SetCoords(fCoords);
}

void TDElementT::SetTime(double t, double t_1)
{
	fT=t;
	fT_1=t_1;

	fGreenFunction->SetTime(t, t_1);
    
}

void TDElementT::SetNQ(int nq)
{
	fNQ=nq;
	fShape->SetNQ(nq);
}

void TDElementT::SetMaterial(double E, double nu, double density)
{
	fE=E;
	fNu=nu;
	fDensity=density;

	fGreenFunction->SetMaterial(E, nu, density);
}

void TDElementT::SetShape(int type)
{
	fType=type;

	fShape=new ShapeFunctionT(fType);

	int ennd=fShape->GetENND();
	int sd=3;
	int G_length=ennd*sd*sd;

	fGE_1=new double[G_length];
	fGE_2=new double[G_length];


	fHE_1=new double[G_length];
	fHE_2=new double[G_length];

}

void TDElementT::GetGHLinear(const double** GE_1, const double** GE_2, const double** HE_1, const double** HE_2,
				 const double** C_1, const double** C_2)
{
	*GE_1=fGE_1;
	*GE_2=fGE_2;
	*HE_1=fHE_1;
	*HE_2=fHE_2;
	*C_1=fC_1;
	*C_2=fC_2;
}

void TDElementT::GetGeomBMatrix(const double** GeomB)
{
	*GeomB=fGeomB;
}


void TDElementT::FormGeomBMatrix()
{
	int sd=3;
	int ennd=fShape->GetENND();

	MathOperationT::VecSet(sd*sd, fGeomB, 0);

	double* direction=new double[3];

	double* KelvinTraction=new double[3*3];

	//MathOperationT::PrintVector(3, fCollocation, "fCollocation");


	for(int q_1=0; q_1<fNQ; q_1++)
	{
		for (int q_2=0; q_2<fNQ; q_2++)
		{
			fShape->Evaluate(q_1, q_2);
			const double* x=fShape->GetX();
			const double* n=fShape->GetOutNormal();
			const double J=fShape->GetJacobian();
			const double w=fShape->GetWeight();


			MathOperationT::VecPlus(3, 1, x, -1, fCollocation, direction);


			//MathOperationT::PrintVector(3, x, "x");
			//MathOperationT::PrintVector(3, n, "n");
			//MathOperationT::PrintVector(3, direction, "direction");


			fGreenFunction->SetSpatial(direction, n);

			fGreenFunction->KelvinTraction();


			const double* T_Kelvin=fGreenFunction->GetKelvinTraction();
			MathOperationT::transMatrix(3, 3, T_Kelvin, KelvinTraction);

			//MathOperationT::PrintMatrix(3, KelvinTraction, "KelvinTraction");


			MathOperationT::VecPlus(3*3, 1, fGeomB, -w*J, KelvinTraction, fGeomB);

			//MathOperationT::PrintMatrix(3, 3, fGeomB, "fGeomB");


		}
	}

	//MathOperationT::PrintMatrix(3, 3, fGeomB, "fGeomB");


	delete[] direction;
	delete[] KelvinTraction;
}


GreenFunctionT* TDElementT::GetGreenFunction()
{
	return fGreenFunction;
}











