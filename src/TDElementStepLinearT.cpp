/*
 * TDElementStepLinearT.cpp
 *
 *  Created on: Sep 21, 2015
 *      Author: xiaoyong
 */

#include "TDElementStepLinearT.h"

#include "MathOperationT.h"
#include "iostream"

using namespace std;
using namespace TD_BEM;
using namespace GreenFunction;

TDElementStepLinearT::TDElementStepLinearT(int type):
		TDElementT(type)
{
	//do nothing
}

TDElementStepLinearT::~TDElementStepLinearT()
{
}


void TDElementStepLinearT::FormGH_Regular()
{

	int ennd=fShape->GetENND();
	int sd=3;
	int G_length=ennd*sd*sd;
	int C_length=sd*sd;

	MathOperationT::VecSet(G_length, fGE_1, 0);
	MathOperationT::VecSet(G_length, fGE_2, 0);

	MathOperationT::VecSet(G_length, fHE_1, 0);
	MathOperationT::VecSet(G_length, fHE_2, 0);

	MathOperationT::VecSet(C_length, fC_1, 0);
	MathOperationT::VecSet(C_length, fC_2, 0);

	double* NMatrix=new double[3*ennd*3];
	double* direction=new double[3];
	double* matrix_temp=new double[3*(ennd*3)];

	//MathOperationT::VecSet(3*ennd*3, matrix_temp,0);
	double* DT_Kelvin_T=new double[3*3];

	for(int q_1=0; q_1<fNQ; q_1++)
	{
		for (int q_2=0; q_2<fNQ; q_2++)
		{
			fShape->Evaluate(q_1, q_2);
			const double* N=fShape->GetShapeValues();
			const double* x=fShape->GetX();
			const double* n=fShape->GetOutNormal();
			const double J=fShape->GetJacobian();
			const double w=fShape->GetWeight();

			MathOperationT::VecSet(3*ennd*3, NMatrix, 0);
			for(int a=0; a<ennd; a++)
			{
				for(int d=0; d<sd; d++)
				{
					NMatrix[d*3*ennd+d+3*a]=N[a];
				}
			}

			//MathOperationT::PrintMatrix(3, 3*ennd, NMatrix,"NMatrix");

			MathOperationT::VecPlus(3, 1, x, -1, fCollocation, direction);

			fGreenFunction->SetSpatial(direction, n);

			fGreenFunction->Compute(4);

			double* DU_Stokes_1;
			double* DU_Stokes_2;
			double* DT_Stokes_1;
			double* DT_Stokes_2;

			fGreenFunction->GetStokesLinear(&DU_Stokes_1, &DU_Stokes_2, &DT_Stokes_1, &DT_Stokes_2);
			//MathOperationT::PrintMatrix(3, DU_Stokes_1, "DU_Stokes_1");
			//MathOperationT::PrintMatrix(3, DU_Stokes_2, "DU_Stokes_2");
			//MathOperationT::PrintMatrix(3, DT_Stokes_1, "DT_Stokes_1");
			//MathOperationT::PrintMatrix(3, DT_Stokes_2, "DT_Stokes_2");

			MathOperationT::multATB(3,3,DU_Stokes_1, 3, 3*ennd, NMatrix, matrix_temp);
			MathOperationT::VecPlus(3*ennd*3, 1, fGE_1, w*J, matrix_temp, fGE_1);

			MathOperationT::multATB(3,3,DU_Stokes_2, 3, 3*ennd, NMatrix, matrix_temp);
			MathOperationT::VecPlus(3*ennd*3, 1, fGE_2, w*J, matrix_temp, fGE_2);

			MathOperationT::multATB(3,3,DT_Stokes_1, 3, 3*ennd, NMatrix, matrix_temp);
			MathOperationT::VecPlus(3*ennd*3, 1, fHE_1, w*J, matrix_temp, fHE_1);

			MathOperationT::multATB(3,3,DT_Stokes_2, 3, 3*ennd, NMatrix, matrix_temp);
			MathOperationT::VecPlus(3*ennd*3, 1, fHE_2, w*J, matrix_temp, fHE_2);

			if (fT==fT_1)
			{
				const double* DT_Kelvin=fGreenFunction->GetKelvinLinear();
				MathOperationT::transMatrix(3,3, DT_Kelvin, DT_Kelvin_T);
				MathOperationT::VecPlus(3*3, 1, fC_1, -w*J, DT_Kelvin_T, fC_1);
				MathOperationT::VecPlus(3*3, 1, fC_2, -w*J, DT_Kelvin_T, fC_2);
			}

			//MathOperationT::PrintMatrix(3, ennd*3, fGE_1, "G_1");
			//MathOperationT::PrintMatrix(3, ennd*3, fGE_2, "G_2");
			//MathOperationT::PrintMatrix(3, ennd*3, fHE_1, "H_1");
			//MathOperationT::PrintMatrix(3, ennd*3, fHE_2, "H_2");
			//MathOperationT::PrintMatrix(3, fC_1, "C_1");
			//MathOperationT::PrintMatrix(3, fC_2, "C_2");

		}
	}

	delete[] NMatrix;
	delete[] direction;
	delete[] matrix_temp;
	delete[] DT_Kelvin_T;
}


void TDElementStepLinearT::FormGH_Singular(int singular)
{
	//cout<< "FormGH_Singular_TL is called\n";
	int ennd=fShape->GetENND();
	int sd=3;
	int G_length=ennd*sd*sd;
	int C_length=sd*sd;

	MathOperationT::VecSet(G_length, fGE_1, 0);
	MathOperationT::VecSet(G_length, fGE_2, 0);

	MathOperationT::VecSet(G_length, fHE_1, 0);
	MathOperationT::VecSet(G_length, fHE_2, 0);

	MathOperationT::VecSet(C_length, fC_1, 0);
	MathOperationT::VecSet(C_length, fC_2, 0);

	if (fType==1)
	{
		FormGH_Singular_1(singular, fShape, fGE_1, fGE_2, fHE_1, fHE_2, fC_1, fC_2);
	}
	else if(fType==2)
	{
		SubElement(singular);

		fSubShape->SetNQ(fNQ);

		int ennd=fSubShape->GetENND();
		int sd=3;
		int G_length=ennd*sd*sd;

		//allocate space for sub element matrices;
		double* GE_1=new double[G_length];
		double* GE_2=new double[G_length];
		double* HE_1=new double[G_length];
		double* HE_2=new double[G_length];
        double* C_1=new double[sd*sd];
        double* C_2=new double[sd*sd];
		/* loop over sub elements */
		for(int se_i=0; se_i<fSubNum; se_i++)
		{
			fSubShape->SetCoords(fSubNodes[se_i]);
			FormGH_Singular_1(fSubCID[se_i], fSubShape, GE_1, GE_2, HE_1, HE_2, C_1, C_2);

			//MathOperationT::PrintMatrix(sd, ennd*sd, GE_1, "GE_1");

			AssembleSubElement(fSubIDs[se_i], GE_1, GE_2, HE_1, HE_2);
		}

		delete[] GE_1;
		delete[] GE_2;
		delete[] HE_1;
		delete[] HE_2;
        delete[] C_1;
        delete[] C_2;

	}
	else
	{
		cout<<"TDElementT::FormGH_Singular_TL, unsupported element type\n";
	}

}

void TDElementStepLinearT::SubElement(int singular)
{

	if(fType==2)
	{
		switch (singular)
		{
		case 0:
			fSubNum=1;
			fSubType=1;
			fSubCID[0]=singular;

			fSubIDs[0][0]=0;
			fSubIDs[0][1]=1;
			fSubIDs[0][2]=2;
			fSubIDs[0][3]=3;

			MathOperationT::MemCopy(12, fCoords, fSubNodes[0]);
			break;

		case 1:
			fSubNum=1;
			fSubType=1;

			fSubCID[0]=singular;

			fSubIDs[0][0]=0;
			fSubIDs[0][1]=1;
			fSubIDs[0][2]=2;
			fSubIDs[0][3]=3;

			MathOperationT::MemCopy(12, fCoords, fSubNodes[0]);
			break;

		case 2:
			fSubNum=1;
			fSubType=1;

			fSubCID[0]=singular;

			fSubIDs[0][0]=0;
			fSubIDs[0][1]=1;
			fSubIDs[0][2]=2;
			fSubIDs[0][3]=3;

			MathOperationT::MemCopy(12, fCoords, fSubNodes[0]);
			break;

		case 3:
			fSubNum=1;
			fSubType=1;

			fSubCID[0]=singular;

			fSubIDs[0][0]=0;
			fSubIDs[0][1]=1;
			fSubIDs[0][2]=2;
			fSubIDs[0][3]=3;

			MathOperationT::MemCopy(12, fCoords, fSubNodes[0]);
			break;

		case 4:
			fSubNum=2;
			fSubType=1;

			fSubCID[0]=0;
			fSubCID[1]=0;

			fSubIDs[0][0]=4;
			fSubIDs[0][1]=6;
			fSubIDs[0][2]=3;
			fSubIDs[0][3]=0;

			fSubIDs[1][0]=4;
			fSubIDs[1][1]=1;
			fSubIDs[1][2]=2;
			fSubIDs[1][3]=6;

			for(int e=0; e<fSubNum; e++)
			{
				for(int i=0; i<4; i++)
				{
					int id=fSubIDs[e][i];
					for(int d=0; d<3; d++)
					{
						fSubNodes[e][3*i+d]=fCoords[3*id+d];
					}
				}
			}

			break;

		case 5:
			fSubNum=2;
			fSubType=1;

			fSubCID[0]=0;
			fSubCID[1]=0;

			fSubIDs[0][0]=5;
			fSubIDs[0][1]=2;
			fSubIDs[0][2]=3;
			fSubIDs[0][3]=7;

			fSubIDs[1][0]=5;
			fSubIDs[1][1]=7;
			fSubIDs[1][2]=0;
			fSubIDs[1][3]=1;

			for(int e=0; e<fSubNum; e++)
			{
				for(int i=0; i<4; i++)
				{
					int id=fSubIDs[e][i];
					for(int d=0; d<3; d++)
					{
						fSubNodes[e][3*i+d]=fCoords[3*id+d];
					}
				}
			}

			break;
		case 6:
			fSubNum=2;
			fSubType=1;

			fSubCID[0]=0;
			fSubCID[1]=0;

			fSubIDs[0][0]=6;
			fSubIDs[0][1]=3;
			fSubIDs[0][2]=0;
			fSubIDs[0][3]=4;

			fSubIDs[1][0]=6;
			fSubIDs[1][1]=4;
			fSubIDs[1][2]=1;
			fSubIDs[1][3]=2;

			for(int e=0; e<fSubNum; e++)
			{
				for(int i=0; i<4; i++)
				{
					int id=fSubIDs[e][i];
					for(int d=0; d<3; d++)
					{
						fSubNodes[e][3*i+d]=fCoords[3*id+d];
					}
				}
			}

			break;
		case 7:
			fSubNum=2;
			fSubType=1;

			fSubCID[0]=0;
			fSubCID[1]=0;

			fSubIDs[0][0]=7;
			fSubIDs[0][1]=0;
			fSubIDs[0][2]=1;
			fSubIDs[0][3]=5;

			fSubIDs[1][0]=7;
			fSubIDs[1][1]=5;
			fSubIDs[1][2]=2;
			fSubIDs[1][3]=3;

			for(int e=0; e<fSubNum; e++)
			{
				for(int i=0; i<4; i++)
				{
					int id=fSubIDs[e][i];
					for(int d=0; d<3; d++)
					{
						fSubNodes[e][3*i+d]=fCoords[3*id+d];
					}
				}
			}

			break;

		default:
			throw "TDElementT::SubElement singular is out of range";
		} /*end switch*/
	}/* end if */
	else
	{
		throw "TDElementT::SubElement unsupported type of element";
	}


}


void TDElementStepLinearT::FormGH_Singular_1(int singular, ShapeFunctionT* Shape,
								   double* GE_1, double* GE_2, double* HE_1, double* HE_2, double* C_1, double* C_2)
{
	int ennd=Shape->GetENND();
	int sd=3;
	int G_length=ennd*sd*sd;
	int C_length=sd*sd;

	MathOperationT::VecSet(G_length, GE_1, 0);
	MathOperationT::VecSet(G_length, GE_2, 0);

	MathOperationT::VecSet(G_length, HE_1, 0);
	MathOperationT::VecSet(G_length, HE_2, 0);
    
    MathOperationT::VecSet(C_length, C_1, 0);
    MathOperationT::VecSet(C_length, C_2, 0);

	double eta_c[2];

	switch (singular){
	case 0:
		eta_c[0]=-1;
		eta_c[1]=-1;
		break;
	case 1:
		eta_c[0]=1;
		eta_c[1]=-1;
		break;
	case 2:
		eta_c[0]=1;
		eta_c[1]=1;
		break;
	case 3:
		eta_c[0]=-1;
		eta_c[1]=1;
		break;
	default:
		throw "TDElementT::FormGH_Singular_TL, collocation id is out of range ";
	}

	double* NMatrix=new double[3*ennd*3];
	double* direction=new double[3];
	double* matrix_temp=new double[3*(ennd*3)];

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	// integrate over triangle A
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	for(int q_1=0; q_1<fNQ; q_1++)
	{
		for (int q_2=0; q_2<fNQ; q_2++)
		{
			double eta_tri_1, eta_tri_2; //Shape coordinate of the Gaussian point
			double w; //weight of the current Gaussian point

			Shape->GetShapeCoords(q_1, q_2, eta_tri_1, eta_tri_2, w);

			double J_tri=0.5*(1-eta_tri_1);	//Jacobian of the triangle

			//Shape coordinates of the real element
			double eta_1=eta_c[0]*(1-0.5*(1-eta_tri_1)*(1+eta_c[0]*eta_c[1]*eta_tri_2));
			double eta_2=eta_c[1]*eta_tri_1;

			Shape->Evaluate(eta_1, eta_2);
			const double* N=Shape->GetShapeValues();
			const double* x=Shape->GetX();
			const double* n=Shape->GetOutNormal();
			const double J=Shape->GetJacobian();

			//form shape function matrix
			MathOperationT::VecSet(3*ennd*3, NMatrix, 0);
			for(int a=0; a<ennd; a++)
			{
				for(int d=0; d<sd; d++)
				{
					NMatrix[d*3*ennd+d+3*a]=N[a];
				}
			}

			direction[0]=x[0]-fCollocation[0];
			direction[1]=x[1]-fCollocation[1];
			direction[2]=x[2]-fCollocation[2];

			fGreenFunction->SetSpatial(direction, n);

			fGreenFunction->Compute(4);

			double* DU_Stokes_1;
			double* DU_Stokes_2;
			double* DT_Stokes_1;
			double* DT_Stokes_2;

			fGreenFunction->GetStokesLinear(&DU_Stokes_1, &DU_Stokes_2, &DT_Stokes_1, &DT_Stokes_2);

			//MathOperationT::PrintVector(3, direction, "direction");
			//MathOperationT::PrintVector(3, n, "n");
			//MathOperationT::PrintMatrix(3, DT_Kelvin, "DT_Kelvin");
			//MathOperationT::PrintMatrix(3, DU_Stokes_1, "fDU_Stokes_1");
			//MathOperationT::PrintMatrix(3, DT_Stokes_1, "fDT_Stokes_1");
			//MathOperationT::PrintMatrix(3, DT_Stokes_2, "fDT_Stokes_2");

			MathOperationT::multATB(3,3,DU_Stokes_1, 3, 3*ennd, NMatrix, matrix_temp);
			MathOperationT::VecPlus(3*ennd*3, 1, GE_1, w*J_tri*J, matrix_temp, GE_1);

			MathOperationT::multATB(3,3,DU_Stokes_2, 3, 3*ennd, NMatrix, matrix_temp);
			MathOperationT::VecPlus(3*ennd*3, 1, GE_2, w*J_tri*J, matrix_temp, GE_2);

			MathOperationT::multATB(3,3,DT_Stokes_1, 3, 3*ennd, NMatrix, matrix_temp);
			MathOperationT::VecPlus(3*ennd*3, 1, HE_1, w*J_tri*J, matrix_temp, HE_1);

			MathOperationT::multATB(3,3,DT_Stokes_2, 3, 3*ennd, NMatrix, matrix_temp);
			MathOperationT::VecPlus(3*ennd*3, 1, HE_2, w*J_tri*J, matrix_temp, HE_2);

			if (fT==fT_1)
			{
				const double* DT_Kelvin=fGreenFunction->GetKelvinLinear();
		double DT_Kelvin_T[9];
                MathOperationT::transMatrix(3,3, DT_Kelvin, DT_Kelvin_T);
                MathOperationT::VecPlus(3*3, 1, C_1, -w*J*J_tri, DT_Kelvin_T, C_1);
                MathOperationT::VecPlus(3*3, 1, C_2, -w*J*J_tri, DT_Kelvin_T, C_2);
			}

		}/* loop over q_1 */
	}/* loop over q_2*/

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	// integrate over triangle B
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	for(int q_1=0; q_1<fNQ; q_1++)
	{
		for (int q_2=0; q_2<fNQ; q_2++)
		{
			double eta_tri_1, eta_tri_2; //Shape coordinate of the Gaussian point
			double w; //weight of the current Gaussian point

			Shape->GetShapeCoords(q_1, q_2, eta_tri_1, eta_tri_2, w);

			//Jacobian of the triangle
			double J_tri=0.5*(1-eta_tri_1);

			//Shape coordinates of the real element
			double eta_1=eta_c[0]*eta_tri_1;
			double eta_2=eta_c[1]*(1-0.5*(1-eta_tri_1)*(1-eta_c[0]*eta_c[1]*eta_tri_2));

			Shape->Evaluate(eta_1, eta_2);
			const double* N=Shape->GetShapeValues();
			const double* x=Shape->GetX();
			const double* n=Shape->GetOutNormal();
			const double J=Shape->GetJacobian();

			//form shape function matrix
			MathOperationT::VecSet(3*ennd*3, NMatrix, 0);
			for(int a=0; a<ennd; a++)
			{
				for(int d=0; d<sd; d++)
				{
					NMatrix[d*3*ennd+d+3*a]=N[a];
				}
			}

			//MathOperationT::PrintMatrix(3, 3*ennd, NMatrix, "NMatrix");
			direction[0]=x[0]-fCollocation[0];
			direction[1]=x[1]-fCollocation[1];
			direction[2]=x[2]-fCollocation[2];

			fGreenFunction->SetSpatial(direction, n);

			fGreenFunction->Compute(4);

			const double* DT_Kelvin=fGreenFunction->GetKelvinLinear();
			double* DU_Stokes_1;
			double* DU_Stokes_2;
			double* DT_Stokes_1;
			double* DT_Stokes_2;

			fGreenFunction->GetStokesLinear(&DU_Stokes_1, &DU_Stokes_2, &DT_Stokes_1, &DT_Stokes_2);

			MathOperationT::multATB(3,3,DU_Stokes_1, 3, 3*ennd, NMatrix, matrix_temp);
			MathOperationT::VecPlus(3*ennd*3, 1, GE_1, w*J_tri*J, matrix_temp, GE_1);

			MathOperationT::multATB(3,3,DU_Stokes_2, 3, 3*ennd, NMatrix, matrix_temp);
			MathOperationT::VecPlus(3*ennd*3, 1, GE_2, w*J_tri*J, matrix_temp, GE_2);

			MathOperationT::multATB(3,3,DT_Stokes_1, 3, 3*ennd, NMatrix, matrix_temp);
			MathOperationT::VecPlus(3*ennd*3, 1, HE_1, w*J_tri*J, matrix_temp, HE_1);

			MathOperationT::multATB(3,3,DT_Stokes_2, 3, 3*ennd, NMatrix, matrix_temp);
			MathOperationT::VecPlus(3*ennd*3, 1, HE_2, w*J_tri*J, matrix_temp, HE_2);

			if (fT==fT_1)
			{
				const double* DT_Kelvin=fGreenFunction->GetKelvinLinear();
		double DT_Kelvin_T[9];
                MathOperationT::transMatrix(3,3, DT_Kelvin, DT_Kelvin_T);
                MathOperationT::VecPlus(3*3, 1, C_1, -w*J*J_tri, DT_Kelvin_T, C_1);
                MathOperationT::VecPlus(3*3, 1, C_2, -w*J*J_tri, DT_Kelvin_T, C_2);
			}

		}/* loop over q_1 */
	}/* loop over q_2*/

	delete[] NMatrix;
	delete[] direction;
	delete[] matrix_temp;
}




void TDElementStepLinearT::AssembleSubElement(const int* ids, const double* GE_1, const double* GE_2,
													const double* HE_1, const double* HE_2)
{
//MathOperationT::PrintVector(4, ids, "ids");

	int ennd_father=8;
	int ennd_son=4;

	for(int ni=0; ni<4; ni++)
	{
		int id=ids[ni];
		for(int ri=0; ri<3; ri++)
		{
			for(int ci=0; ci<3; ci++)
			{
				fGE_1[ri*(ennd_father*3)+3*id+ci]+=GE_1[ri*(ennd_son*3)+3*ni+ci];
				fGE_2[ri*(ennd_father*3)+3*id+ci]+=GE_2[ri*(ennd_son*3)+3*ni+ci];
				fHE_1[ri*(ennd_father*3)+3*id+ci]+=HE_1[ri*(ennd_son*3)+3*ni+ci];
				fHE_2[ri*(ennd_father*3)+3*id+ci]+=HE_2[ri*(ennd_son*3)+3*ni+ci];
			}
		}
	}
}


