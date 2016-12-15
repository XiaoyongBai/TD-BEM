/*
 * ElastoDynamicsT.cpp
 *
 *  Created on: May 28, 2015
 *      Author: xiaoyong
 */

#include "ElastoDynamicsT.h"
#include "MathOperationT.h"
#include "TDElementT.h"
#include "TDElementHeavisideLinearT.h"
#include "TDElementRampHeaviLinearT.h"
#include "TDElementStepLinearT.h"

#include <iostream>

using namespace TD_BEM;
using namespace GreenFunction;
using namespace std;

ElastoDynamicsT::ElastoDynamicsT(ModelManagerT* Model)
{
	fModel=Model;

	int etype=fModel->GetElementType();

    fElement=new TDElementStepLinearT(etype);

	fElement->SetNQ(fModel->GetNQ());

	double E, nu, rho;
	fModel->GetMaterialConstants(E, nu, rho);
	fElement->SetMaterial(E, nu, rho);

	fG_1=NULL;
	fG_2=NULL;
	fH_1=NULL;
	fH_2=NULL;

	SetLowerUpper(0, fModel->GetNND()-1);

	cout << "ElastoDynamicsT is constructed\n";
}

ElastoDynamicsT::~ElastoDynamicsT()
{
	delete[] fG_1;
	delete[] fG_2;

	delete[] fH_1;
	delete[] fH_2;

	if(fElement)
	{
		delete fElement;
		fElement=NULL;
	}
}

void ElastoDynamicsT::SetLowerUpper(int Lower, int Upper)
{

	fLower=Lower;
	fUpper=Upper;

	int NumNodes=fUpper-fLower+1;

	int sd=fModel->GetSD();
	int nnd=fModel->GetNND();

	int nrow=NumNodes*sd;
	int ncol=nnd*sd;

	if(fG_1!=NULL)
	{
		delete [] fG_1;
		fG_1=NULL;
	}

	if(fG_2!=NULL)
	{
		delete [] fG_2;
		fG_2=NULL;
	}

	if(fH_1!=NULL)
	{
		delete [] fH_1;
		fH_1=NULL;
	}

	if(fH_2!=NULL)
	{
		delete [] fH_2;
		fH_2=NULL;
	}

	fG_1=new double[nrow*ncol];
	fG_2=new double[nrow*ncol];


	fH_1=new double[nrow*ncol];
	fH_2=new double[nrow*ncol];
}



void ElastoDynamicsT::G_H_Driver(double t, double t_1)
{

	int Interpolation=fModel->GetTimeInterpolation();

	int sd=fModel->GetSD();
	int nnd=fModel->GetNND();
	int elnnd=fModel->GetENND();
	//int etype=fModel->GetElementType();
	int nel=fModel->GetNEL();
	//int nq=fModel->GetNQ();

	int num_dof=sd*nnd;
	int loc_dof=sd*(fUpper-fLower+1);

	/*initialize G and H*/
	MathOperationT::VecSet(loc_dof*num_dof, fG_1, 0);
	MathOperationT::VecSet(loc_dof*num_dof, fG_2, 0);
	MathOperationT::VecSet(loc_dof*num_dof, fH_1, 0);
	MathOperationT::VecSet(loc_dof*num_dof, fH_2, 0);


	const int* IEN=fModel->GetIEN();
	const double* Coords=fModel->GetCoords();

	//MathOperationT::PrintMatrix(nnd, sd, Coords, "coordinates");

	double* ele_nodes= new double[sd*elnnd]; //Element coordinates
	double* collocation= new double[sd]; //coordinates of the collocation point

	//pass parameters to element
	fElement->SetCoords(ele_nodes);
	fElement->SetCollocation(collocation);
	fElement->SetTime(t,t_1);

	for(int n_i=fLower; n_i<=fUpper; n_i++)
	{
		for(int d_i=0; d_i<sd; d_i++)
			collocation[d_i]=Coords[n_i*sd+d_i];

		//cout << "n_i is " << n_i << " fLower is " << fLower << " fUpper is " << fUpper << endl;
		//MathOperationT::PrintVector(sd, collocation, "collocation");

		for(int e_i=0; e_i<nel; e_i++)
		{
			//******
			//extract coordinates for element nodes
			for(int a=0; a<elnnd; a++)
			{
				int gid=IEN[e_i*elnnd+a];

				for(int di=0;di<sd;di++)
					ele_nodes[a*sd+di]=Coords[gid*sd+di];
			}
            //MathOperationT::PrintMatrix(elnnd, sd, ele_nodes, "ele_nodes");
            
            //*******
            // determine if the element is singular
            int singular=-1;
            
            for (int a=0; a<elnnd; a++)
            {
                double r2=pow(ele_nodes[a*sd]-collocation[0],2)+pow(ele_nodes[a*sd+1]-collocation[1],2)+pow(ele_nodes[a*sd+2]-collocation[2],2);
                
                if (r2<1e-9)
                    singular=a;
                
            }
            //cout << "singular is " << singular << endl;


				if(singular==-1)
					fElement->FormGH_Regular();
				else if(singular>=0 && singular<elnnd)
					fElement->FormGH_Singular(singular);
				else
					throw "ElastoDynamicsT::G_H_Driver, singular is not correctly determined";

				const double* GE_1;
				const double* GE_2;
				const double* HE_1;
				const double* HE_2;
				const double* C_1;
				const double* C_2;
				fElement->GetGHLinear(&GE_1, &GE_2, &HE_1, &HE_2, &C_1, &C_2);

				//cout <<"element " << e_i+1 << endl;
				//MathOperationT::PrintMatrix(3, elnnd*3, GE_1, "local G_1");
				//MathOperationT::PrintMatrix(3, elnnd*3, GE_2, "local G_2");
				//MathOperationT::PrintMatrix(3, elnnd*3, HE_1, "local H_1");
				//MathOperationT::PrintMatrix(3, elnnd*3, HE_2, "local H_2");
				//MathOperationT::PrintMatrix(3, C_1, "local C_1");
		                AssembleGH(fG_1, GE_1, n_i, e_i);
                    		AssembleGH(fG_2, GE_2, n_i, e_i);
				AssembleGH(fH_1, HE_1, n_i, e_i);
				AssembleGH(fH_2, HE_2, n_i, e_i);
				AssembleC(fH_1, C_1, n_i);
				AssembleC(fH_2, C_2, n_i);


		} //end loop over elements

		int Localid;

		if(fModel->GetIfExt())
		{
				if (t==t_1)
				{
					Localid=n_i-fLower;
					fH_1[(Localid*sd+0)*(nnd*sd) + (n_i*sd)+0]+=0.5*t_1;
					fH_1[(Localid*sd+1)*(nnd*sd) + (n_i*sd)+1]+=0.5*t_1;
					fH_1[(Localid*sd+2)*(nnd*sd) + (n_i*sd)+2]+=0.5*t_1;

					fH_2[(Localid*sd+0)*(nnd*sd) + (n_i*sd)+0]+=0.5*t_1;
					fH_2[(Localid*sd+1)*(nnd*sd) + (n_i*sd)+1]+=0.5*t_1;
					fH_2[(Localid*sd+2)*(nnd*sd) + (n_i*sd)+2]+=0.5*t_1;
				}
		}
	} //end loop over nodes


	delete[] ele_nodes;
	delete[] collocation;
}



void ElastoDynamicsT::AssembleGH(double* GlobalMatrix, const double* ElementMatrix, int cid, int eid)
{
	int LocalCid=cid-fLower;

	int sd=fModel->GetSD();
	int nnd=fModel->GetNND();
	int elnnd=fModel->GetENND();
	const int* IEN=fModel->GetIEN();

	for(int ni=0; ni<elnnd; ni++)
	{
		int nid=IEN[eid*elnnd+ni]; //node id

		for(int i=0; i<sd; i++)
		{
			for(int j=0; j<sd; j++)
			{
				GlobalMatrix[(i+LocalCid*sd)*(nnd*sd) + (nid*sd)+j]+= ElementMatrix[i*(elnnd*sd)+ (ni*sd)+j];
			}

		}

	}

}

void ElastoDynamicsT::AssembleC(double* GlobalMatrix, const double* CMatrix, int cid)
{
	int LocalCid=cid-fLower;

	int sd=fModel->GetSD();
	int nnd=fModel->GetNND();
	//int elnnd=fModel->GetENND();

	for(int i=0; i<sd; i++)
	{
		for(int j=0; j<sd; j++)
		{
			int P=LocalCid*sd+i;
			int Q=cid*sd+j;

			GlobalMatrix[P*(nnd*sd)+Q]+=CMatrix[i*sd+j];
		}
	}

}


void ElastoDynamicsT::GetGHLinear(double** G_1, double** G_2, double** H_1, double** H_2)
{
	*G_1=fG_1;
	*G_2=fG_2;
	*H_1=fH_1;
	*H_2=fH_2;
}











