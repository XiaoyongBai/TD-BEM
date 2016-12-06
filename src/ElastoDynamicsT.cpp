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

	int interpolation=fModel->GetTimeInterpolation();

	switch(interpolation)
	{
	case 1:
		//fElement=new TDElementHeavisideConstantT(etype);
		break;
	case 2:
		fElement=new TDElementHeavisideLinearT(etype);
		break;
	case 3:
		//fElement=new TDElementStepConstantT(etype);
		break;
	case 4:
		fElement=new TDElementStepLinearT(etype);
		break;
	case 5:
		//fElement=new TDElementRampHeaviLinearT(etype);
		break;
	default:
		throw "ElastoDynamicT::ElastoDynamicsT, Unsupported element type";
	}

	fElement->SetNQ(fModel->GetNQ());

	double E, nu, rho;
	fModel->GetMaterialConstants(E, nu, rho);
	fElement->SetMaterial(E, nu, rho);

	fG_1=NULL;
	fG_2=NULL;
	fH_1=NULL;
	fH_2=NULL;

	fGeometry_B=NULL;

	SetLowerUpper(0, fModel->GetNND()-1);

	cout << "ElastoDynamicsT is constructed\n";
}

ElastoDynamicsT::~ElastoDynamicsT()
{

	delete[] fG_1;
	delete[] fG_2;

	delete[] fH_1;
	delete[] fH_2;

	if (fGeometry_B) delete[] fGeometry_B;

	if(fElement)
	{
		delete fElement;
		fElement=NULL;
	}
	//delete fModel;
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


	if (fModel->GetFormulation()==3)
	{
		if(fGeometry_B!=NULL)
		{
			delete[] fGeometry_B;
			fGeometry_B=NULL;
		}

		fGeometry_B=new double[nrow*ncol];
	}

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
            

			switch (Interpolation){
			case 2:
				{
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

					//MathOperationT::PrintMatrix(loc_dof, num_dof, fG_1, "Global G_1");
					//MathOperationT::PrintMatrix(loc_dof, num_dof, fG_2, "Global G_2");
					//MathOperationT::PrintMatrix(loc_dof, num_dof, fH_1, "Global H_1");
					//MathOperationT::PrintMatrix(loc_dof, num_dof, fH_2, "Global H_2");

					//MathOperationT::PrintMatrix(loc_dof, num_dof, fH_1, "Global H_1");
					//MathOperationT::PrintMatrix(loc_dof, num_dof, fH_2, "Global H_2");
				}
				break;

			case 4:
				{
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

			}
				break;

			case 5:
			{
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
				//MathOperationT::PrintMatrix(3, C_2, "local C_2");

				AssembleGH(fG_1, GE_1, n_i, e_i);
				AssembleGH(fG_2, GE_2, n_i, e_i);
				AssembleGH(fH_1, HE_1, n_i, e_i);
				AssembleGH(fH_2, HE_2, n_i, e_i);

				if(singular==-1)
				{
					AssembleC(fH_1, C_1, n_i);
					AssembleC(fH_2, C_2, n_i);
				}

				//MathOperationT::PrintMatrix(loc_dof, num_dof, fG_1, "Global G_1");
				//MathOperationT::PrintMatrix(loc_dof, num_dof, fG_2, "Global G_2");
				//MathOperationT::PrintMatrix(loc_dof, num_dof, fH_1, "Global H_1");
				//MathOperationT::PrintMatrix(loc_dof, num_dof, fH_2, "Global H_2");

				//MathOperationT::PrintMatrix(loc_dof, num_dof, fH_1, "Global H_1");
				//MathOperationT::PrintMatrix(loc_dof, num_dof, fH_2, "Global H_2");

			}
				break;

			default:
				throw "ElastoDynamicsT::G_H_Driver, unsupported time interpolation";
			} //end switch interpolation


		} //end loop over elements

		int Localid;

		if(fModel->GetIfExt())
		{
			switch(Interpolation)
			{
			case 2:
				Localid=n_i-fLower;
				fH_1[(Localid*sd+0)*(nnd*sd) + (n_i*sd)+0]+=0.5*t_1;
				fH_1[(Localid*sd+1)*(nnd*sd) + (n_i*sd)+1]+=0.5*t_1;
				fH_1[(Localid*sd+2)*(nnd*sd) + (n_i*sd)+2]+=0.5*t_1;

				fH_2[(Localid*sd+0)*(nnd*sd) + (n_i*sd)+0]+=0.5*t_1;
				fH_2[(Localid*sd+1)*(nnd*sd) + (n_i*sd)+1]+=0.5*t_1;
				fH_2[(Localid*sd+2)*(nnd*sd) + (n_i*sd)+2]+=0.5*t_1;
				break;

			case 4:
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
				break;
			case 5:
				if (t==t_1)
				{
					double g1, g2;
					GreenFunctionT* GF;
					GF= fElement->GetGreenFunction();
					GF->Int_CubicB_Linear(g1, g2);

					Localid=n_i-fLower;
					fH_1[(Localid*sd+0)*(nnd*sd) + (n_i*sd)+0]+=g1;
					fH_1[(Localid*sd+1)*(nnd*sd) + (n_i*sd)+1]+=g1;
					fH_1[(Localid*sd+2)*(nnd*sd) + (n_i*sd)+2]+=g1;

					fH_2[(Localid*sd+0)*(nnd*sd) + (n_i*sd)+0]+=g2;
					fH_2[(Localid*sd+1)*(nnd*sd) + (n_i*sd)+1]+=g2;
					fH_2[(Localid*sd+2)*(nnd*sd) + (n_i*sd)+2]+=g2;
				}
				break;
			default:
				throw "ElastodynamicsT::G_H_Driver, unsupported time interpolation";
			}

		}
	} //end loop over nodes


	//MathOperationT::PrintMatrix(9, 693, fGeometry_B, "fGeometry_B");

	//MathOperationT::PrintMatrix(72, 72, fH_1, "fH_1");




	if(fModel->GetFormulation()==3)
	{
		switch(Interpolation)
		{
			case 2:
				MathOperationT::VecPlus(loc_dof*num_dof, 1, fH_1, 0.5*t_1, fGeometry_B, fH_1);

				MathOperationT::VecPlus(loc_dof*num_dof, 1, fH_2, 0.5*t_1, fGeometry_B, fH_2);
			break;

			case 4:
				if (t==t_1)
				{
					double alpha=0.5*t_1;

					MathOperationT::VecPlus(loc_dof*num_dof, 1, fH_1, alpha, fGeometry_B, fH_1);

					MathOperationT::VecPlus(loc_dof*num_dof, 1, fH_2, alpha, fGeometry_B, fH_2);

				}
			break;

			default:
				throw "ElastodynamicsT::G_H_Driver, unsupported time interpolation";
		    }

		}//end if;


	delete[] ele_nodes;
	delete[] collocation;

	/*MathOperationT::PrintMatrix(loc_dof, num_dof, fG_1, "Global G_1");
	MathOperationT::PrintMatrix(loc_dof, num_dof, fG_2, "Global G_2");

	MathOperationT::PrintMatrix(loc_dof, num_dof, fH_1, "Global H_1");
	MathOperationT::PrintMatrix(loc_dof, num_dof, fH_2, "Global H_2");
    */
}






void ElastoDynamicsT::FormGeomBMatrix()
{
	cout <<"Form the Geometry B Matrix\n";

	int sd=fModel->GetSD();
	int nnd=fModel->GetNND();

	const double* Coords=fModel->GetCoords();


	int elnnd_aux=fModel->GetENND_Aux();
	//int etype_aux=fModel->GetElementType_Aux();
	int nel_aux=fModel->GetNEL_Aux();
	//int nnd_aux=fModel->GetNND_Aux();
	//int nq=fModel->GetNQ();

	int num_dof=sd*nnd;
	int loc_dof=sd*(fUpper-fLower+1);

	/*initialize G and H*/
	MathOperationT::VecSet(loc_dof*num_dof, fGeometry_B, 0);


	const int* IEN_Aux=fModel->GetIEN_Aux();
	const double* Coords_Aux=fModel->GetCoords_Aux();

	//MathOperationT::PrintMatrix(nel_aux, elnnd_aux, IEN_Aux, "IEN_Aux");

	//MathOperationT::PrintMatrix(nnd_aux, sd, Coords_Aux, "coordinates_Aux");

	double* ele_nodes_aux= new double[sd*elnnd_aux]; //Element coordinates
	double* collocation= new double[sd]; //coordinates of the collocation point

	//pass parameters to element
	fElement->SetCoords(ele_nodes_aux);
	fElement->SetCollocation(collocation);


	for(int n_i=fLower; n_i<=fUpper; n_i++)
	{
		for(int d_i=0; d_i<sd; d_i++)
			collocation[d_i]=Coords[n_i*sd+d_i];

		//cout << "n_i is " << n_i << " fLower is " << fLower << " fUpper is " << fUpper << endl;
		//MathOperationT::PrintVector(sd, collocation, "collocation");

		for(int e_i=0; e_i<nel_aux; e_i++)
		{

			//******
			//extract coordinates for element nodes
			for(int a=0; a<elnnd_aux; a++)
			{
				int gid=IEN_Aux[e_i*elnnd_aux+a];

				for(int di=0;di<sd;di++)
					ele_nodes_aux[a*sd+di]=Coords_Aux[gid*sd+di];
			}

			//MathOperationT::PrintMatrix(elnnd_aux, sd, ele_nodes_aux, "ele_nodes");

			const double* GeomB;

			fElement->FormGeomBMatrix();
			fElement->GetGeomBMatrix(&GeomB);

			//cout <<"element " << e_i+1 << endl;
			//MathOperationT::PrintMatrix(3, 3, GeomB, "GeomB");

			int Localid=n_i-fLower;
			for (int i=0; i<sd; i++)
			{
				for (int j=0; j<sd; j++)
				{
					int r_num=3*Localid+i;
					int c_num=3*n_i+j;
					fGeometry_B[r_num*num_dof+c_num]+=GeomB[i*sd+j];
				}

			}

		}

	} //end loop over nodes

	//MathOperationT::PrintMatrix(30, 693, fGeometry_B, "Global fGeometry_B");


	delete[] ele_nodes_aux;
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











