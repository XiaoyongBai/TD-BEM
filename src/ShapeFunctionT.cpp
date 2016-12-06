/*
 * ShapeFunctionT.cpp
 *
 *  Created on: May 29, 2015
 *      Author: xiaoyong
 */

#include "ShapeFunctionT.h"
#include "MathOperationT.h"
#include <iostream>
#include <cmath>

using namespace TD_BEM;
using namespace GreenFunction;
using namespace std;

ShapeFunctionT::ShapeFunctionT(int type)
{
	fType=type;

	fW=1;
	fJ=1;
	fCoords=NULL;

	switch (fType) {
	case 1:
		fENND=4;
		break;
	case 2:
		fENND=8;
		break;
	default:
		throw "ShapeFunctionT::ShapeFunctionT, fType is unsupported" ;
	}

	int sd=3;

	fX=new double[sd];
	fOutNormal=new double[sd];
	fN=new double[fENND];
	fDN=new double[fENND*2];

	/* initialize the Gaussian points*/
	fNQ=6;
	fShapeCoords=new double[fNQ];
	fWeights=new double[fNQ];

    double shape_coords[6]={
       		0.6612093864662645, -0.6612093864662645, -0.2386191860831969,
              0.2386191860831969, -0.9324695142031521,  0.9324695142031521};

    double weights[6]={
  		  	0.3607615730481386,0.3607615730481386,0.4679139345726910,
              0.4679139345726910,0.1713244923791704,0.1713244923791704};

    MathOperationT::MemCopy(fNQ, shape_coords, fShapeCoords);
    MathOperationT::MemCopy(fNQ, weights, fWeights);

}

ShapeFunctionT::~ShapeFunctionT()
{
	delete[] fX;
	delete[] fOutNormal;
	delete[] fN;
	delete[] fDN;
	delete[] fShapeCoords;
	delete[] fWeights;
}

void ShapeFunctionT::SetNQ(int nq)
{
	fNQ=nq;

	if (fShapeCoords!=NULL) delete[] fShapeCoords;
	if (fWeights!=NULL) delete[] fWeights;

	fShapeCoords=new double[fNQ];
	fWeights=new double[fNQ];

	switch (fNQ){

	case 6:
	{
	      double shape_coords[6]={
	         		0.6612093864662645, -0.6612093864662645, -0.2386191860831969,
	                0.2386191860831969, -0.9324695142031521,  0.9324695142031521};

	      double weights[6]={
	    		  	0.3607615730481386,0.3607615730481386,0.4679139345726910,
	                0.4679139345726910,0.1713244923791704,0.1713244923791704};

	      MathOperationT::MemCopy(fNQ, shape_coords, fShapeCoords);
	      MathOperationT::MemCopy(fNQ, weights, fWeights);
	}
	      break;


    case 10:
    {
           double shape_coords[10]={
        		    -0.148874339, 0.148874339,-0.433395394, 0.433395394,-0.679409568,
                    0.679409568,-0.865063367, 0.865063367,-0.973906529, 0.973906529};
           double weights[10]={
        		   0.295524225,0.295524225,0.269266719,0.269266719,0.219086363,
        		   0.219086363,0.149451349,0.149451349,0.066671344,0.066671344};

           MathOperationT::MemCopy(fNQ, shape_coords, fShapeCoords);
           MathOperationT::MemCopy(fNQ, weights, fWeights);
    }

           break;

    case 20:
    {
    		double shape_coords[20]={
    				-0.076526521, 0.076526521,-0.227785851, 0.227785851,
    				-0.373706089, 0.373706089,-0.510867002, 0.510867002,
    				-0.636053681, 0.636053681,-0.746331906, 0.746331906,
    				-0.839116972, 0.839116972,-0.912234428, 0.912234428,
    				-0.963971927, 0.963971927, 0.993128599, 0.993128599};
    		double weights[20]={
    				0.152753387,  0.152753387, 0.149172986, 0.149172986,
    				0.142096109,  0.142096109, 0.131688638, 0.131688638,
    				0.118194532,  0.118194532, 0.10193012,  0.10193012,
    				0.083276742,  0.083276742, 0.062672048, 0.062672048,
    				0.04060143,   0.04060143,  0.017614007, 0.017614007};

            MathOperationT::MemCopy(fNQ, shape_coords, fShapeCoords);
            MathOperationT::MemCopy(fNQ, weights, fWeights);
    }
    		break;

    case 40:
    {
    	double shape_coords[40]={
    			-0.0387724175060508, 0.0387724175060508, -0.1160840706752552, 0.1160840706752552,
    			-0.1926975807013711, 0.1926975807013711, -0.2681521850072537, 0.2681521850072537,
    			-0.3419940908257585, 0.3419940908257585, -0.4137792043716050, 0.4137792043716050,
    			-0.4830758016861787, 0.4830758016861787, -0.5494671250951282, 0.5494671250951282,
    			-0.6125538896679802, 0.6125538896679802, -0.6719566846141796, 0.6719566846141796,
    			-0.7273182551899271, 0.7273182551899271, -0.7783056514265194, 0.7783056514265194,
    			-0.8246122308333117, 0.8246122308333117, -0.8659595032122595, 0.8659595032122595,
    			-0.9020988069688743, 0.9020988069688743, -0.9328128082786765, 0.9328128082786765,
    			-0.9579168192137917, 0.9579168192137917, -0.9772599499837743, 0.9772599499837743,
    			-0.9907262386994570, 0.9907262386994570, -0.9982377097105593, 0.9982377097105593};

    	double weights[40]={
    			0.0775059479784248, 0.0775059479784248, 0.0770398181642480, 0.0770398181642480,
    			0.0761103619006262, 0.0761103619006262, 0.0747231690579683, 0.0747231690579683,
    			0.0728865823958041, 0.0728865823958041, 0.0706116473912868, 0.0706116473912868,
    			0.0679120458152339, 0.0679120458152339, 0.0648040134566010, 0.0648040134566010,
    			0.0613062424929289, 0.0613062424929289, 0.0574397690993916, 0.0574397690993916,
    			0.0532278469839368, 0.0532278469839368,	0.0486958076350722, 0.0486958076350722,
    			0.0438709081856733, 0.0438709081856733,	0.0387821679744720, 0.0387821679744720,
    			0.0334601952825478, 0.0334601952825478,	0.0279370069800234, 0.0279370069800234,
    			0.0222458491941670, 0.0222458491941670, 0.0164210583819079, 0.0164210583819079,
    			0.0104982845311528, 0.0104982845311528,	0.0045212770985332, 0.0045212770985332};

        MathOperationT::MemCopy(fNQ, shape_coords, fShapeCoords);
        MathOperationT::MemCopy(fNQ, weights, fWeights);

    }
    		break;
    default:
    	throw "ShapeFunctionT::SetNQ, Unsupported number of Gaussian points";
	}

}


int ShapeFunctionT::GetNQ()
{
	return fNQ;
}

int ShapeFunctionT::GetENND()
{
	return fENND;
}

const double* ShapeFunctionT::GetShapeValues()
{
	return fN;
}

double ShapeFunctionT::GetJacobian()
{
	return fJ;
}

double ShapeFunctionT::GetWeight()
{
	return fW;
}

const double* ShapeFunctionT::GetX()
{
	return fX;
}

const double* ShapeFunctionT::GetOutNormal()
{
	return fOutNormal;
}

void ShapeFunctionT::GetShapeCoords(int q_1, int q_2, double& eta_1, double& eta_2, double& w)
{
	eta_1=fShapeCoords[q_1];
	eta_2=fShapeCoords[q_2];
	w=fWeights[q_1]*fWeights[q_2];
}

void ShapeFunctionT::GetShapeCoords(int q_1, double& eta_1, double& w)
{
    eta_1=fShapeCoords[q_1];
    w=fWeights[q_1];
}


void ShapeFunctionT::SetCoords(double* coords)
{
	fCoords=coords;

	//MathOperationT::PrintMatrix(4, 3, fCoords, "shape fCoords, for sub element" );

}

void ShapeFunctionT::Evaluate(int q_1, int q_2)
{
	Evaluate(fShapeCoords[q_1], fShapeCoords[q_2]);
	//MathOperationT::PrintVector(fNQ, fWeights, "fWeights");
	fW=fWeights[q_1]*fWeights[q_2];
}

void ShapeFunctionT::Evaluate(double eta_1, double eta_2)
{
	int sd=3;

	switch(fType){
	case 1:
	{
		//%%%%compute shape functions
		fN[0]=0.25*(1-eta_1)*(1-eta_2);
		fN[1]=0.25*(1+eta_1)*(1-eta_2);
		fN[2]=0.25*(1+eta_1)*(1+eta_2);
		fN[3]=0.25*(1-eta_1)*(1+eta_2);

		MathOperationT::VecSet(sd, fX, 0);

//MathOperationT::PrintMatrix(4, 3, fCoords, "shape fCoords, for sub element" );


		for(int a=0; a<fENND; a++)
		{
			for(int d=0; d<sd; d++)
			{
				fX[d]+=fN[a]*fCoords[a*sd+d];
			}
		}


		fDN[0]=-0.25*(1-eta_2);
		fDN[1]= 0.25*(1-eta_2);
		fDN[2]= 0.25*(1+eta_2);
		fDN[3]=-0.25*(1+eta_2);

		fDN[4]=-0.25*(1-eta_1);
		fDN[5]=-0.25*(1+eta_1);
		fDN[6]= 0.25*(1+eta_1);
		fDN[7]= 0.25*(1-eta_1);

		double D_X_1[3]={0,0,0};
		double D_X_2[3]={0,0,0};

		for(int a=0; a<fENND; a++)
		{
			for(int d=0; d<sd; d++)
			{
				D_X_1[d]+=fDN[a]*fCoords[a*sd+d];
				D_X_2[d]+=fDN[a+fENND]*fCoords[a*sd+d];
			}
		}

		MathOperationT::VecCross(D_X_1, D_X_2, fOutNormal);
		fJ=MathOperationT::VecNormal(sd, fOutNormal);

		for(int d=0; d<sd; d++)
			fOutNormal[d] /=fJ;


		//MathOperationT::PrintVector(fENND, fN, "fN");
		//MathOperationT::PrintVector(sd, fX, "fX" );
		//MathOperationT::PrintVector(fENND*2, fDN, "fDN");
		//MathOperationT::PrintVector(3, D_X_1, "dx1");
		//MathOperationT::PrintVector(3, D_X_2, "dx2");
		//MathOperationT::PrintVector(3, fOutNormal, "fOutNormal");
		//cout << "fJ is " << fJ <<endl;

		break;
	}
	case 2:
	{
		fN[0]=0.25*(1-eta_1)*(1-eta_2)*(-eta_1-eta_2-1);
		fN[1]=0.25*(1+eta_1)*(1-eta_2)*( eta_1-eta_2-1);
		fN[2]=0.25*(1+eta_1)*(1+eta_2)*( eta_1+eta_2-1);
		fN[3]=0.25*(1-eta_1)*(1+eta_2)*(-eta_1+eta_2-1);
		fN[4]=0.5*(1-pow(eta_1,2))*(1-eta_2);
		fN[5]=0.5*(1+eta_1)*(1-pow(eta_2,2));
		fN[6]=0.5*(1-pow(eta_1,2))*(1+eta_2);
		fN[7]=0.5*(1-eta_1)*(1-pow(eta_2,2));

		MathOperationT::VecSet(sd, fX, 0);
		for(int a=0; a<fENND; a++)
		{
			for(int d=0; d<sd; d++)
			{
				fX[d]+=fN[a]*fCoords[a*sd+d];
			}
		}

        fDN[0]= -0.25*(1 - eta_1)*(1 - eta_2) - 0.25*(1 - eta_2)*(-1 - eta_1 - eta_2);
        fDN[1]= 0.25*(1 + eta_1)*(1 - eta_2) + 0.25*(1 - eta_2)*(-1 + eta_1 - eta_2);
        fDN[2]=0.25*(1 + eta_1)*(1 + eta_2) + 0.25*(1 + eta_2)*(-1 + eta_1 + eta_2);
        fDN[3]=-0.25*(1 - eta_1)*(1 + eta_2) - 0.25*(1 + eta_2)*(-1 - eta_1 + eta_2);
        fDN[4]=-eta_1*(1 - eta_2);
        fDN[5]=0.5*(1 - pow(eta_2,2));
        fDN[6]=-eta_1*(1 + eta_2);
        fDN[7]=-0.5*(1 - pow(eta_2,2));

        fDN[8]=-0.25*(1 - eta_1)*(1 - eta_2) - 0.25*(1 - eta_1)*(-1 - eta_1 - eta_2);
        fDN[9]=-0.25*(1 + eta_1)*(1 - eta_2) - 0.25*(1 + eta_1)*(-1 + eta_1 - eta_2);
        fDN[10]=0.25*(1 + eta_1)*(1 + eta_2) + 0.25*(1 + eta_1)*(-1 + eta_1 + eta_2);
        fDN[11]=0.25*(1 - eta_1)*(1 + eta_2) + 0.25*(1 - eta_1)*(-1 - eta_1 + eta_2);
        fDN[12]=-0.5*(1 - pow(eta_1,2));
        fDN[13]=-(1 + eta_1)*eta_2;
        fDN[14]=0.5*(1 - pow(eta_1,2));
        fDN[15]=-(1 - eta_1)*eta_2;

		double D_X_1[3]={0,0,0};
		double D_X_2[3]={0,0,0};

		for(int a=0; a<fENND; a++)
		{
			for(int d=0; d<sd; d++)
			{
				D_X_1[d]+=fDN[a]*fCoords[a*sd+d];
				D_X_2[d]+=fDN[a+fENND]*fCoords[a*sd+d];
			}
		}

		MathOperationT::VecCross(D_X_1, D_X_2, fOutNormal);
		fJ=MathOperationT::VecNormal(sd, fOutNormal);

		for(int d=0; d<sd; d++)
			fOutNormal[d] /=fJ;

		break;
	}


	default:
		throw "ShapeFunctionT::Evaluate(double, double), unsupported Element Type";

	}
}



