/*
 * TDElementT.h
 *
 *  Created on: May 29, 2015
 *      Author: xiaoyong
 */

#ifndef TDELEMENTT_H_
#define TDELEMENTT_H_

#include "ShapeFunctionT.h"
#include "GreenFunctionT.h"

namespace TD_BEM{

class TDElementT
{
public:
	TDElementT(int type);
	virtual ~TDElementT();

	/**
	 * @Function=set parameters */
	void SetCollocation(double* collocation);
	void SetCoords(double* coords);
	virtual void SetTime(double t, double t_1);
	void SetNQ(int nq);

	void SetMaterial(double E, double nu, double density);

	/**
	 * @Function=create a shape function and allocate matrices*/
	void SetShape(int type);

	/**
	 * @Function=form G and H matrices for an element, when linear time interpolation is used */
	virtual void FormGH_Regular()=0;
	virtual void FormGH_Singular(int singular)=0;

	/**
	 * @Function=form B for an element
	 */
	virtual void FormGeomBMatrix();

	void GetGeomBMatrix(const double** GeomB);

	void GetGHLinear(const double** GE_1, const double** GE_2, const double** HE_1, const double** HE_2,
					 const double** C_1, const double** C_2);

	// get Green's function
	GreenFunctionT* GetGreenFunction(void);

protected:
/*
	 * Material constants
	 */
	double fE;
	double fNu;
	double fDensity;

	/**
	 * fType=type of the element
	 * fNQ=number of Gaussian points
	 * fCollocation=coordinates of collocation points
	 * fCoords=coordinates of nodes
	 *
	 * fT, fT_1=time
	 *
	 * fGE_1,..=temporally store element G matrices
	 * fHE_1,..=temporally store element H matrices
	 * fHC_1,..=temporally store element C matrices */

	int fType;
	int fNQ;

	double* fCollocation;
	double* fCoords;

	double fT;
	double fT_1;

	double* fGE_1;
	double* fGE_2;


	double* fHE_1;
	double* fHE_2;


	double* fC_1;
	double* fC_2;

	double* fGeomB;

	ShapeFunctionT* fShape;
	GreenFunctionT* fGreenFunction;
    
	/*
	 * fSubNum=number of sub elements
	 * fSubType=type of sub elements
	 * fSubCID=local collocation elements for each sub element
	 * fSubNodes=node coordinates of each sub element
	 * fSubIDs=node id (in quadratic elements) of each sub element
	 */
	int fSubNum;
	int fSubType;
	int* fSubCID;

	double** fSubNodes;

	int** fSubIDs;

	ShapeFunctionT* fSubShape;

};

} /* name space */

#endif /* TDELEMENTT_H_ */
