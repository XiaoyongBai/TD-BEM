/*
 * ShapeFunctionT.h
 *
 *  Created on: May 29, 2015
 *      Author: xiaoyong
 */

#ifndef SHAPEFUNCTIONT_H_
#define SHAPEFUNCTIONT_H_

namespace TD_BEM{

class ShapeFunctionT
{
public:
	ShapeFunctionT(int type);

	/**
	 * @Function=set number of Gaussian points */
	void SetNQ(int nq);
	int GetNQ(void);

	int GetENND(void);

	const double* GetShapeValues(void);
	double GetJacobian(void);
	double GetWeight(void);
	const double* GetX(void);
	const double* GetOutNormal(void);
	/**
	 * @Function=get shape coordinates of a given Gaussian point */
	void GetShapeCoords(int q_1, int q_2, double& eta_1, double& eta_2, double& w);
    
    void GetShapeCoords(int q_1, double& eta_1, double& w);

	void SetCoords(double* coords);

	/**
	 * @Function=evaluate shape functions, Jacobian, etc, at the given Gaussian points
	 *
	 * @Input:
	 * 		q_1, q_2=the number of Gaussian points*/
	void Evaluate(int q_1, int q_2);

	/**
	 * @Function=evaluate shape functions, Jacobian, etc, at a point with specified shape coords
	 *
	 * @Input:
	 * 		eta_1=shape coordinate in the first direction
	 * 		eta_2=shape coordinate in the second direction */
	void Evaluate(double eta_1, double eta_2);

	~ShapeFunctionT();

private:
	/**
	 * fType=type of the element
	 * fENND=number of nodes in this element
	 * fNQ=number of Gaussian points
	 * fCoords=global nodal coordinates
	 * fShapeCoords=Shape coordinates of nodals
	 * fWeights=vector of weights
	 */
	int fType;
	int fENND;
	int fNQ;

	double* fCoords;
	double* fShapeCoords;
	double* fWeights;

	/**
	 * fX=coordinates of current Gaussian point, obtained by interpolation
	 * fOutNormal=out normal of current Gaussian point
	 * fN=values of shape functions at current Gaussian point
	 * fDN=derivatives of shape functions with respect to shape coordinates
	 * fW=weight of current Gaussian point
	 * fJ=Jacobian of current Gaussian point
	 */
	double* fX;
	double* fOutNormal;
	double* fN;
	double* fDN;
	double 	fW;
	double fJ;

};



}/* name space */




#endif /* SHAPEFUNCTIONT_H_ */
