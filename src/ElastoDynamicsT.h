/*
 * ElastoDynamicsT.h
 *
 *  Created on: May 28, 2015
 *      Author: xiaoyong
 */

#ifndef ELASTODYNAMICST_H_
#define ELASTODYNAMICST_H_

#include "ModelManagerT.h"
#include "TDElementT.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscsys.h"

namespace TD_BEM{

class ElastoDynamicsT
{
public:
	ElastoDynamicsT(ModelManagerT* Model);
	~ElastoDynamicsT();

	/**
	 * @Function=set fLower and fUpper */
	void SetLowerUpper(int Lower, int upper);

	/**
	 * @Function=drive to formulate G and H matrices, for linear interpolation
	 */
	void G_H_Driver(double t, double t_1);


	/**
	 * @Function=form the B matrix, when solve half space problem with full space Green's function
	 */
	void FormGeomBMatrix(void);


	/**
	 * @Function=Assemble element G and H into global matrix
	 *
	 * @Input:
	 * 		cid=id of the collocation node
	 * 		eid=id of element */
	void AssembleGH(double* GlobalMatrix, const double* ElementMatrix, int cid, int eid);

	void AssembleC(double* GlobalMatirx, const double* CMatrix, int cid);

	void GetGHLinear(double** G_1, double** G_2, double** H_1, double** H_2);

	//void GetGHLinear(double** G_1, double** G_2, double** H_1, double** H_2);


private:
	/**
	 * fModel=the model manager that contains information of mesh
	 * fLower, fUpper= Id range of nodes that are dealt with by the local processor
	 * G_1, ..= temporally store the G matrices formed in current time step
	 * H_1, ..= temporally store the H matrices formed in current time step
	 */
	ModelManagerT* fModel;
	TDElementT* fElement;



	int fLower, fUpper;

	double* fG_1;
	double* fG_2;

	double* fH_1;
	double* fH_2;

	double* fGeometry_B;

};

} /* name space */

#endif /* ELASTODYNAMICST_H_ */
