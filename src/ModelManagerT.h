/*
 * ModelManagerT.h
 *
 *  Created on: May 28, 2015
 *      Author: xiaoyong
 */

#ifndef MODELMANAGERT_H_
#define MODELMANAGERT_H_

#include "InputT.h"

namespace TD_BEM{


class ModelManagerT
{
public:

	ModelManagerT();
	~ModelManagerT();

	int GetFormulation(void);
	void SetFormulation(int formulation);

	int GetNumStep(void);
	void SetNumStep(int numstep);

	double GetDT(void);
	void SetDT(double dt);

	const int GetTimeInterpolation(void);
	void SetTimeInterpolation(int interpolation);

	const int GetNQ(void);
	void SetNQ(int nq);

	/**
	 * @Functions to get model parameters
	 * 		GetSD=return spatial dimensions
	 * 		GetNND=return number of nodes
	 * 		GetElementType=return type of the element
	 * 		GetNEL=return number of elements
	 * 		GetIEN=return element connectivities
	 * 		GetCoords=return nodes' coordinates*/
	int GetSD(void);
	void SetSD(int sd);

	int GetIfExt(void);
	void SetIfExt(int if_ext);

	int GetNND(void);
	void SetNND(int nnd);

	int GetElementType(void);
	void SetElementType(int et);

	int GetENND(void);
	void SetENND(int ENND);

	int GetNEL(void);
	void SetNEL(int nel);

	const int* GetIEN(void);
	void SetIEN(int* IEN);

	const double* GetCoords(void);
	void SetCoords(double* coords);



	int GetNND_Aux(void);
	void SetNND_Aux(int nnd);

	int GetElementType_Aux(void);
	void SetElementType_Aux(int et);

	int GetENND_Aux(void);
	void SetENND_Aux(int ENND);

	int GetNEL_Aux(void);
	void SetNEL_Aux(int nel);

	const int* GetIEN_Aux(void);
	void SetIEN_Aux(int* IEN);

	const double* GetCoords_Aux(void);
	void SetCoords_Aux(double* coords);



	void GetMaterialConstants(double& E, double& nu, double& rho);
	void SetMaterialConstants(double E, double nu, double rho);

	void BoundaryCondition_Global(int& UBC_num, int** UBC_DOFs, double** UBC_Values,
						   int& FBC_num, int** FBC_DOFs, double** FBC_Values);
	void SetBCs(int UBC_num, int* UBC_DOFs, double* UBC_values,
				int FBC_num, int* FBC_DOFs, double* FBC_values);


	int GetNumInterface(void);
	void SetNumInterface(int num_interface);

	const int* GetInterface(void);
	void SetInterface(int* Interface);

	void ReadInput(const char* name);

	/*distribute boundary conditions into processors*/
	void Distributed_BC(int low, int high);

	void BoundaryCondition_Local(int& UBC_num, int** UBC_DOFs, double** UBC_Values,
						   int& FBC_num, int** FBC_DOFs, double** FBC_Values);

private:
	/**
	 * algorithm parameters
	 */
	int fFormulation;
	double fDT;
	int fNumStep;
	int fTimeInterpolation;
	int fNQ;

	/*
	 * fSD= spatial dimension (currently only three dimension case is supported)
	 * fNND= number of nodes in the model
	 * fENND= number of node per element
	 * fNEL= number of elements in the model
	 * fET= element type
	 * fIEN= element connectivities
	 * fCoords= nodes coordinates
	 * fTimeInterpolation= type of time interpolation, for example, constant, linear, b-spline
	 * fNQ= number of Gaussian quadrature point per element */
	int fSD;
	int fIfExt;
	int fNND;
	int fENND;
	int fNEL;
	int fET;
	int* fIEN;
	double* fCoords;



	int fNND_Aux;
	int fENND_Aux;
	int fNEL_Aux;
	int fET_Aux;
	int* fIEN_Aux;
	double* fCoords_Aux;




	/*
	 * fE=Young's modulus
	 * fNu=Possion's ratio
	 * fDensity=density
	 */
	double fE;
	double fNu;
	double fDensity;

	/*
	 * fUBC_num=number of displacement boundary condition
	 * fUBC_DOFs=the specified displacement dofs
	 * fUBC_Values=the specified displacement values
	 * fFBC_num=number of traction boundary condition
	 * fFBC_DOFs=the specified traction dofs
	 * fFBC_Values=the specified traction values
	 */
	int fUBC_num_global;
	int* fUBC_DOFs_global;
	double* fUBC_Values_global;
	int fFBC_num_global;
	int* fFBC_DOFs_global;
	double* fFBC_Values_global;


	int fUBC_num_local;
	int* fUBC_DOFs_local;
	double* fUBC_Values_local;
	int fFBC_num_local;
	int* fFBC_DOFs_local;
	double* fFBC_Values_local;

	int fNumInterface;
	int* fInterface;

	InputT fInput;
};


} /* name space */



#endif /* MODELMANAGERT_H_ */
