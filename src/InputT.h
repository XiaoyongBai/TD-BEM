/*
 * InputT.h
 *
 *  Created on: Jun 5, 2015
 *      Author: xiaoyong
 */

#ifndef INPUTT_H_
#define INPUTT_H_

#include <cstring>
#include <iostream>
#include <fstream>

namespace TD_BEM{

class InputT
{
public:
	InputT();
	~InputT();

	void ReadInput(const char* name);

	bool ReadAlgorithm(char* buff, int line);
	bool ReadGeometry(std::ifstream& fin, char* buff, int& line);
	bool ReadGeometry_Aux(std::ifstream& fin, char* buff, int& line);
	bool ReadMaterial(char* buff, int line);
	bool ReadBCs(std::ifstream& fin, char* buff, int& line);
	bool ReadInterface(std::ifstream& fin, char* buff, int& line);

	/*
	 * Get algorithm parameters
	 */
	int GetFormulation(void);

	int GetNumStep(void);

	double GetDT(void);

	int GetTimeInterpolation(void);

	int GetNQ(void);


	/*
	 * Geometric information
	 */

	int GetSD(void);

	int GetIfExt(void);

	int GetNND(void);

	int GetElementType(void);

	int GetENND(void);

	int GetNEL(void);

	int* GetIEN(void);

	double* GetCoords(void);

	int GetNND_Aux(void);

	int GetElementType_Aux(void);

	int GetENND_Aux(void);

	int GetNEL_Aux(void);

	int* GetIEN_Aux(void);

	double* GetCoords_Aux(void);


	void GetMaterialConstants(double& E, double& nu, double& rho);

	void BoundaryCondition(int& UBC_num, int** UBC_DOFs, double** UBC_Values,
							   int& FBC_num, int** FBC_DOFs, double** FBC_Values);



	void GetInterface(int& num_interface, int** interface);


private:
	/*
	 * algorithm coefficients
	 */
	int fFormulation;
	int fInterpolation;
	int fNumStep;
	double fDT;
	int fNQ;

	/*
	 * material parameters
	 */
	double fE;
	double fNu;
	double fDensity;

	/*
	 * geoemtry information
	 */
	int fSD;
	int fIfExt;
	int fNND;
	int fNEL;
	int fEType;
	int fENND;
	double* fCoords;
	int* fIEN;

	int fNND_Aux;
	int fNEL_Aux;
	int fEType_Aux;
	int fENND_Aux;
	double* fCoords_Aux;
	int* fIEN_Aux;

	int fNumInterface;
	int* fInterface;

	/*
	 * boundary conditions
	 */
	int fUBC_num;
    	int* fUBC_DOFs;
    	double* fUBC_Values;

	int fFBC_num;
    	int* fFBC_DOFs;
    	double* fFBC_Values;
};


}/*namespace*/

#endif /* INPUTT_H_ */










