/*
 * InputT.cpp
 *
 *  Created on: Jun 5, 2015
 *      Author: xiaoyong
 */

#include "InputT.h"

#include <fstream>
#include <cstring>
#include "MathOperationT.h"

#include "petscsys.h"

using namespace std;
using namespace TD_BEM;

const int MAX_CHARS_LINE=300;
const int MAX_TOKEN_LINE=100;
const char* const DELIMITER = " \t=;,";

InputT::InputT()
{

	fFormulation=1;
	fInterpolation=1;
	fNumStep=1;
	fDT=0;
	fNQ=6;


	fE=1;
	fNu=0;
	fDensity=1;

	fSD=3;
	fIfExt=0;
	fNND=0;
	fNEL=0;
	fEType=1;
	fENND=4;
	fIEN=NULL;
	fCoords=NULL;

	fNND_Aux=0;
	fNEL_Aux=0;
	fEType_Aux=1;
	fENND_Aux=4;

	fCoords_Aux=NULL;
	fIEN_Aux=NULL;

	fUBC_num=0;
	fUBC_DOFs=NULL;
	fUBC_Values=NULL;

	fFBC_num=0;
	fFBC_DOFs=NULL;
	fFBC_Values=NULL;

	fNumInterface=0;
	fInterface=NULL;
}

InputT::~InputT()
{
	delete[] fIEN;
	delete[] fCoords;

	if(!fIEN_Aux) delete[] fIEN_Aux;
	if(!fCoords_Aux) delete[] fCoords_Aux;

	delete[] fUBC_DOFs;
	delete[] fUBC_Values;

	delete[] fFBC_DOFs;
	delete[] fFBC_Values;

	if(!fInterface) delete[] fInterface;

}

void InputT::ReadInput(const char* fileName)
{
	cout << "reading " << fileName << endl;

	ifstream fin;
	fin.open(fileName);

	if(!fin.good())
	{
		cout << fileName << " not found \n";
		throw "file not found";
	}

	int line=0;

	while(!fin.eof())
	{
		line++;

		char buff[MAX_CHARS_LINE];
		char buff_temp[MAX_CHARS_LINE];

		fin.getline(buff, MAX_CHARS_LINE);

		PetscStrcpy(buff_temp, buff);
		if(ReadAlgorithm(buff_temp, line)) continue;

		PetscStrcpy(buff_temp, buff);
		if(ReadMaterial(buff_temp, line)) continue;

		PetscStrcpy(buff_temp, buff);
		if(ReadGeometry(fin, buff_temp, line)) continue;


		if(fFormulation==3)
		{
			PetscStrcpy(buff_temp, buff);
			if(ReadGeometry_Aux(fin, buff_temp, line)) continue;
		}


		PetscStrcpy(buff_temp, buff);
		if(ReadBCs(fin, buff_temp, line)) continue;



		char* token[MAX_TOKEN_LINE]={};
		token[0]=strtok(buff, DELIMITER);

		if(!token[0])	continue;

		cout << "key word " << token[0] << endl;
		throw  "InputT::ReadInput, unaccepted key word";

	}

	fin.close();

	cout << "model reading over\n";

}


bool InputT::ReadAlgorithm(char* buff, int line)
{
	char* token[MAX_TOKEN_LINE]={};

	token[0]=strtok(buff, DELIMITER);

	PetscBool flg;

	/************************************************
	 ** Read formulation to use
	 ************************************************/
	PetscStrcmp(token[0], "formulation", &flg);
	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "if_ext is not set properly";
		}
		else
		{
			fFormulation=atoi(token[1]);
			cout<< "fFormulation=" << fFormulation << endl;
		}

		return 1;
	}

	/************************************************
	 ** Read number of steps
	 ************************************************/
	PetscStrcmp(token[0], "numstep", &flg);

	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "number of step to run is not set properly";
		}
		else
		{
			fNumStep=atoi(token[1]);
			cout<< "fNumStep=" << fNumStep << endl;
		}

		return 1;
	}

	/************************************************
	 ** Read dt
	 ************************************************/
	PetscStrcmp(token[0], "dt", &flg);

	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "time step size is not set properly";
		}
		else
		{
			fDT=atof(token[1]);
			cout<< "fDT=" << fDT << endl;
		}

		return 1;
	}

	/************************************************
	 ** Read interpolation
	 ************************************************/
	PetscStrcmp(token[0], "interpolation", &flg);

	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "interpolation is not set properly";
		}
		else
		{
			fInterpolation=atof(token[1]);
			cout<< "fInterpolation=" << fInterpolation << endl;
		}

		return 1;
	}

	/************************************************
	 ** Read interpolation
	 ************************************************/
	PetscStrcmp(token[0], "nq", &flg);

	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "number of Gaussian points (nq) is not set properly";
		}
		else
		{
			fNQ=atof(token[1]);
			cout<< "fNQ=" << fNQ << endl;
		}

		return 1;
	}

	return 0;
}

bool InputT::ReadMaterial(char* buff, int line)
{
	char* token[MAX_TOKEN_LINE]={};

	token[0]=strtok(buff, DELIMITER);

	PetscBool flg;

	/************************************************
	 ** Read modulus
	 ************************************************/
	PetscStrcmp(token[0], "E", &flg);

	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "Young's modulus is not set properly";
		}
		else
		{
			fE=atof(token[1]);
			cout<< "fE=" << fE << endl;
		}

		return 1;
	}

	/************************************************
	 ** Read Poisson Ratio
	 ************************************************/
	PetscStrcmp(token[0], "nu", &flg);

	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "possion's ratio is not set properly";
		}
		else
		{
			fNu=atof(token[1]);
			cout<< "fNu=" << fNu << endl;
		}

		return 1;
	}

	/************************************************
	 ** Read material density
	 ************************************************/
	PetscStrcmp(token[0], "density", &flg);

	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "density is not set properly";
		}
		else
		{
			fDensity=atof(token[1]);
			cout<< "fDensity=" << fDensity << endl;
		}

		return 1;
	}

	return 0;
}


bool InputT::ReadGeometry(ifstream& fin, char* buff, int& line)
{
	char* token[MAX_TOKEN_LINE]={};

	token[0]=strtok(buff, DELIMITER);

	PetscBool flg;

	/************************************************
	 ** Read spatial dimension
	 ************************************************/
	PetscStrcmp(token[0], "sd", &flg);
	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "spatial dimension is not set properly";
		}
		else
		{
			fSD=atoi(token[1]);
			cout<< "fSD=" << fSD << endl;
		}

		return 1;
	}

	/************************************************
	 ** Read if external problem
	 ************************************************/
	PetscStrcmp(token[0], "if_ext", &flg);
	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "if_ext is not set properly";
		}
		else
		{
			fIfExt=atoi(token[1]);
			cout<< "fIfExt=" << fIfExt << endl;
		}

		return 1;
	}




	/************************************************
	 ** Read number of nodes
	 ************************************************/
	PetscStrcmp(token[0], "nnd", &flg);
	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "number of nodes is not set properly";
		}
		else
		{
			fNND=atoi(token[1]);
			cout<< "fNND=" << fNND << endl;
		}

		return 1;
	}

	/************************************************
	 ** Read type of element
	 ************************************************/
	PetscStrcmp(token[0], "EType", &flg);
	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "number of nodes is not set properly";
		}
		else
		{
			fEType=atoi(token[1]);
			if(fEType==1)
				fENND=4;
			else if(fEType==2)
				fENND=8;
			else
			{
				throw "InputT::ReadInput, unsupported element type";
			}
			cout<< "fEType=" << fEType << endl;
		}

		return 1;
	}

	/************************************************
	 ** Read number of elements
	 ************************************************/
	PetscStrcmp(token[0], "nel", &flg);
	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "number of nodes is not set properly";
		}
		else
		{
			fNEL=atoi(token[1]);
			cout<< "fNel=" << fNEL << endl;
		}

		return 1;
	}


	/************************************************
	 ** Read node coordinates
	 ************************************************/
	PetscStrcmp(token[0], "nodes", &flg);
	if(flg)
	{
		fCoords=new double[fSD*fNND];

		char bf[MAX_CHARS_LINE];

		for(int nid=0; nid<fNND; nid++)
		{
			fin.getline(bf, MAX_CHARS_LINE);

			token[0]=strtok(bf, DELIMITER);
			fCoords[nid*fSD+0]=atof(token[0]);

			token[1]=strtok(NULL, DELIMITER);
			fCoords[nid*fSD+1]=atof(token[1]);

			token[2]=strtok(NULL, DELIMITER);
			fCoords[nid*fSD+2]=atof(token[2]);
		}


		//MathOperationT::PrintMatrix(fNND, fSD, fCoords, "Nodes");

		return 1;
	}

	/************************************************
	 ** Read element connections
	 ************************************************/
	PetscStrcmp(token[0], "IEN", &flg);
	if(flg)
	{
		fIEN=new int[fNEL*fENND];

		char buff[MAX_CHARS_LINE];

		for(int eid=0; eid<fNEL; eid++)
		{
			fin.getline(buff, MAX_CHARS_LINE);
			line++;

			token[0]=strtok(buff, DELIMITER);
			fIEN[eid*fENND+0]=atoi(token[0])-1;

			for(int ni=1; ni<fENND; ni++)
			{
				token[ni]=strtok(NULL, DELIMITER);
				fIEN[eid*fENND+ni]=atoi(token[ni])-1;
			}
		}

		//MathOperationT::PrintMatrix(fNEL, fENND, fIEN, "fIEN");

		return 1;
	}

	return 0;
}



bool InputT::ReadGeometry_Aux(ifstream& fin, char* buff, int& line)
{
	char* token[MAX_TOKEN_LINE]={};

	token[0]=strtok(buff, DELIMITER);

	PetscBool flg;


	/************************************************
	 ** Read number of nodes
	 ************************************************/
	PetscStrcmp(token[0], "nnd_aux", &flg);
	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "number of nodes is not set properly";
		}
		else
		{
			fNND_Aux=atoi(token[1]);
			cout<< "fNND_Aux=" << fNND_Aux << endl;
		}

		return 1;
	}

	/************************************************
	 ** Read type of element
	 ************************************************/
	PetscStrcmp(token[0], "etype_aux", &flg);
	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "number of nodes is not set properly";
		}
		else
		{
			fEType_Aux=atoi(token[1]);
			if(fEType_Aux==1)
				fENND_Aux=4;
			else if(fEType==2)
				fENND_Aux=8;
			else
			{
				throw "InputT::ReadInput, unsupported element type";
			}
			cout<< "fEType_Aux=" << fEType_Aux << endl;
		}

		return 1;
	}

	/************************************************
	 ** Read number of elements
	 ************************************************/
	PetscStrcmp(token[0], "nel_aux", &flg);
	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "number of nodes is not set properly";
		}
		else
		{
			fNEL_Aux=atoi(token[1]);
			cout<< "fNel_Aux=" << fNEL_Aux << endl;
		}

		return 1;
	}


	/************************************************
	 ** Read node coordinates
	 ************************************************/
	PetscStrcmp(token[0], "nodes_aux", &flg);
	if(flg)
	{
		fCoords_Aux=new double[fSD*fNND_Aux];

		char bf[MAX_CHARS_LINE];

		for(int nid=0; nid<fNND_Aux; nid++)
		{
			fin.getline(bf, MAX_CHARS_LINE);

			token[0]=strtok(bf, DELIMITER);
			fCoords_Aux[nid*fSD+0]=atof(token[0]);

			token[1]=strtok(NULL, DELIMITER);
			fCoords_Aux[nid*fSD+1]=atof(token[1]);

			token[2]=strtok(NULL, DELIMITER);
			fCoords_Aux[nid*fSD+2]=atof(token[2]);
		}


		//MathOperationT::PrintMatrix(fNND_Aux, fSD, fCoords_Aux, "Nodes_aux");

		return 1;
	}

	/************************************************
	 ** Read element connections
	 ************************************************/
	PetscStrcmp(token[0], "IEN_aux", &flg);
	if(flg)
	{
		fIEN_Aux=new int[fNEL_Aux*fENND_Aux];

		char buff[MAX_CHARS_LINE];

		for(int eid=0; eid<fNEL_Aux; eid++)
		{
			fin.getline(buff, MAX_CHARS_LINE);
			line++;

			token[0]=strtok(buff, DELIMITER);
			fIEN_Aux[eid*fENND_Aux+0]=atoi(token[0])-1;

			for(int ni=1; ni<fENND_Aux; ni++)
			{
				token[ni]=strtok(NULL, DELIMITER);
				fIEN_Aux[eid*fENND_Aux+ni]=atoi(token[ni])-1;
			}
		}

		//MathOperationT::PrintMatrix(fNEL_Aux, fENND_Aux, fIEN_Aux, "fIEN_Aux");

		return 1;
	}

	return 0;
}




bool InputT::ReadBCs(std::ifstream& fin, char* buff, int& line)
{
	char* token[MAX_TOKEN_LINE]={};

	token[0]=strtok(buff, DELIMITER);

	PetscBool flg;

	/************************************************
	 ** Read number of displacement boundary
	 ************************************************/
	PetscStrcmp(token[0], "UBC_num", &flg);
	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "number of displacement BC is not set properly";
		}
		else
		{
			fUBC_num=atoi(token[1]);
			cout<< "fUBC_num=" << fUBC_num << endl;
		}

		return 1;
	}


	/************************************************
	 ** Read UBC dofs
	 ************************************************/
	PetscStrcmp(token[0], "UBC_dof", &flg);
	if(flg)
	{
		fUBC_DOFs=new int[fUBC_num];

		char bf[MAX_CHARS_LINE];

		int n_temp=0;

		while(n_temp<fUBC_num)
		{
			fin.getline(bf, MAX_CHARS_LINE);
			line++;

			token[0]=strtok(bf, DELIMITER);
			fUBC_DOFs[n_temp]=atoi(token[0])-1;
			n_temp++;

			int li=1;

			while(n_temp<fUBC_num)
			{
				token[li]=strtok(NULL, DELIMITER);

				if(!token[li])
					break;

				fUBC_DOFs[n_temp]=atoi(token[li])-1;
				li++;
				n_temp++;
			}
		}

		//MathOperationT::PrintVector(fUBC_num, fUBC_DOFs, "UBC_dof");

		return 1;
	}



	/************************************************
	 ** Read UBC values
	 ************************************************/
	PetscStrcmp(token[0], "UBC_values", &flg);
	if(flg)
	{
		fUBC_Values=new double[fUBC_num];

		char bf[MAX_CHARS_LINE];

		int n_temp=0;

		while(n_temp<fUBC_num)
		{
			fin.getline(bf, MAX_CHARS_LINE);
			line++;

			token[0]=strtok(bf, DELIMITER);
			fUBC_Values[n_temp]=atof(token[0]);
			n_temp++;

			int li=1;

			while(n_temp<fUBC_num)
			{
				token[li]=strtok(NULL, DELIMITER);

				if(!token[li])
					break;

				fUBC_Values[n_temp]=atof(token[li]);
				li++;
				n_temp++;
			}
		}

		//MathOperationT::PrintVector(fUBC_num, fUBC_Values, "UBC_Values");

		return 1;
	}



	/************************************************
	 ** Read number of traction boundary
	 ************************************************/
	PetscStrcmp(token[0], "FBC_num", &flg);
	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "number of displacement BC is not set properly";
		}
		else
		{
			fFBC_num=atoi(token[1]);
			cout<< "fFBC_num=" << fFBC_num << endl;
		}

		return 1;
	}


	/************************************************
	 ** Read FBC dofs
	 ************************************************/
	PetscStrcmp(token[0], "FBC_dof", &flg);
	if(flg)
	{
		fFBC_DOFs=new int[fFBC_num];

		char bf[MAX_CHARS_LINE];

		int n_temp=0;

		while(n_temp<fFBC_num)
		{
			fin.getline(bf, MAX_CHARS_LINE);
			line++;

			token[0]=strtok(bf, DELIMITER);
			fFBC_DOFs[n_temp]=atoi(token[0])-1;
			n_temp++;

			int li=1;

			while(n_temp<fFBC_num)
			{
				token[li]=strtok(NULL, DELIMITER);

				if(!token[li])
					break;

				fFBC_DOFs[n_temp]=atoi(token[li])-1;
				li++;
				n_temp++;
			}
		}

		//MathOperationT::PrintVector(fFBC_num, fFBC_DOFs, "FBC_dof");

		return 1;
	}



	/************************************************
	 ** Read FBC values
	 ************************************************/
	PetscStrcmp(token[0], "FBC_values", &flg);
	if(flg)
	{
		fFBC_Values=new double[fFBC_num];

		char bf[MAX_CHARS_LINE];

		int n_temp=0;

		while(n_temp<fFBC_num)
		{
			fin.getline(bf, MAX_CHARS_LINE);
			line++;

			token[0]=strtok(bf, DELIMITER);
			fFBC_Values[n_temp]=atof(token[0]);
			n_temp++;

			int li=1;

			while(n_temp<fFBC_num)
			{
				token[li]=strtok(NULL, DELIMITER);

				if(!token[li])
					break;

				fFBC_Values[n_temp]=atof(token[li]);
				li++;
				n_temp++;
			}
		}

		//MathOperationT::PrintVector(fFBC_num, fFBC_Values, "FBC_Values");

		return 1;
	}


	/************************************************
	 ** Read number of interface nodes
	 ************************************************/
	PetscStrcmp(token[0], "num_interface", &flg);
	if(flg)
	{
		token[1]=strtok(NULL, DELIMITER);

		if(!token[1])
		{
			cout<< "Error occurs as reading line " << line << endl;
			throw "number of interface nodes is not set properly";
		}
		else
		{
			fNumInterface=atoi(token[1]);
			cout<< "fNumInterface=" << fNumInterface << endl;
		}

		return 1;
	}







		/************************************************
			** Read inteface dofs
		************************************************/
		PetscStrcmp(token[0], "interface", &flg);
		if(flg)
		{
			fInterface=new int[fNumInterface];

				char bf[MAX_CHARS_LINE];

				int n_temp=0;

				while(n_temp<fNumInterface)
				{
					fin.getline(bf, MAX_CHARS_LINE);
					line++;

					token[0]=strtok(bf, DELIMITER);
					fInterface[n_temp]=atoi(token[0])-1;
					n_temp++;

					int li=1;

					while(n_temp<fNumInterface)
					{
						token[li]=strtok(NULL, DELIMITER);

						if(!token[li])
							break;

						fInterface[n_temp]=atoi(token[li])-1;
						li++;
						n_temp++;
					}
				}

				//MathOperationT::PrintVector(fNumInterface, fInterface, "fInterface");

				return 1;
			}


	return 0;
}


int InputT::GetFormulation()
{
	return fFormulation;
}

int InputT::GetNumStep(void)
{
	return fNumStep;
}

double InputT::GetDT(void)
{
	return fDT;
}

int InputT::GetTimeInterpolation(void)
{
	return fInterpolation;
}

int InputT::GetNQ(void)
{
	return fNQ;
}


int InputT::GetSD(void)
{
	return fSD;
}

int InputT::GetIfExt(void)
{
	return fIfExt;
}

int InputT::GetNND(void)
{
	return fNND;
}

int InputT::GetElementType(void)
{
	return fEType;
}

int InputT::GetENND(void)
{
	return fENND;
}

int InputT::GetNEL(void)
{
	return fNEL;
}

int* InputT::GetIEN(void)
{
	return fIEN;
}

double* InputT::GetCoords(void)
{
	return fCoords;
}


int InputT::GetNND_Aux(void)
{
	return fNND_Aux;
}

int InputT::GetElementType_Aux(void)
{
	return fEType_Aux;
}

int InputT::GetENND_Aux(void)
{
	return fENND_Aux;
}

int InputT::GetNEL_Aux(void)
{
	return fNEL_Aux;
}

int* InputT::GetIEN_Aux(void)
{
	return fIEN_Aux;
}

double* InputT::GetCoords_Aux(void)
{
	return fCoords_Aux;
}


void InputT::GetMaterialConstants(double& E, double& nu, double& rho)
{
	E=fE;
	nu=fNu;
	rho=fDensity;
}

void InputT::BoundaryCondition(int& UBC_num, int** UBC_DOFs, double** UBC_Values,
						   int& FBC_num, int** FBC_DOFs, double** FBC_Values)
{
    UBC_num=fUBC_num;
    FBC_num=fFBC_num;
    *UBC_DOFs=fUBC_DOFs;
    *UBC_Values=fUBC_Values;
    *FBC_DOFs=fFBC_DOFs;
    *FBC_Values=fFBC_Values;
}


void InputT::GetInterface(int& num_interface, int** interface)
{
	num_interface=fNumInterface;
	*interface=fInterface;
}


