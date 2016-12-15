/*
 * SolverT.cpp
 *
 *  Created on: Jun 2, 2015
 *      Author: xiaoyong
 */

#include "SolverT.h"
#include "MathOperationT.h"

#include <iostream>
#include <fstream>
#include <iomanip>

#include <string>
#include <sstream>
#include <algorithm>

using namespace std;
using namespace TD_BEM;
using namespace GreenFunction;

SolverT::SolverT(ModelManagerT* model)
{
	fModel=model;

	fComm=PETSC_COMM_WORLD;

	MPI_Comm_size(fComm, &fSize);
	MPI_Comm_rank(fComm, &fRank);

	LocalIndexing();

	fIerr=MatCreateDense(fComm, PETSC_DECIDE, PETSC_DECIDE, fGblNumDof, fGblNumDof, PETSC_NULL, &fG_Compute);
	fIerr=MatCreateDense(fComm, PETSC_DECIDE, PETSC_DECIDE, fGblNumDof, fGblNumDof, PETSC_NULL, &fH_Compute);
	fIerr=VecCreateMPI(fComm,  PETSC_DECIDE, fGblNumDof, &fLoad);
	fIerr=VecCreateMPI(fComm,  PETSC_DECIDE, fGblNumDof, &fRHS);

	//SetLoad();

	fIerr=KSPCreate(fComm, &fKsp);

	fMaxStep=10;
	fBeta=0;
	fDt=model->GetDT();
	fNumStep=model->GetNumStep();

    ComputeMaxStep();
    
	cout<<"Solver is constructed\n";
}

SolverT::~SolverT()
{
	delete[] fIndexRow;
	delete[] fIndexCol;

	/**************************************************************************
	 **
	 ** Release the memory
	 **
	 * ************************************************************************/
	for (int i=0; i<fNumStep+1; i++)
	{
		VecDestroy(&fDis_DB[i]);
		VecDestroy(&fTrac_DB[i]);
	}

	delete[] fDis_DB;
	delete[] fTrac_DB;


	for (int i=0; i<fMaxStep+1; i++)
	{
		MatDestroy(&fG_DB[i]);
		MatDestroy(&fH_DB[i]);
	}

	delete[] fG_DB;
	delete[] fH_DB;


	for (int i=0; i<3; i++)
	{
		MatDestroy(&fG_DB_bgn[i]);
		MatDestroy(&fH_DB_bgn[i]);
	}
	delete[] fG_DB_bgn;
	delete[] fH_DB_bgn;

	MatDestroy(&fG_Compute);
	MatDestroy(&fH_Compute);

	VecDestroy(&fLoad);
	VecDestroy(&fRHS);

	KSPDestroy(&fKsp);
}


void SolverT::SetTimeControl(double dt, int Numsteps)
{
	fDt=dt;
	fNumStep=Numsteps;

	/**********************************************************
	 **
	 ** Allocate field variables and matrices
	 **
	 **********************************************************/
	fDis_DB= new Vec[fNumStep+1];
	fTrac_DB= new Vec[fNumStep+1];

	Vec* temp_d=fDis_DB;
	Vec* temp_t=fTrac_DB;

	for (int i=0; i<fNumStep+1; i++)
	{
		fIerr=VecCreateMPI(PETSC_COMM_WORLD,  PETSC_DECIDE, fGblNumDof, temp_d);
		fIerr=VecCreateMPI(PETSC_COMM_WORLD,  PETSC_DECIDE, fGblNumDof, temp_t);
		VecSet(*temp_d, 0);
		VecSet(*temp_t, 0);
		temp_d++;
		temp_t++;
	}


	fG_DB=new Mat[fMaxStep+1];
	fH_DB=new Mat[fMaxStep+1];

	Mat* temp_G=fG_DB;
	Mat* temp_H=fH_DB;

	for (int i=0; i<fMaxStep+1; i++)
	{
		fIerr=MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, fGblNumDof, fGblNumDof, PETSC_NULL, temp_G);
		fIerr=MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, fGblNumDof, fGblNumDof, PETSC_NULL, temp_H);
		temp_G++;
		temp_H++;
	}

	fG_DB_bgn=new Mat[3];
	fH_DB_bgn=new Mat[3];

	temp_G=fG_DB_bgn;
	temp_H=fH_DB_bgn;

	for (int i=0; i<3; i++)
	{
		fIerr=MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, fGblNumDof, fGblNumDof, PETSC_NULL, temp_G);
		fIerr=MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, fGblNumDof, fGblNumDof, PETSC_NULL, temp_H);
		temp_G++;
		temp_H++;
	}
}


void SolverT::LocalIndexing()
{
	int sd=fModel->GetSD();
	int nnd=fModel->GetNND();
	fGblNumDof=sd*nnd;

	/**********************************************************
	 **
	 ** Create local indexing
	 **
	 **********************************************************/
	Vec NodeIndex;
	VecCreateMPI(fComm, PETSC_DECIDE, nnd, &NodeIndex);
	VecSet(NodeIndex, 0);

	VecGetOwnershipRange(NodeIndex, &fNodeLow, &fNodeHigh);


	fNodeHigh=fNodeHigh-1;

	fDofLow=fNodeLow*3;
	fDofHigh=fNodeHigh*3+2;

	fLocNumDof=(fNodeHigh-fNodeLow+1)*sd;
	fIndexRow=new int[fLocNumDof];
	for(int mi=fNodeLow*sd; mi<(fNodeHigh+1)*sd; mi++)
		fIndexRow[mi-fNodeLow*sd]=mi;

	fIndexCol=new int[fGblNumDof];
	for(int ni=0; ni<fGblNumDof; ni++)
		fIndexCol[ni]=ni;


	VecDestroy(&NodeIndex);
	//MathOperationT::PrintVector(fLocNumDof, fIndexRow, "indexRow");
	//MathOperationT::PrintVector(fGblNumDof, fIndexCol, "indexCol");
}

void SolverT::Initialize()
{
	//do nothing;
}


void SolverT::SetCurrStep(int curr_step)
{
	fCurrStep=curr_step;
}

/**
 * Determine the maximum number of G and H needed to be computed.
 * MaxStep=MaxDistance/Cs.
*/
void SolverT::ComputeMaxStep(void)
{
    double E, nu, rho;
    fModel->GetMaterialConstants(E, nu, rho);
    
    double Mu=E/(2.0*(1.0+nu));
    double Cs=sqrt(Mu/rho);
    
    double nnd=fModel->GetNND();
    const double* nodes=fModel->GetCoords();
    
    double max_dis=0; //maximum distance between any two nodes
    
    double node_1[3], node_2[3];
    for (int i=0; i<nnd-1; i++) {
        
        node_1[0]=nodes[i*3+0];
        node_1[1]=nodes[i*3+1];
        node_1[2]=nodes[i*3+2];
        
        for (int j=i+1; j<nnd; j++) {
            node_2[0]=nodes[j*3+0];
            node_2[1]=nodes[j*3+1];
            node_2[2]=nodes[j*3+2];
            
            double temp=0;
            for (int d=0; d<3; d++) {
                temp+=pow(node_1[d]-node_2[d], 2);
            }
            temp=sqrt(temp);
            
            if (temp > max_dis) {
                max_dis=temp;
            }
        }
    }
    
    
    fMaxStep=ceil(max_dis/(Cs*fDt))+1;
}

int SolverT::GetMaxStep()
{
    return fMaxStep;
}


void SolverT::UpdateGHLinear(const double* G_0, const double* G_1, const double* H_0, const double* H_1)
{
	//do nothing;
}



void SolverT::ExchangeColumns()
{
	int UBC_num;
	int* UBC_DOFs;
	double* UBC_Values;
	int FBC_num;
	int* FBC_DOFs;
	double* FBC_Values;

	fModel->BoundaryCondition_Global(UBC_num, &UBC_DOFs, &UBC_Values, FBC_num, &FBC_DOFs, &FBC_Values);

	Mat G_temp, H_temp;
	fIerr=MatConvert(fG_Compute, MATSAME, MAT_INITIAL_MATRIX, &G_temp);
	fIerr=MatConvert(fH_Compute, MATSAME, MAT_INITIAL_MATRIX, &H_temp);

	//rearrange the matrices
	Vec column_vec;
	fIerr=VecCreateMPI(fComm, PETSC_DECIDE, fGblNumDof, &column_vec);
	int column_low, column_high;
	int matrix_low, matrix_high;
	VecGetOwnershipRange(column_vec, &column_low, &column_high);

	MatGetOwnershipRange(H_temp, &matrix_low, &matrix_high);
	if(column_low != matrix_low || column_high != matrix_high)
	{
		throw "SolverT::solver, the matrix and column are partitioned in different manners";
	}

	MatGetOwnershipRange(G_temp, &matrix_low, &matrix_high);
	if(column_low != matrix_low || column_high != matrix_high)
	{
		throw "SolverT::solver, the matrix and column are partitioned in different manners";
	}

	int column_length=column_high-column_low;
//cout<< "rank" << fRank << "column_length is " << column_length <<endl;
	int* LocalIndexing=new int[column_length];
	for(int ii=column_low; ii<column_high; ii++)
		LocalIndexing[ii-column_low]=ii;

//MathOperationT::PrintVector(column_length, LocalIndexing, "localIndexing");

	for(int u_i=0; u_i<UBC_num; u_i++)
	{
		int u_dof=UBC_DOFs[u_i];
//cout << "u_dof is " << u_dof << endl;
		fIerr=MatGetColumnVector(H_temp, column_vec, u_dof);
		fIerr=VecScale(column_vec, -1);
		double* temp;
		fIerr=VecGetArray(column_vec, &temp);

		fIerr=MatSetValues(fG_Compute, column_length, LocalIndexing, 1, &u_dof, temp, INSERT_VALUES);

		fIerr=VecRestoreArray(column_vec, &temp);


		fIerr=MatGetColumnVector(G_temp, column_vec, u_dof);
		fIerr=VecScale(column_vec, -1);
		fIerr=VecGetArray(column_vec, &temp);

		fIerr=MatSetValues(fH_Compute, column_length, LocalIndexing, 1, &u_dof, temp, INSERT_VALUES);

		fIerr=VecRestoreArray(column_vec, &temp);

		//fIerr=VecView(column_vec, PETSC_VIEWER_STDOUT_WORLD);
	}

	fIerr=MatAssemblyBegin(fG_Compute, MAT_FINAL_ASSEMBLY);
	fIerr=MatAssemblyEnd(fG_Compute, MAT_FINAL_ASSEMBLY);
	fIerr=MatAssemblyBegin(fH_Compute, MAT_FINAL_ASSEMBLY);
	fIerr=MatAssemblyEnd(fH_Compute, MAT_FINAL_ASSEMBLY);

	fIerr=MatDestroy(&G_temp);
	fIerr=MatDestroy(&H_temp);
	fIerr=VecDestroy(&column_vec);

	delete[] LocalIndexing;
}

void SolverT::SetLoad()
{
	int UBC_num;
	int* UBC_DOFs;
	double* UBC_Values;
	int FBC_num;
	int* FBC_DOFs;
	double* FBC_Values;

	fModel->BoundaryCondition_Local(UBC_num, &UBC_DOFs, &UBC_Values, FBC_num, &FBC_DOFs, &FBC_Values);

	//****************************************
	double curr_time=(fCurrStep-1)*fDt;
    double amp=LoadTimeVariation(curr_time);
    
	double* FBC_Values_temp=new double[FBC_num];
	MathOperationT::MemCopy(FBC_num, FBC_Values, FBC_Values_temp);

	MathOperationT::VecScale(FBC_num, FBC_Values_temp, amp);

	//MathOperationT::PrintVector(FBC_num, FBC_Values_temp, "FBC_temp");

	VecSet(fLoad, 0);

	fIerr=VecSetValues(fLoad, FBC_num, FBC_DOFs, FBC_Values_temp, INSERT_VALUES);
	fIerr=VecSetValues(fLoad, UBC_num, UBC_DOFs, UBC_Values, INSERT_VALUES);

	fIerr=VecAssemblyBegin(fLoad);
	fIerr=VecAssemblyEnd(fLoad);

	//VecView(fLoad, PETSC_VIEWER_STDOUT_WORLD);

	delete[] FBC_Values_temp;
}

double SolverT::LoadTimeVariation(double t)
{
    double amp=1.0;
    
    //half-sine load
    double period=0.02; //period of the sin function
    
    if (t<=period/2.0) {
        amp=sin(2.0*M_PI*t/period);
    }
    
    //Heaviside load
    amp=1.0;
    
    return amp;
}



void SolverT::GetNodeLowHigh(int* low, int* high)
{
	*low=fNodeLow;
	*high=fNodeHigh;
}

void SolverT::GetDofLowHigh(int* low, int* high)
{
	*low=fDofLow;
	*high=fDofHigh;
}


void SolverT::OutPutTecPlot()
{
	ofstream myfile;

	Vec dis_collection;
	VecScatter ctx;

	VecScatterCreateToAll(fDis_DB[fCurrStep-1],&ctx,&dis_collection);
	VecScatterBegin(ctx,fDis_DB[fCurrStep-1],dis_collection,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(ctx,fDis_DB[fCurrStep-1],dis_collection,INSERT_VALUES,SCATTER_FORWARD);

	double* temp;
	fIerr=VecGetArray(dis_collection, &temp);

	int nnd=fModel->GetNND();
	int nel=fModel->GetNEL();

	const int* IEN=fModel->GetIEN();
	const double* Coords=fModel->GetCoords();


	if(fRank==0)
	{
		if (fCurrStep==1)
		{
			myfile.open("BEM_result.tecplot", fstream::out|fstream::trunc);

			/* write the Head */
			myfile<<"TITLE = \"BEM results\" " << endl;
			myfile<<"VARIABLES= \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\" " <<endl;
		}
		else
		{
			myfile.open("BEM_result.tecplot", fstream::out|fstream::app);
		}


		myfile<<"ZONE  ";
		myfile<<"N=" <<nnd << "," ;
		myfile<<" E=" <<nel << "," << "  DATAPACKING=POINT,";
		myfile<<" ZONETYPE=FEQUADRILATERAL" << endl;
		myfile<<"\t STRANDID=1, " << "SOLUTIONTIME=" << fDt*(fCurrStep-1) << endl;

		myfile<<setprecision(6);


		for(int ni=0; ni<nnd; ni++)
		{
			myfile<< setw(10)<< Coords[3*ni] <<setw(10) <<Coords[3*ni+1] <<setw(10) <<Coords[3*ni+2];
			myfile<< setw(15)<< temp[3*ni] <<setw(15) <<temp[3*ni+1] <<setw(15) <<temp[3*ni+2];
			myfile<<endl;
		}

		for(int ei=0; ei<nel; ei++)
		{
			myfile<< setw(8)<< IEN[4*ei]+1 <<setw(8) <<IEN[4*ei+1]+1 <<setw(8) <<IEN[4*ei+2]+1<<setw(8) <<IEN[4*ei+3]+1;
			myfile<<endl;
		}

		myfile.close();

	}

	fIerr=VecRestoreArray(dis_collection, &temp);

	VecScatterDestroy(&ctx);
	VecDestroy(&dis_collection);
}


void SolverT::OutPutVTK()
{
    int ele_type=fModel->GetElementType();
    
    if (ele_type==1) {
        OutPutVTK_Linear();
    }else if (ele_type==2){
        OutPutVTK_Quadratic();
    }else{
        throw "SolverT::OutPutVTK, unsupported element type";
    }
}


void SolverT::OutPutVTK_Linear()
{
	ofstream myfile;

	Vec dis_collection;
	VecScatter ctx;

	VecScatterCreateToAll(fDis_DB[fCurrStep-1],&ctx,&dis_collection);
	VecScatterBegin(ctx,fDis_DB[fCurrStep-1],dis_collection,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(ctx,fDis_DB[fCurrStep-1],dis_collection,INSERT_VALUES,SCATTER_FORWARD);

	int nnd=fModel->GetNND();
	int nel=fModel->GetNEL();

	const int* IEN=fModel->GetIEN();
	const double* Coords=fModel->GetCoords();


	if(fRank==0)
	{
		string filename="BEM_result_";

		ostringstream oo;
		oo<<fCurrStep;
		filename+=oo.str();

		filename+=".vtk";

		const char* name=filename.c_str();

		myfile.open(name, fstream::out|fstream::trunc);

		// write the Head
		myfile<<"# vtk DataFile Version 3.0 " << endl;
		myfile<<"Boundary element result"<<endl;
		myfile<<"ASCII"<<endl;
		myfile<<"DATASET UNSTRUCTURED_GRID " <<endl;
		myfile<<"POINTS " << nnd << " float"<<endl;

		myfile<<setprecision(6);

		for(int ni=0; ni<nnd; ni++)
		{
			myfile<< setw(15)<<Coords[3*ni] <<setw(15) <<Coords[3*ni+1] <<setw(15) <<Coords[3*ni+2];
			//myfile<< setw(15)<< temp[3*ni] <<setw(15) <<temp[3*ni+1] <<setw(15) <<temp[3*ni+2];
			myfile<<endl;
		}

		myfile<<"CELLS  " << nel << " " << nel*5 <<endl;

		for(int ei=0; ei<nel; ei++)
		{
			myfile<< setw(5)<<4 << setw(5)<< IEN[4*ei] << setw(5) <<IEN[4*ei+1] << setw(5) <<IEN[4*ei+2]<< setw(5) <<IEN[4*ei+3];
			myfile<<endl;
		}

		myfile<<"CELL_TYPES "<<nel<<endl;
		for(int ei=0; ei<nel; ei++)
		{
			myfile<<9<<endl;
		}


		double* temp;
		fIerr=VecGetArray(dis_collection, &temp);

		myfile<<"POINT_DATA " << nnd << endl;
		myfile<<"VECTORS " << "Displacement " << "double"<<endl;
		for(int ni=0; ni<nnd; ni++)
		{
     		myfile<< setw(15)<< temp[3*ni] <<setw(15) <<temp[3*ni+1] <<setw(15) <<temp[3*ni+2];
			myfile<<endl;
		}

		fIerr=VecRestoreArray(dis_collection, &temp);


		myfile.close();
	}
    
    VecScatterDestroy(&ctx);
    VecDestroy(&dis_collection);
}

void SolverT::OutPutVTK_Quadratic()
{
    ofstream myfile;
    
    Vec dis_collection;
    VecScatter ctx;
    
    VecScatterCreateToAll(fDis_DB[fCurrStep-1],&ctx,&dis_collection);
    VecScatterBegin(ctx,fDis_DB[fCurrStep-1],dis_collection,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx,fDis_DB[fCurrStep-1],dis_collection,INSERT_VALUES,SCATTER_FORWARD);
    
    int nnd=fModel->GetNND();
    int nel=fModel->GetNEL();
    
    const int* IEN=fModel->GetIEN();
    const double* Coords=fModel->GetCoords();
    
    if(fRank==0)
    {
        string filename="BEM_result_";
        
        ostringstream oo;
        oo<<fCurrStep;
        filename+=oo.str();
        
        filename+=".vtk";
        
        const char* name=filename.c_str();
        
        myfile.open(name, fstream::out|fstream::trunc);
        
        // write the Head
        myfile<<"# vtk DataFile Version 3.0 " << endl;
        myfile<<"Boundary element result"<<endl;
        myfile<<"ASCII"<<endl;
        myfile<<"DATASET UNSTRUCTURED_GRID " <<endl;
        myfile<<"POINTS " << nnd << " float"<<endl;
        
        myfile<<setprecision(6);
        
        for(int ni=0; ni<nnd; ni++)
        {
            myfile<< setw(15)<<Coords[3*ni] <<setw(15) <<Coords[3*ni+1] <<setw(15) <<Coords[3*ni+2];
            myfile<<endl;
        }
        
        myfile<<"CELLS  " << nel << " " << nel*9 <<endl;
        
        for(int ei=0; ei<nel; ei++)
        {
            myfile<< setw(5)<<8 << setw(5)<< IEN[8*ei] << setw(5) <<IEN[8*ei+1] << setw(5) <<IEN[8*ei+2]<< setw(5) <<IEN[8*ei+3] << setw(5) <<IEN[8*ei+4] << setw(5) <<IEN[8*ei+5] << setw(5) <<IEN[8*ei+6] <<
                              setw(5) <<IEN[8*ei+7] << endl;
            myfile<<endl;
        }
        
        myfile<<"CELL_TYPES "<<nel<<endl;
        for(int ei=0; ei<nel; ei++)
        {
            myfile<<23<<endl;
        }
        
        
        double* temp;
        fIerr=VecGetArray(dis_collection, &temp);
        
        myfile<<"POINT_DATA " << nnd << endl;
        myfile<<"VECTORS " << "Displacement " << "double"<<endl;
        for(int ni=0; ni<nnd; ni++)
        {
            myfile<< setw(15)<< temp[3*ni] <<setw(15) <<temp[3*ni+1] <<setw(15) <<temp[3*ni+2];
            myfile<<endl;
        }
        
        fIerr=VecRestoreArray(dis_collection, &temp);
        
        
        myfile.close();
        
    }

	VecScatterDestroy(&ctx);
	VecDestroy(&dis_collection);
}


/***********************************************************************************************
 * Objective : This function write one global G or H matrices for the later usage,  namely,
 *             the eigenvalue analysis or other computations
 * Step 1: each processor write its local G or H
 * Step 2: The root processor read all the G or H, emerge them and delete them
 ***********************************************************************************************/
void SolverT::WriteGH(const char* name_prefix, int step, int num_row, int num_column, double* data, int comm_size, int comm_rank)
{
    //All processor write the matrices simutaneously
    string file_name;
    ostringstream converter;
    converter << name_prefix << step <<"_p_"<< comm_rank<<".txt";
    file_name=converter.str();
    
    ofstream file_writer;
    file_writer.open(file_name.c_str(), ios_base::out);
    file_writer << num_row << " "<<num_column <<endl;
    for (int i=0; i<num_row; i++) {
        for (int j=0; j<num_column; j++) {
            file_writer<< setw(15) << setprecision(5) << std::scientific<<data[i*fGblNumDof+j];
        }
        file_writer<<endl;
    }
    file_writer.close();
    
    //rank 0 merges all the matrices
    string global_file;
    ostringstream global_converter;
    global_converter << name_prefix << step <<".txt";
    global_file=global_converter.str();
    
    ofstream global_writer;
    global_writer.open(global_file.c_str(), ios_base::out);

    MPI_Barrier(fComm);
    
    if (comm_rank==0) {
        for (int ri=0; ri<comm_size; ri++) {
            string read_name;
            ostringstream read_converter;
            read_converter << name_prefix << step << "_p_" << ri <<".txt";
            read_name=read_converter.str();
            
            ifstream file_reader;
            file_reader.open(read_name.c_str(), ios_base::in);
            
            int nr, nc; //number of rows, columns
            file_reader >> nr;
            file_reader >> nc;
            
            for (int i=0; i<nr; i++) {
                for (int j=0; j<nc; j++) {
                    double temp;
                    file_reader>>temp;
                    global_writer<< setw(15) << setprecision(5) << std::scientific<<temp;
                }
                global_writer << endl;
            }
            
            file_reader.close();
            
            std::remove(read_name.c_str());
        }
    }
    
    global_writer.close();
    
    MPI_Barrier(fComm);
}





















