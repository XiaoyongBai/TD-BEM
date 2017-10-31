#include <iostream>
#include "petscksp.h"
#include "petscmat.h"
#include "petscsys.h"
#include "mpi.h"

#include "GreenFunctionT.h"
#include "ModelManagerT.h"
#include "ElastoDynamicsT.h"
#include "MathOperationT.h"
#include "SolverT.h"
#include "SolverHeavisideLinearT.h"
#include "SolverStepLinearT.h"
#include "SolverMRT.h"
#include "cstring"
#include "fstream"

using namespace std;
using namespace TD_BEM;


int main(int argc, char **args)
{

	PetscInitialize(&argc, &args, (char*)0, NULL);

	int size, rank;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	if(argc <=1)
	{
		try{
			throw "Input file not specified, try ./TD_BEM input.txt" ;
		}
		catch(const char* msg){
			cerr<<msg<<endl;
			return 0;
		}
	}
	else if(argc==2)
	{
		ModelManagerT* model=new ModelManagerT;

		model->ReadInput(args[1]);

		SolverT* solver;
        	solver=new SolverStepLinearT(model);

		int dof_low, dof_high;
		solver->GetDofLowHigh(&dof_low, &dof_high);
		model->Distributed_BC(dof_low, dof_high);
	
        	int num_step=model->GetNumStep();
        	double time_step=model->GetDT();
        
		solver->SetTimeControl(time_step, num_step);

		/*********************************************************************
		 **
		 **  Initialize the traction for the zeroth time step
		 **
		 *********************************************************************/
		solver->Initialize();

		/***********************************************************************
		 **
		 **  Computation
		 **
		 **********************************************************************/
		try{
			ElastoDynamicsT* ED=new ElastoDynamicsT(model);

			int node_low, node_high;
			solver->GetNodeLowHigh(&node_low, &node_high);

			ED->SetLowerUpper(node_low, node_high);

			double t_curr, t_1;

            		double a1, a2;
            		ifstream alpha_reader;
            		alpha_reader.open("alpha.txt");
            		alpha_reader>>a1>>a2;
            		alpha_reader.close();
            
			if (rank==0) {
				cout <<"***********************************" <<endl;
				cout <<"alpha_1="<<a1<<"  alpha_2="<<a2<<endl;
				cout <<"***********************************" <<endl;
		 	}
            
			for (int curr_step=1; curr_step<=num_step; curr_step++)
			{
				t_curr=curr_step*time_step;
				t_1=time_step;

                		int MaxStep=solver->GetMaxStep();
                		if(curr_step <=  MaxStep)
				{
					ED->G_H_Driver(t_curr, t_1);
				}

				double* G_0;
				double* G_1;
				double* H_0;
				double* H_1;

				ED->GetGHLinear(&G_0, &G_1, &H_0, &H_1);

				//write the matrices
				if(curr_step<=MaxStep)
				{
				    int num_row = dof_high-dof_low+1;
				    int num_column = 3*(model->GetNND());

				    //solver->WriteGH("G1_step", curr_step-1, num_row, num_column, G_1, size, rank);
				    //solver->WriteGH("G2_step", curr_step-1, num_row, num_column, G_0, size, rank);
				    //solver->WriteGH("H1_step", curr_step-1, num_row, num_column, H_1, size, rank);
				    //solver->WriteGH("H2_step", curr_step-1, num_row, num_column, H_0, size, rank);
				}

				solver->SetCurrStep(curr_step);

						solver->UpdateGHLinear(G_0, G_1, H_0, H_1);

				solver->SetLoad();
                
				solver->Solve(a1, a2);

				solver->RecordResult();
			}

			delete ED;

		}catch(const char* msg){
			cerr<<msg<<endl;
			return 0;
		}

		delete solver;
		delete model;

	}


	PetscFinalize();

    	cout << "The computation terminated properly" << endl;

	return 0;

}
	



