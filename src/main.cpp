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

void Model_Compability_check(ModelManagerT*, ModelManagerT*);

int main(int argc, char **args)
{

	PetscInitialize(&argc, &args, (char*)0, NULL);

	int size, rank;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	/**********************************************************
	 **
	 ** //get geometry information
	 **
	 **********************************************************/
	if(argc <=1)
	{
		try
		{
			throw "input file not specified, try ../geometry/TD_BEM model.txt";
		}
		catch(const char* msg){
			cerr<<msg<<endl;
			return 0;
		}
	}
	else if(argc==2) //case of single region (one layer)
	{
		ModelManagerT* model=new ModelManagerT;

		model->ReadInput(args[1]);

		 // Set time controls
		int num_step=model->GetNumStep();
		double time_step=model->GetDT();

		SolverT* solver;

		switch(model->GetTimeInterpolation())
		{
			case 1:
				//TODO
				//solver=new SolverConstantT(model);
				break;
			case 2:
				solver=new SolverHeavisideLinearT(model);
				break;
			case 3:
				//solver=new SolverConstantT(model);
				break;
			case 4:
				solver=new SolverStepLinearT(model);
				break;
			case 5:
				solver=new SolverStepLinearT(model);
				break;
			default:
				throw "main::create solver, unsupported interpolation type";
		}


		int low, high;
		solver->GetDofLowHigh(&low, &high);
		model->Distributed_BC(low, high);
	
		double beta=0.1;
		int temp_max=ceil(50/beta);
		int MaxStep=min(temp_max, num_step);
		solver->SetMaxStep(MaxStep);

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
		int interpolation=model->GetTimeInterpolation();

		try{
			ElastoDynamicsT* ED=new ElastoDynamicsT(model);

			int node_low, node_high;
			solver->GetNodeLowHigh(&node_low, &node_high);

			ED->SetLowerUpper(node_low, node_high);

			if(model->GetFormulation()==3)	ED->FormGeomBMatrix();

			double t_curr, t_1;

			for (int curr_step=1; curr_step<=num_step; curr_step++)
			{
				t_curr=curr_step*time_step;
				t_1=time_step;


				if(interpolation==1)
				{
					//waitting for implementation
				}
				else if(interpolation==2)
				{
					if(curr_step<=MaxStep)
					{
						ED->G_H_Driver(t_curr, t_1);
					}

					double* G_0;
					double* G_1;
					double* H_0;
					double* H_1;

					ED->GetGHLinear(&G_0, &G_1, &H_0, &H_1);


					solver->SetCurrStep(curr_step);

					solver->UpdateGHLinear(G_0, G_1, H_0, H_1);
				}
				else if(interpolation==4 || interpolation==5)
				{
					if(curr_step<=MaxStep)
					{
						ED->G_H_Driver(t_curr, t_1);
					}

					double* G_0;
					double* G_1;
					double* H_0;
					double* H_1;

					ED->GetGHLinear(&G_0, &G_1, &H_0, &H_1);

					//MathOperationT::PrintMatrix(9, 693, G_0, "G_0");
					//MathOperationT::PrintMatrix(9, 693, G_1, "G_1");
					//MathOperationT::PrintMatrix(9, 693, H_0, "H_0");
					//MathOperationT::PrintMatrix(9, 693, H_1, "H_1");


					solver->SetCurrStep(curr_step);

					solver->UpdateGHLinear(G_0, G_1, H_0, H_1);
				}
				else
				{
					throw "main:unsupported solving procedure";
				}

				solver->SetLoad();

				solver->Solve();

				solver->RecordResult();
			}

			delete ED;

		}catch(const char* msg){
			cerr<<msg<<endl;
			return 0;
		}

		delete solver;
		delete model;

	}//end case of single region
	else if(argc==3) //case of two regions
	{
		ModelManagerT* model_1=new ModelManagerT;
		ModelManagerT* model_2=new ModelManagerT;

		model_1->ReadInput(args[1]);
		model_2->ReadInput(args[2]);

		Model_Compability_check(model_1, model_2);

		/**********************************************************
		**
		** Set time controls
		**
		**********************************************************/
		int num_step=model_1->GetNumStep();
		double time_step=model_1->GetDT();

		SolverMRT* solver;


		switch(model_1->GetTimeInterpolation())
		{
			case 1:
				//solver=new SolverConstantT(model);
				break;
			case 2:
				//
				break;
			case 3:
				//solver=new SolverConstantT(model);
				break;
			case 4:
				solver=new SolverMRT(model_1, model_2);
				break;
			default:
				throw "main::create solver, unsupported interpolation type";
		}


		//int low_1, high_1, low_2, high_2;
		//solver->GetDofLowHigh(&low_1, &high_1, &low_2, &high_2);
		//model_1->Distributed_BC(low_1, high_1);
		//model_2->Distributed_BC(low_2, high_2);

		int MaxStep=min(num_step, 45);
		solver->SetMaxStep(MaxStep);

		solver->SetTimeControl(time_step, num_step);

		/*********************************************************************
		 **
		 **  Initialize the traction for the zeroth time step
		 **
		 **  In case of traction being initially zero, this step is skipped
		 *********************************************************************/
		//solver->Initialize();

		/***********************************************************************
		 **
		 **  Computation
		 **
		 **********************************************************************/
		int interpolation=model_1->GetTimeInterpolation();

		try{
				ElastoDynamicsT* ED_1=new ElastoDynamicsT(model_1);
				ElastoDynamicsT* ED_2=new ElastoDynamicsT(model_2);

				int node_low_1, node_high_1, node_low_2, node_high_2;
				solver->GetNodeLowHigh(&node_low_1, &node_high_1, &node_low_2, &node_high_2);

				ED_1->SetLowerUpper(node_low_1, node_high_1);
				ED_2->SetLowerUpper(node_low_2, node_high_2);

				if(model_1->GetFormulation()==3)
				{
					ED_1->FormGeomBMatrix();
					ED_2->FormGeomBMatrix();
				}

				double t_curr, t_1;

				for (int curr_step=1; curr_step<=num_step; curr_step++)
				{
					t_curr=curr_step*time_step;
					t_1=time_step;

					if(interpolation==1)
					{
						//waitting for implementation
					}
					else if(interpolation==2)
					{
						//waitting for implementation
					}
					else if(interpolation==4)
					{
						if(curr_step<=MaxStep)
						{
							ED_1->G_H_Driver(t_curr, t_1);
							ED_2->G_H_Driver(t_curr, t_1);
						}

						double* G_0, *g_0;
						double* G_1, *g_1;
						double* H_0, *h_0;
						double* H_1, *h_1;

						ED_1->GetGHLinear(&G_0, &G_1, &H_0, &H_1);
						ED_2->GetGHLinear(&g_0, &g_1, &h_0, &h_1);
						//MathOperationT::PrintMatrix(9, 693, G_0, "G_0");

						solver->SetCurrStep(curr_step);

						solver->UpdateGHLinear(G_0, G_1, H_0, H_1, g_0, g_1, h_0, h_1);
					}
					else
					{
						throw "main:unsupported solving procedure";
					} //end if

					solver->SetLoad();

					solver->Solve();

					solver->RecordResult();

				} //end for loop

				delete ED_1;
				delete ED_2;

			}catch(const char* msg){
				cerr<<msg<<endl;
				return 0;
			} //end try.



		delete model_1;
		delete model_2;

		delete solver;

	} //end if arg==3


	PetscFinalize();

cout << "The computation terminated properly" << endl;

return 0;

}





void Model_Compability_check(ModelManagerT* model_1, ModelManagerT* model_2)
{
	try{
		if(model_1->GetDT() != model_2->GetDT())
			throw "Time step size for the two regions doesn't match \n";

		if(model_1->GetFormulation() != model_2->GetFormulation())
			throw "Formulation doesn't match\n";

		if(model_1->GetTimeInterpolation() != model_2->GetTimeInterpolation())
		{
			throw "Time interpolation mismatch\n";
		}

		int num_interface_1=model_1->GetNumInterface();
		int num_interface_2=model_2->GetNumInterface();
		if( num_interface_1 != num_interface_2)
		{
			cout << "Number of Interface nodes in region 1 is " << num_interface_1 << endl;
			cout << "Number of Interface nodes in region 2 is " << num_interface_2 << endl;
			throw "Interface doesn't match\n";
		}
	}
	catch(const char* msg){
		cerr<<msg<<endl;
	}
}

