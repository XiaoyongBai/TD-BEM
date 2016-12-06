#include "math.h"
#include "MathOperationT.h"

using namespace std;
using namespace GreenFunction;

double MathOperationT::Determinant(int MD, double* matrix)
{

	if(MD==2)
	{
		return matrix[0]*matrix[3] - matrix[1]*matrix[2];
	}

	else if(MD==3)
	{
		return matrix[0]*(matrix[4]*matrix[8]-matrix[5]*matrix[7])
					-matrix[1]*(matrix[3]*matrix[8]-matrix[5]*matrix[6])
					+matrix[2]*(matrix[3]*matrix[7]-matrix[4]*matrix[6]);
	}

	else
	{
		return 0;
	}

}

void MathOperationT::Inverse(int MD, double* matrix, double* inverseMatrix)
{

  double det= Determinant(MD, matrix);

  if (det<0.0)
  {
		//exception;
	}	

	if(MD==2)
	{
		double* pthis = inverseMatrix;

		*pthis++ = matrix[0]/det;	
		*pthis++ = -matrix[1]/det;
		*pthis++ = -matrix[2]/det;
		*pthis++ = matrix[3]/det;		
	}

	else if(MD==3)
	{

		double* pthis = inverseMatrix;

		*pthis++ = (matrix[4]*matrix[8] - matrix[5]*matrix[7])/det; 
		*pthis++ = (matrix[2]*matrix[7] - matrix[1]*matrix[8])/det;
		*pthis++ = (matrix[1]*matrix[5] - matrix[2]*matrix[4])/det; 

		*pthis++ = (matrix[5]*matrix[6] - matrix[3]*matrix[8])/det; 
		*pthis++ = (matrix[0]*matrix[8] - matrix[2]*matrix[6])/det; 
		*pthis++ = (matrix[2]*matrix[3] - matrix[0]*matrix[5])/det; 

		*pthis++ = (matrix[3]*matrix[7] - matrix[4]*matrix[6])/det; 
		*pthis++ = (matrix[1]*matrix[6] - matrix[0]*matrix[7])/det; 
		*pthis++ = (matrix[0]*matrix[4] - matrix[1]*matrix[3])/det; 
	}

}


void MathOperationT::PrintMatrix(int MD, const double* matrix, std::string minfo)
{
if(MD==3)
	{	
	cout<<"MATRIX\t"<<minfo<<"::"<<endl;
	cout<<matrix[0]<<'\t'<<matrix[1]<<'\t'<<matrix[2]<<endl;
	cout<<matrix[3]<<'\t'<<matrix[4]<<'\t'<<matrix[5]<<endl;
	cout<<matrix[6]<<'\t'<<matrix[7]<<'\t'<<matrix[8]<<endl;
	}	
}


void MathOperationT::PrintMatrix(int NRow, int NCol, const double* matrix, std::string minfo)
{
	cout << "MATRIX\t"<<minfo<<"::"<<endl;
	const double* temp=matrix;

	for(int i=0; i<NRow; i++)
	{
		for(int j=0;j<NCol;j++)
		{
			cout<< *temp << "\t";
			temp++;
		}

		cout << endl;
	}
}

void MathOperationT::PrintMatrix(int NRow, int NCol, const int* matrix, std::string minfo)
{
	cout << "MATRIX\t"<<minfo<<"::"<<endl;
	const int* temp=matrix;

	for(int i=0; i<NRow; i++)
	{
		for(int j=0;j<NCol;j++)
		{
			cout<< *temp << "\t";
			temp++;
		}

		cout << endl;
	}
}


void MathOperationT::PrintVector(int VD, const double* Vec, std::string vinfo)
{
	cout << "Vector \t" << vinfo << "::" << endl;

	for(int i=0;i<VD; i++)
		cout<< Vec[i] << "\t" << endl;
}

void MathOperationT::PrintVector(int VD, const int* Vec, std::string vinfo)
{
	cout << "Vector \t" << vinfo << "::" << endl;

	for(int i=0;i<VD; i++)
		cout<< Vec[i] << "\t" << endl;
}

void MathOperationT::transMatrix(int nRow, int nCol, const double* matrix, double* resultMatrix)
{
	//resultMatrix=new double[nCol*nRow];

	double* temp=resultMatrix;

	for (int i=0; i<nCol; i++)
	{
		for (int j=0; j<nRow; j++)
		{
			*temp=matrix[j*nCol+i];
			temp++;
		}
	}
}

void MathOperationT::multAB(int nRowA, int nColA, double* matrixA,
		   int nRowB, int nColB, double* matrixB, double* matrixC)
{
	double* temp=matrixC;

	for(int i=0; i<nRowA; i++)
	{
		for(int j=0; j<nColB; j++)
		{
			*temp=0;

			for(int l=0; l<nColA; l++)
			{
				*temp+=matrixA[i*nColA+l]*matrixB[l*nColB+j];
			}

			temp++;
		}
	}

}



void MathOperationT::multAB(int nRowA, int nColA, complex<double>* matrixA,
                            int nRowB, int nColB, double* matrixB, complex<double>* matrixC)
{
    complex<double>* temp=matrixC;
    
    for(int i=0; i<nRowA; i++)
    {
        for(int j=0; j<nColB; j++)
        {
            *temp=0;
            
            for(int l=0; l<nColA; l++)
            {
                *temp+=matrixA[i*nColA+l]*matrixB[l*nColB+j];
            }
            
            temp++;
        }
    }
    
}



void MathOperationT::multAB(int nRowA, int nColA, double* matrixA,
                            int nRowB, int nColB, complex<double>* matrixB, complex<double>* matrixC)
{
    complex<double>* temp=matrixC;
    
    for(int i=0; i<nRowA; i++)
    {
        for(int j=0; j<nColB; j++)
        {
            *temp=0;
            
            for(int l=0; l<nColA; l++)
            {
                *temp+=matrixA[i*nColA+l]*matrixB[l*nColB+j];
            }
            
            temp++;
        }
    }
    
}




void MathOperationT::multATB(int nRowA, int nColA, const double* matrixA,
			int nRowB, int nColB, const double* matrixB, double* matrixC)
{
	double* temp=matrixC;

	for(int i=0; i<nColA; i++)
	{
		for(int j=0; j<nColB; j++)
		{
			*temp=0;

			for(int l=0; l<nRowA; l++)
			{
				*temp+=matrixA[l*nColA+i]*matrixB[l*nColB+j];
			}

			temp++;
		}
	}

}







void MathOperationT::multABT(int MD,double* matrixA, double* matrixB, double* resultMatrix)
{

}

void MathOperationT::multATBA(int MD,double* matrixA, double* matrixB, double* resultMatrix)
{

}

void MathOperationT::multABAT(int MD,double* matrixA, double* matrixB, double* resultMatrix)
{

}



void MathOperationT::MatExtraction(int nRowA, int nColA, const double* matrixA, int subRow, int* Row, int subCol, int* Col, double* matrixB)
{
    for(int ri=0; ri<subRow; ri++)
    {
        int A_row=Row[ri];
        
        if(A_row>=nRowA)
        {
            cout<< "Attempting to extract Row " << A_row << "while Maximum Row number is " << nRowA <<endl;
            throw "In MathOperationT::MatExtraction, row is out of range\n";
        }
        
        for(int ci=0; ci<subCol; ci++)
        {
            int A_col=Col[ci];
            
            if(A_col>=nColA)
            {
                cout<< "Attempting to extract from Column " << A_col << "while Maximum Column number is " << nColA <<endl;
                throw "In MathOperationT::MatExtraction, column is out of range\n";
            }
            
            int A_position=A_row*nColA+A_col;
            int B_position=ri*subCol+ci;
            matrixB[B_position]=matrixA[A_position];
        }
    }
}


void MathOperationT::VecPlus(int VD, double a, double* VecA, double b, double* VecB, double* VecC)
{
	double* A=VecA;
	double* B=VecB;
	double* C=VecC;

	for(int i=0; i<VD; i++)
	{
		*C=a*(*A)+b*(*B);
		A++;
		B++;
		C++;
	}
}

void MathOperationT::VecPlus(int VD, double a, const double* VecA, double b, const double* VecB, double* VecC)
{
	const double* A=VecA;
	const double* B=VecB;
	double* C=VecC;

	for(int i=0; i<VD; i++)
	{
		*C=a*(*A)+b*(*B);
		A++;
		B++;
		C++;
	}
}


void MathOperationT::VecPlus(int VD, double a, double* VecA, double b, double* VecB)
{
    double* A=VecA;
    double* B=VecB;
    
    
    for(int i=0; i<VD; i++)
    {
        *A=a*(*A)+b*(*B);
        A++;
        B++;
        
    }
}



double MathOperationT::VecDot(int VD, double* VecA, double* VecB)
{
	double Dot=0;

	double* temp1=VecA;
	double* temp2=VecB;

	for (int i=0; i<VD; i++)
	{
		Dot += (*temp1)*(*temp2);
		temp1++;
		temp2++;
	}

	return Dot;
}


void MathOperationT::VecCross(double* A, double* B, double* C)
{
	C[0]=A[1]*B[2]-A[2]*B[1];
	C[1]=A[2]*B[0]-A[0]*B[2];
	C[2]=A[0]*B[1]-A[1]*B[0];

}




double MathOperationT::VecNormal(int VD, const double* VecA)
{
	double normal=0;

	const double* temp=VecA;

	for(int i=0; i<VD; i++)
	{
		normal+=pow((*temp),2);
		temp++;
	}

	normal=sqrt(normal);

	return normal;
}


void MathOperationT::MemCopy(int VD, double* A, double* B)
{
	double* temp_1=A;
	double* temp_2=B;

	for(int i=0; i<VD; i++)
	{
		*temp_2=*temp_1;
		temp_1++;
		temp_2++;
	}
}

void MathOperationT::MemCopy(int VD, int* A, int* B)
{
	int* temp_1=A;
	int* temp_2=B;

	for(int i=0; i<VD; i++)
	{
		*temp_2=*temp_1;
		temp_1++;
		temp_2++;
	}
}


void MathOperationT::VecSet(int VD, double* VecA, double a)
{
	double* temp=VecA;

	for (int i=0; i<VD; i++)
	{
		*temp=a;
		temp++;
	}
}


void MathOperationT::VecSet(int VD, int* VecA, int a)
{
	int* temp=VecA;

	for (int i=0; i<VD; i++)
	{
		*temp=a;
		temp++;
	}
}

void MathOperationT::VecScale(int VD, double* VecA, double a)
{
	for(int i=0; i<VD; i++)
		VecA[i]=VecA[i]*a;
}

