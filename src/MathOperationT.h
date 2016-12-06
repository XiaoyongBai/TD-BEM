#ifndef MATH_OPERATION_T
#define MATH_OPERATION_T

#include <iostream>
#include <string>
#include "complex"

namespace GreenFunction{

class MathOperationT
{
public:
	
	/**
 	* @function=compute the value of determinant of a given matrix
	*
	* @input:
	* 		MD=matrix dimension. Only 2 and 3 are supported
	* 		matrix=pointer to the matrix
	* @output:
	*       return determinant of the matrix*/
	static double Determinant(int MD, double* matrix);

	/**
 	* @function=compute the value of determinant of a given matrix
	*
	*  @input
	*  		MD=matrix dimension. Only 2 and 3 are supported
	* 		matrix=pointer to the matrix
	*       inverseMatrix=pointer to the inverse */
	static void Inverse(int MD, double* matrix, double* inverseMatrix);

	/**
	 * @Function=print matrix with format
	 *
	 * @Input:
	 *		MD=matrix dimension
	 *		matrix=pointer to the matrix
	 *		minfo=information of the matrix */
	static void PrintMatrix(int MD, const double* matrix, std::string minfo);

	/**
	 * @function=print matrix with format
	 *
	 * @Input:
	 *		NRow=Number of rows of the matrix
	 *		NCol=Number of columns of the matrix
	 *		matrix=pointer to the matrix
	 *		minfo=information of the matrix */
	static void PrintMatrix(int NRow, int NCol, const double* matrix, std::string minfo);
	static void PrintMatrix(int NRow, int NCol, const int* matrix, std::string minfo);


	static void PrintVector(int VD, const double* Vec, std::string minfo);
	static void PrintVector(int VD, const int* Vec, std::string vinfo );

	/**
	 * @Function=transpose the matrix
	 */
	static void transMatrix(int nRow, int nCol, const double* matrix, double* resultMatrix);

	/**
	 * @Function=matrix multiply, C=AB
	 *
	 * @Input:
	 * 		nRowA, nColA=number of row and column of matrix A
	 * 		nRowB, nColB=number of row and column of matrix B*/
	static void multAB(int nRowA, int nColA, double* matrixA,
					   int nRowB, int nColB, double* matrixB, double* matrixC);
    
    static void multAB(int nRowA, int nColA, double* matrixA,
                       int nRowB, int nColB, std::complex<double>* matrixB, std::complex<double>* matrixC);
    
    static void multAB(int nRowA, int nColA, std::complex<double>* matrixA,
                       int nRowB, int nColB, double* matrixB, std::complex<double>* matrixC);

	/**
	 * @Function=matrix multiply, C=ATB
	 *
	 * @Input:
	 * 		nRowA, nColA=number of row and column of matrix A
	 * 		nRowB, nColB=number of row and column of matrix B*/
	static void multATB(int nRowA, int nColA, const double* matrixA,
			  			int nRowB, int nColB, const double* matrixB, double* matrixC);

	//A*BT������ˣ���֧��MD=3
	static void multABT(int MD,double* matrixA, double* matrixB, double* resultMatrix);
	//AT*B*A������ˣ���֧��MD=3
	static void multATBA(int MD,double* matrixA, double* matrixB, double* resultMatrix);
	//A*B*AT������ˣ���֧��MD=3
	static void multABAT(int MD,double* matrixA, double* matrixB, double* resultMatrix);

    //Extract a sub matrix B from a matrix A
    static void MatExtraction(int nRowA, int nColA, const double* matrixA, int subRow, int* Row, int subCol, int* Col, double* matrixB);
    
    
    //Extract a sub vector B from a vector A
    template <typename Type>
    static void VecExtraction(int lengthA, const Type* vecA, int lengthB, int* position, Type* vecB);
    
	/**
	 * @Function=compute dot product of two vectors
	 *
	 * @Input
	 * 		VD=length of the vectors
	 * 		VecA, VecB=the input vectors
	 *
	 * @Output
	 * 		the dot product*/
	static double VecDot(int VD, double* VecA, double* VecB);

	/**
	 * @VecC=a*VecA+b*VecB;
	 */
	static void VecPlus(int VD, double a, double* VecA, double b, double* VecB, double* VecC);
	static void VecPlus(int VD, double a, const double* VecA, double b, const double* VecB, double* VecC);

    /**
     * @VecA+=a*VecA+b*VecB;
     */
    static void VecPlus(int VD, double a, double* VecA, double b, double* VecB);
    
	/**
	 * @Function=compute cross product of two vectors
	 *
	 * @Input
	 * 		VD=length of the vectors
	 * 		VecA, VecB=the input vectors
	 * 		VecC=result*/
	static void VecCross(double* VecA, double* VecB, double* VecC);

	/**
	 * @Function=copy the value of A to B*/
	static void MemCopy(int VD, int* A, int* B);
	static void MemCopy(int VD, double* A, double* B);

	/**
	 * @function=compute L2 normal of a vector
	 *
	 * @input
	 * 		VD=length of the vector
	 * 		VecA=the input vector
	 *
	 * @output
	 * 		normal of the vector */
	static double VecNormal(int VD, const double* VecA);


	/**
	 * @Function=set all components to a constant
	 *
	 * @Input:
	 * 		VD=length of the vector
	 * 		VecA=the vector
	 * 		a=the value we want to set*/
	static void VecSet(int VD, double* VecA, double a);
	static void VecSet(int VD, int* VecA, int a);

	/**
	 * @Function=multiply all components with the same scalar
	 *
	 * VecA=a*VecA;
	 */
	static void VecScale(int VD, double* VecA, double a);
};

    
    
    using namespace GreenFunction;
    using namespace std;
    
    template<typename Type>
    void MathOperationT::VecExtraction(int lengthA, const Type* vecA, int lengthB, int* position, Type* vecB)
    {
        int* temp=position;
        
        for(int pi=0; pi<lengthB; pi++)
        {
            int A_position=*temp;
            
            if(A_position>=lengthA)
            {
                cout<< "Attempting to extract in position " << A_position << "while length is " << lengthA <<endl;
                throw "In MathOperationT::VecExtraction, position is out of range\n";
            }
            
            vecB[pi]=vecA[A_position];
            
            temp++;
        }
    }
    
}/* name space */

#endif
