/*
 * GreenFunctionT.h
 *
 *  Created on: May 27, 2015
 *      Author: xiaoyong
 *
 *    Function: to compute time integration of Green's function, given time period
 *              and the direction vector
 */

#ifndef GREENFUNCTIONT_H_
#define GREENFUNCTIONT_H_

namespace TD_BEM{


class GreenFunctionT
{
public:
	GreenFunctionT();
	~GreenFunctionT();

	/**
	 * @Function=set material constants
	 *
	 * @input:
	 * 		E=Young's modulus
	 * 		nu=Possion's ratio
	 * 		density=density*/
	void SetMaterial(double E, double nu, double density);

	/**
	 * @Function=compute and set wave speed using fE and fMu*/
	void WaveSpeed(void);

	/* *
	 * @Function=set spatial parameters
	 *
	 * @input:
	 * 		Direction=direction vector from source to obersevation point
	 * 		Normal=outnormal of the surface */
	void SetSpatial(const double* Direction, const double* Normal);

	/**
	 * @Function=set time parameters
	 *
	 * @input:
	 * 		t=observation time
	 * 		t_1=upper bound of the time integration*/
	void SetTime(double t, double t_1);

	/**
	 * @Function=carry out the integration
	 *
	 * @input:
	 * 		type=type of the green's function
	 * 			1:Stokes' function with Heaviside time function, constant time interpolation
	 * 			2:Stokes' function with Heaviside time function, linear time interpolation
	 * 			3:Stokes' function with Step time function, constant time interpolation
	 * 			4:Stokes' function with Step time function, linear time interpolation */
	void Compute(int type);

	/**
	 * @Function=Kelvin's solution of Traction
	 */
	void KelvinTraction(void);

	/**
	 * @Function=Integrate Kelvin's solution over time for constant interpolation
	 */
	void Kelvin_Constant(void);

	/**
	 * @Function=Integrate Kelvin's solution over time for linear interpolation
	 */
	void Kelvin_Linear(void);

	/**
	 * @Function=Integrate Stokes's solution with Heaviside time function, for constant interpolation
	 */
	void Stokes_Heaviside_Constant(void);

	/**
	 * @Function=Integrate Stokes's solution with Heaviside time function, for linear interpolation
	 */
	void Stokes_Heaviside_Linear(double* DU_1, double* DU_2, double* DT_1, double* DT_2);

	/**
	 * @Function=Integrate Stokes's solution with Step time function, for constant interpolation
	 */
	void Stokes_Step_Constant(void);

	/**
	 * @Function=Integrate Stokes's solution with Step time function, for linear interpolation
	 */
	void Stokes_Step_Linear(void);

	/**
	 * @Function=Integrate Stokes's solution with ramped Heaviside time function, for linear interpolation
	 */
	void Int_Stokes_RampedHeaviside_Linear(double* DU_1, double* DU_2, double* DT_1, double* DT_2);
	void Int_Stokes_CubicBSpline_Linear(double* DU_1, double* DU_2, double* DT_1, double* DT_2);

	/**
	* @Function=compute Stokes's solution to Ramped Heaviside time function
	*/
	void Stokes_RampedHeaviside(double* U, double* T, double t);

	void Stokes_CubicBSpline(double* U, double* T, double t);

	/**
	 * @Function=Integrate multiplication of Kelvin's solution, ramped Heaviside,
	 *           and linear time interpolation
	 */
	void Int_Kelvin_RampedHeaviside_Linear(void);
	void Int_Kelvin_CubicBSpline_Linear(void);

	//convolution of ramped Heaviside and linear shape function for time
	void Int_RampedHeaviside_Linear(double& g1, double& g2);
	void Int_CubicB_Linear(double& g1, double& g2);


	/**
	 * @Function=return 1 if i==j, otherwise return 0
	 */
	int KronDelta(int i, int j);
	double Heaviside(double t);
	void CubicB1(double& B1, double& DB1, double t, double T);
	double CubicB1_int(double r, double C1, double C2, double t, double T);

	/**
	 * @Function=get pointer to the Green's functions*/
	const double* GetKelvinTraction();
	const double* GetKelvinConstant();
	const double* GetKelvinLinear();
	void GetKelvinRampedLinear(double** DT_1, double** DT_2);

	void GetStokesConstant(double** DU_1, double** DT_1);
	void GetStokesLinear(double** DU_1, double** DU_2, double** DT_1, double** DT_2);

	double Gaussian_Abs(int num_q, int q_i);
	double Gaussian_Weight(int num_q, int q_i);

private:
	/* @
	 * fE=Young's modulus
	 * fNu=Possion's ratio
	 * fMu=Shear modulus
	 * fDensity=density
	 * fCp=Compression wave speed
	 * fCs=Shear wave speed
	 * fR_v=Directional vector from the source point to observation point
	 * fR=Distance from the source point to observation point
	 * fN=Out normal of the surface on which the observation point sits
	 * fT=Time of observation
	 * fT_1=Upper time bound of the time integration
	 */
	double fE, fNu, fMu, fDensity;
	double fCp, fCs;
	double *fR_v;
	double fR;
	double *fN;
	double fT, fT_1;

	double fRampingSlope;
	double fBWidth;  //width of BSline function


	/**
	 * fDT_Kelvin_1, ..=Time integrations of Kelvin's tractions
	 * fDU_Stokes_1, ..=Time integrations of Stokes' displacements
	 * fDT_Stokes_1, ..=Time integrations of Stokes' tractions */
	double* fKelvinTraction;

	double* fDT_Kelvin_1;
	double* fDT_Kelvin_2;

	double* fDU_Stokes_1;
	double* fDU_Stokes_2;
	double* fDU_Stokes_3;
	double* fDU_Stokes_4;

	double* fDT_Stokes_1;
	double* fDT_Stokes_2;
	double* fDT_Stokes_3;
	double* fDT_Stokes_4;
};



};/* name space */

#endif /* GREENFUNCTIONT_H_ */
