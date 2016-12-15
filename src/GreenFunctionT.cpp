/*
 * GreenFunctionT.cpp
 *
 *  Created on: May 27, 2015
 *      Author: xiaoyong
 */

#include <iostream>
#include "GreenFunctionT.h"
#include "cmath"
#include "MathOperationT.h"

using namespace std;
using namespace TD_BEM;
using namespace GreenFunction;


GreenFunctionT::GreenFunctionT()
{
	fE=1;
	fNu=0;
	fMu=1;
	fDensity=1;

	fCp=1;
	fCs=0.5;

	fT=1;
	fT_1=1;

	fR=1;

	fRampingSlope=1/fT_1;

	fR_v=new double[3];
	fN=new double[3];

	/*allocation for Green's functions*/
	fKelvinTraction=new double[9];

	fDT_Kelvin_1=new double[9];
	fDT_Kelvin_2=new double[9];

	fDU_Stokes_1=new double[9];
	fDU_Stokes_2=new double[9];
	fDU_Stokes_3=new double[9];
	fDU_Stokes_4=new double[9];

	fDT_Stokes_1=new double[9];
	fDT_Stokes_2=new double[9];
	fDT_Stokes_3=new double[9];
	fDT_Stokes_4=new double[9];
}

GreenFunctionT::~GreenFunctionT()
{
	delete[] fR_v;
	delete[] fN;

	delete[] fKelvinTraction;

	delete[] fDT_Kelvin_1;
	delete[] fDT_Kelvin_2;

	delete[] fDU_Stokes_1;
	delete[] fDU_Stokes_2;
	delete[] fDU_Stokes_3;
	delete[] fDU_Stokes_4;

	delete[] fDT_Stokes_1;
	delete[] fDT_Stokes_2;
	delete[] fDT_Stokes_3;
	delete[] fDT_Stokes_4;
}

void GreenFunctionT::WaveSpeed()
{
	double L;
	fMu=fE/(2*(1+fNu));
	L=fE*fNu/((1+fNu)*(1-2*fNu));

	fCp=sqrt((L+2*fMu)/fDensity);
	fCs=sqrt(fMu/fDensity);
}

void GreenFunctionT::SetMaterial(double E, double nu, double density)
{
	fE=E;
	fNu=nu;
	fDensity=density;

	WaveSpeed();
}

void GreenFunctionT::SetSpatial(const double* Direction, const double* Normal)
{
	fR=MathOperationT::VecNormal(3, Direction);

	for(int i=0; i<3; i++)
		fR_v[i]=Direction[i]/fR;

	double R=MathOperationT::VecNormal(3, Normal);

	for(int i=0; i<3; i++)
		fN[i]=Normal[i]/R;

	//MathOperationT::PrintVector(3, fR_v, "fR_v");
	//MathOperationT::PrintVector(3, fN, "fN");
}

void GreenFunctionT::SetTime(double t, double t_1)
{
	fT=t;
	fT_1=t_1;
	fRampingSlope=2/fT_1;
	fBWidth=0.7*fT_1;
}

void GreenFunctionT::Compute(int type)
{
	switch (type){
	case 1:
		Kelvin_Constant();
		Stokes_Heaviside_Constant();
		break;
	case 2:
		Kelvin_Linear();
		Stokes_Heaviside_Linear(fDU_Stokes_1, fDU_Stokes_2, fDT_Stokes_1,fDT_Stokes_2);
		break;
	case 3:
		Kelvin_Constant();
		Stokes_Step_Constant();
		break;
	case 4:
		Kelvin_Linear();
		Stokes_Step_Linear();
		break;
	case 5:
		Int_Stokes_CubicBSpline_Linear(fDU_Stokes_1, fDU_Stokes_2, fDT_Stokes_1,fDT_Stokes_2);
		Int_Kelvin_CubicBSpline_Linear();
		break;
	default:
		cout<<"GreenFunctionT::Compute unsupported type " <<type <<endl;
	}
}

void GreenFunctionT::KelvinTraction(void)
{
	double* temp=fKelvinTraction;

	//Compute Kelvin's traction

	double rdn=MathOperationT::VecDot(3, fR_v, fN);

	double divider=-(8*M_PI)*(1-fNu)*fR*fR;
	for(int i=0; i<3; i++)
	{
		for(int k=0; k<3; k++)
		{
			int dik=KronDelta(i,k);
			double rirk=fR_v[i]*fR_v[k];
			double rink=fR_v[i]*fN[k];
			double rkni=fR_v[k]*fN[i];

			*temp = ((1-2*fNu)*dik+3*rirk)*rdn-(1-2*fNu)*(rkni-rink);
			*temp /= divider;

			temp++;
		}
	}

	//MathOperationT::PrintMatrix(3, fKelvinTraction, "fKelvinTraction");
}



void GreenFunctionT::Kelvin_Constant(void)
{
	double* temp=fDT_Kelvin_1;

	double rdn=MathOperationT::VecDot(3, fR_v, fN);

	//Integrate traction over time
	double divider=-(8*M_PI)*(1-fNu)*fR*fR;
	for(int i=0; i<3; i++)
	{
		for(int k=0; k<3; k++)
		{
			int dik=KronDelta(i,k);
			double rirk=fR_v[i]*fR_v[k];
			double rink=fR_v[i]*fN[k];
			double rkni=fR_v[k]*fN[i];

			*temp = ((1-2*fNu)*dik+3*rirk)*rdn-(1-2*fNu)*(rkni-rink);
			*temp /= divider;
			*temp *= fT_1;

			temp++;
		}
	}

	//MathOperationT::PrintMatrix(3, fDT_Kelvin_1, "fDT_Kelvin");
}

void GreenFunctionT::Kelvin_Linear(void)
{
	double* temp=fDT_Kelvin_1;

	double rdn=MathOperationT::VecDot(3, fR_v, fN);

	//Integrate traction over time
	double divider=-(8*M_PI)*(1-fNu)*fR*fR;
	for(int i=0; i<3; i++)
	{
		for(int k=0; k<3; k++)
		{
			int dik=KronDelta(i,k);
			double rirk=fR_v[i]*fR_v[k];
			double rink=fR_v[i]*fN[k];
			double rkni=fR_v[k]*fN[i];

			*temp = ((1-2*fNu)*dik+3*rirk)*rdn-(1-2*fNu)*(rkni-rink);
			*temp /= divider;
			*temp *= 0.5*fT_1;

			temp++;
		}
	}

	//MathOperationT::PrintMatrix(3, fDT_Kelvin_1, "fDT_Kelvin");
}


void GreenFunctionT::Stokes_Heaviside_Constant(void)
{
	//waitting for implementation
}



void GreenFunctionT::Stokes_Heaviside_Linear(double* DU_1, double* DU_2, double* DT_1, double* DT_2)
{
	/*************************************
	 *  integral of time functions
	 *************************************/
	double I1, I2, I3, I4, I5, I6;
	double J1, J2, J3, J4, J5, J6;

	double b;

	double t_ref=fT-fR/fCp;
	if (fT_1<=t_ref)
	{
		b=fT_1;
		I5=0;
		J5=0;
	}
	else if(0<=t_ref && fT_1>t_ref)
	{
		b=t_ref;
		I5=(fT_1-t_ref)/fT_1;
		J5=t_ref/fT_1;
	}
	else
	{
		b=0;
		I5=0;
		J5=0;
	}

	I1=( -3*pow(b,4) + (8*fT+4*fT_1)*pow(b,3) -(6*pow(fT,2)+12*fT*fT_1)*pow(b,2)+ 12*fT_1*pow(fT,2)*b )/(12*pow(fR, 2))
			-(fT_1*b-pow(b,2)/2.0)/pow(fCp,2);
	I1=I1/fT_1;

	J1=-( -3*pow(b,4) + 8*fT*pow(b,3) - 6*pow(fT,2)*pow(b,2))/(12*pow(fR, 2))
					-(pow(b,2)/2.0)/pow(fCp,2);
	J1=J1/fT_1;

	I3=(fT_1*b-pow(b,2)/2.0)/fT_1;
	J3=pow(b,2)/(2.0*fT_1);


	/*------------------------------------------------
	 -------------------------------------------------*/
	t_ref=fT-fR/fCs;
	if (fT_1<=t_ref)
	{
		b=fT_1;
		I6=0;
		J6=0;
	}
	else if(0<=t_ref && fT_1>t_ref)
	{
		b=t_ref;
		I6=(fT_1-t_ref)/fT_1;
		J6=t_ref/fT_1;
	}
	else
	{
		b=0;
		I6=0;
		J6=0;
	}

	I2=( -3*pow(b,4) + (8*fT+4*fT_1)*pow(b,3) -(6*pow(fT,2)+12*fT*fT_1)*pow(b,2)+ 12*fT_1*pow(fT,2)*b )/(12*pow(fR, 2))
			-(fT_1*b-pow(b,2)/2.0)/pow(fCs,2);
	I2=I2/fT_1;

	J2=-( -3*pow(b,4) + 8*fT*pow(b,3) - 6*pow(fT,2)*pow(b,2))/(12*pow(fR, 2))
					-(pow(b,2)/2.0)/pow(fCs,2);
	J2=J2/fT_1;

	I4=(fT_1*b-pow(b,2)/2.0)/fT_1;
	J4=pow(b,2)/(2.0*fT_1);

	/*------------------------------------------------
	 -------------------------------------------------*/
	double It_1, It_2, It_3, It_4, It_5;
	double Jt_1, Jt_2, Jt_3, Jt_4, Jt_5;

	It_1=0.5*(I1-I2);
	It_2=1/pow(fCp,2)*I3-1/pow(fCs,2)*I4;
	It_3=I4;

	Jt_1=0.5*(J1-J2);
	Jt_2=1.0/pow(fCp,2)*J3-1.0/pow(fCs,2)*J4;
	Jt_3=J4;

	double* temp_1=DU_1;
	double* temp_2=DU_2;

	double divider=4*M_PI*fDensity*fR;
	for (int i=0; i<3; i++)
	{
		for (int k=0; k<3; k++)
		{
			int dik=KronDelta(i,k);

			double rirk=fR_v[i]*fR_v[k];
			double s1=3*rirk-dik;
			double s2=rirk;
			double s3=dik/pow(fCs,2);

			*temp_1=s1*It_1+s2*It_2+s3*It_3;
			*temp_1/=divider;

			*temp_2=s1*Jt_1+s2*Jt_2+s3*Jt_3;
			*temp_2/=divider;

			temp_1++;
			temp_2++;
		}
	}

	//MathOperationT::PrintMatrix(3, fDU_Stokes_1, "fDU_Stokes_1");
	//MathOperationT::PrintMatrix(3, fDU_Stokes_2, "fDU_Stokes_2");

	/*------------------------------------------------
	 -------------------------------------------------*/
	It_1=0.5*(I1-I2);
	It_2=I4-pow(fCs/fCp,2)*I3;
	It_3=I6-pow(fCs/fCp,3)*I5;
	It_4=I3+fR/fCp*I5;
	It_5=I4+fR/fCs*I6;

	Jt_1=0.5*(J1-J2);
	Jt_2=J4-pow(fCs/fCp,2)*J3;
	Jt_3=J6-pow(fCs/fCp,3)*J5;
	Jt_4=J3+fR/fCp*J5;
	Jt_5=J4+fR/fCs*J6;

	/*------------------------------------------------
	 -------------------------------------------------*/
	double rdn=MathOperationT::VecDot(3, fR_v, fN);

	temp_1=DT_1;
	temp_2=DT_2;

	for (int i=0; i<3; i++)
	{
		for (int k=0; k<3; k++)
		{
			int dik=KronDelta(i,k);

			double RR=pow(fR,2);

			double rirk=fR_v[i]*fR_v[k];
			double nirk=fN[i]*fR_v[k];
			double nkri=fN[k]*fR_v[i];

			double s1=-6*pow(fCs,2)*( 5*rirk*rdn - (nirk+nkri+dik*rdn) )/RR;
			double s2=2*(6*rirk*rdn-(nirk+nkri+dik*rdn) )/RR;
			double s3=2*rirk*rdn/(fCs*fR);
			double s4=-nirk*(1-2*pow(fCs/fCp,2))/RR;
			double s5=-(dik*rdn+nkri)/RR;

			*temp_1 = s1*It_1+s2*It_2+s3*It_3+s4*It_4+s5*It_5;
			*temp_1 /= 4.0*M_PI;

			*temp_2 = s1*Jt_1+s2*Jt_2+s3*Jt_3+s4*Jt_4+s5*Jt_5;
			*temp_2 /= 4.0*M_PI;

			temp_1++;
			temp_2++;
		}
	}

	/*MathOperationT::PrintVector(3, fR_v, "direction");
	MathOperationT::PrintVector(3, fN, "normal");
	cout<< "r " << fR << endl;*/
	//cout << "I" << I1 <<" " << I2 <<" "<< I3 <<" "<< I4 <<" "<< I5 <<" " << I6 <<endl;
	//cout << "J" << J1 <<" " << J2 <<" "<< J3 <<" "<< J4 <<" "<< J5 <<" " << J6 <<endl;

	//MathOperationT::PrintMatrix(3, fDT_Stokes_1, "fDT_Stokes_1");
	//MathOperationT::PrintMatrix(3, fDT_Stokes_2, "fDT_Stokes_2");
}


void GreenFunctionT::Stokes_Step_Constant(void)
{
	//waitting for implementation
}

void GreenFunctionT::Stokes_Step_Linear(void)
{
	Stokes_Heaviside_Linear(fDU_Stokes_1, fDU_Stokes_2, fDT_Stokes_1,fDT_Stokes_2);

	double* DU_1_temp=new double[9];
	double* DU_2_temp=new double[9];
	double* DT_1_temp=new double[9];
	double* DT_2_temp=new double[9];

	fT=fT-fT_1;
	Stokes_Heaviside_Linear(DU_1_temp, DU_2_temp,DT_1_temp, DT_2_temp);
	fT=fT+fT_1;

	MathOperationT::VecPlus(9, 1, fDU_Stokes_1, -1, DU_1_temp, fDU_Stokes_1);
	MathOperationT::VecPlus(9, 1, fDU_Stokes_2, -1, DU_2_temp, fDU_Stokes_2);
	MathOperationT::VecPlus(9, 1, fDT_Stokes_1, -1, DT_1_temp, fDT_Stokes_1);
	MathOperationT::VecPlus(9, 1, fDT_Stokes_2, -1, DT_2_temp, fDT_Stokes_2);

	delete[] DU_1_temp;
	delete[] DU_2_temp;
	delete[] DT_1_temp;
	delete[] DT_2_temp;
}


void GreenFunctionT::Stokes_RampedHeaviside(double* U, double* T, double t)
{
	double HC1, HC2, HKC1, HKC2;
	double GC1, GC2, DGC1, DGC2;
	double I1, Dis_I2, Dis_I3;
	double T_I2, T_I3, T_I4, T_I5;

	double k=fRampingSlope;

	HC1=Heaviside(t-fR/fCp);
	HC2=Heaviside(t-fR/fCs);
	HKC1=Heaviside(t-1/k-fR/fCp);
	HKC2=Heaviside(t-1/k-fR/fCs);

	GC1=k*(t-fR/fCp)*(HC1-HKC1)+HKC1;
	GC2=k*(t-fR/fCs)*(HC2-HKC2)+HKC2;
	DGC1=k*(HC1-HKC1);
	DGC2=k*(HC2-HKC2);

	I1=   k*( (pow(t/fR,2)-1/pow(fCp,2))*t/2-(pow(t/fR,3)-1/pow(fCp,3))*fR/3 )*HC1;
	I1=I1-k*( (pow(t/fR,2)-1/pow(fCs,2))*t/2-(pow(t/fR,3)-1/pow(fCs,3))*fR/3 )*HC2;
	I1=I1+((1-k*t)*(pow(t-1/k,2)/pow(fR,2)-1/pow(fCp,2))/2+k*fR*(pow(t-1/k,3)/pow(fR,3)-1/pow(fCp,3))/3)*HKC1;
	I1=I1-((1-k*t)*(pow(t-1/k,2)/pow(fR,2)-1/pow(fCs,2))/2+k*fR*(pow(t-1/k,3)/pow(fR,3)-1/pow(fCs,3))/3)*HKC2;

	Dis_I2=GC1/pow(fCp,2)-GC2/pow(fCs,2);
	Dis_I3=GC2;

	T_I2=GC2-pow(fCs/fCp,2)*GC1;
	T_I3=DGC2-pow(fCs/fCp,3)*DGC1;
	T_I4=GC1+(fR/fCp)*DGC1;
	T_I5=GC2+(fR/fCs)*DGC2;


	double rdn=MathOperationT::VecDot(3, fR_v, fN);


	double divider=4*M_PI*fDensity*fR;

	double* U_temp=U;
	double* T_temp=T;

	for (int i=0; i<3; i++)
	{
		for (int k=0; k<3; k++)
		{
			int dik=KronDelta(i,k);

			double RR=pow(fR,2);

			double rirk=fR_v[i]*fR_v[k];
			double nirk=fN[i]*fR_v[k];
			double nkri=fN[k]*fR_v[i];

			*U_temp = (3*rirk-dik)*I1 + rirk*Dis_I2 + dik/pow(fCs,2)*Dis_I3;
			*U_temp/= divider;

			double s1=-6*pow(fCs,2)*( 5*rirk*rdn - (nirk+nkri+dik*rdn) )/RR;
			double s2=2*(6*rirk*rdn-(nirk+nkri+dik*rdn) )/RR;
			double s3=2*rirk*rdn/(fCs*fR);
			double s4=-nirk*(1-2*pow(fCs/fCp,2))/RR;
			double s5=-(dik*rdn+nkri)/RR;

			*T_temp = s1*I1+s2*T_I2+s3*T_I3+s4*T_I4+s5*T_I5;
			*T_temp /=4.0*M_PI;

			U_temp++;
			T_temp++;
		}
	}

}


void GreenFunctionT::Stokes_CubicBSpline(double* U, double* T, double t)
{
	double GC1, GC2, DGC1, DGC2;
	double I1, Dis_I2, Dis_I3;
	double T_I2, T_I3, T_I4, T_I5;

	CubicB1(GC1, DGC1, t-fR/fCp, fBWidth);
	CubicB1(GC2, DGC2, t-fR/fCs, fBWidth);

	I1=CubicB1_int(fR, fCp, fCs, t, fBWidth);
	Dis_I2=GC1/pow(fCp,2)-GC2/pow(fCs,2);
	Dis_I3=GC2;

	T_I2=GC2-pow(fCs/fCp,2)*GC1;
	T_I3=DGC2-pow(fCs/fCp,3)*DGC1;
	T_I4=GC1+(fR/fCp)*DGC1;
	T_I5=GC2+(fR/fCs)*DGC2;


	double rdn=MathOperationT::VecDot(3, fR_v, fN);


	double divider=4*M_PI*fDensity*fR;

	double* U_temp=U;
	double* T_temp=T;

	for (int i=0; i<3; i++)
	{
		for (int k=0; k<3; k++)
		{
			int dik=KronDelta(i,k);

			double RR=pow(fR,2);

			double rirk=fR_v[i]*fR_v[k];
			double nirk=fN[i]*fR_v[k];
			double nkri=fN[k]*fR_v[i];

			*U_temp = (3*rirk-dik)*I1 + rirk*Dis_I2 + dik/pow(fCs,2)*Dis_I3;
			*U_temp/= divider;

			double s1=-6*pow(fCs,2)*( 5*rirk*rdn - (nirk+nkri+dik*rdn) )/RR;
			double s2=2*(6*rirk*rdn-(nirk+nkri+dik*rdn) )/RR;
			double s3=2*rirk*rdn/(fCs*fR);
			double s4=-nirk*(1-2*pow(fCs/fCp,2))/RR;
			double s5=-(dik*rdn+nkri)/RR;

			*T_temp = s1*I1+s2*T_I2+s3*T_I3+s4*T_I4+s5*T_I5;
			*T_temp /=4.0*M_PI;

			U_temp++;
			T_temp++;
		}
	}


	//MathOperationT::PrintVector(3, fR_v, "direction");
	//MathOperationT::PrintVector(3, fN, "normal");
	//cout<< "r " << fR << endl;
	//cout << "t is " << t << endl;
	//MathOperationT::PrintMatrix(3, U, "U");
	//MathOperationT::PrintMatrix(3, T, "T");
}



void GreenFunctionT::Int_Stokes_RampedHeaviside_Linear(double* DU_1, double* DU_2, double* DT_1, double* DT_2)
{
	MathOperationT::VecSet(9, DU_1, 0);
	MathOperationT::VecSet(9, DU_2, 0);
	MathOperationT::VecSet(9, DT_1, 0);
	MathOperationT::VecSet(9, DT_2, 0);

	int sample_size=40;

	double* U_temp=new double[9];
	double* T_temp=new double[9];

	for(int i=0; i<sample_size; i++)
	{
		double abs= Gaussian_Abs(sample_size, i);
		double w= Gaussian_Weight(sample_size, i);

		double phi_1=0.5*(1-abs);
		double phi_2=0.5*(1+abs);
		double J=fT_1*0.5;

		double t_temp=phi_2*fT_1;

		Stokes_RampedHeaviside(U_temp, T_temp, fT-t_temp);

		MathOperationT::VecPlus(9, 1, DU_1, w*J*phi_1, U_temp);
		MathOperationT::VecPlus(9, 1, DU_2, w*J*phi_2, U_temp);
		MathOperationT::VecPlus(9, 1, DT_1, w*J*phi_1, T_temp);
		MathOperationT::VecPlus(9, 1, DT_2, w*J*phi_2, T_temp);
	}

	delete[] U_temp;
	delete[] T_temp;
}


void GreenFunctionT::Int_Stokes_CubicBSpline_Linear(double* DU_1, double* DU_2, double* DT_1, double* DT_2)
{
	MathOperationT::VecSet(9, DU_1, 0);
	MathOperationT::VecSet(9, DU_2, 0);
	MathOperationT::VecSet(9, DT_1, 0);
	MathOperationT::VecSet(9, DT_2, 0);

	int sample_size=20;

	double* U_temp=new double[9];
	double* T_temp=new double[9];

	for(int i=0; i<sample_size; i++)
	{
		double abs= Gaussian_Abs(sample_size, i);
		double w= Gaussian_Weight(sample_size, i);

		double phi_1=0.5*(1-abs);
		double phi_2=0.5*(1+abs);
		double J=fT_1*0.5;

		double t_temp=phi_2*fT_1;

		Stokes_CubicBSpline(U_temp, T_temp, fT-t_temp);

		MathOperationT::VecPlus(9, 1, DU_1, w*J*phi_1, U_temp);
		MathOperationT::VecPlus(9, 1, DU_2, w*J*phi_2, U_temp);
		MathOperationT::VecPlus(9, 1, DT_1, w*J*phi_1, T_temp);
		MathOperationT::VecPlus(9, 1, DT_2, w*J*phi_2, T_temp);
	}

	delete[] U_temp;
	delete[] T_temp;
}


void GreenFunctionT::Int_Kelvin_RampedHeaviside_Linear(void)
{
	double g1=0;
	double g2=0;

	Int_RampedHeaviside_Linear(g1, g2);

	KelvinTraction();

	for(int i=0; i<9; i++)
	{
		fDT_Kelvin_1[i]=fKelvinTraction[i]*g1;
		fDT_Kelvin_2[i]=fKelvinTraction[i]*g2;
	}

	//MathOperationT::PrintMatrix(3, fDT_Kelvin_1, "fKelvinTraction_1");
	//MathOperationT::PrintMatrix(3, fDT_Kelvin_2, "fKelvinTraction_2");


}

void GreenFunctionT::Int_Kelvin_CubicBSpline_Linear(void)
{
	double g1=0;
	double g2=0;

	Int_CubicB_Linear(g1, g2);

	KelvinTraction();

	for(int i=0; i<9; i++)
	{
		fDT_Kelvin_1[i]=fKelvinTraction[i]*g1;
		fDT_Kelvin_2[i]=fKelvinTraction[i]*g2;
	}

	//MathOperationT::PrintMatrix(3, fDT_Kelvin_1, "fKelvinTraction_1");
	//MathOperationT::PrintMatrix(3, fDT_Kelvin_2, "fKelvinTraction_2");


}



void GreenFunctionT::Int_RampedHeaviside_Linear(double& g1, double& g2)
{
	g1=0;
	g2=0;

	int sample_size=40;


	for(int i=0; i<sample_size; i++)
	{
		double abs= Gaussian_Abs(sample_size, i);
		double w= Gaussian_Weight(sample_size, i);

		double phi_1=0.5*(1-abs);
		double phi_2=0.5*(1+abs);
		double J=fT_1*0.5;

		double t_temp=phi_2*fT_1;

		double t=fT-t_temp;
		double g=fRampingSlope*t*(Heaviside(t)-Heaviside(t-1/fRampingSlope))+Heaviside(t-1/fRampingSlope);

		g1+=phi_1*g*w*J;
		g2+=phi_2*g*w*J;
	}

}

void GreenFunctionT::Int_CubicB_Linear(double& g1, double& g2)
{
	g1=0;
	g2=0;

	int sample_size=20;

	for(int i=0; i<sample_size; i++)
	{
		double abs= Gaussian_Abs(sample_size, i);
		double w= Gaussian_Weight(sample_size, i);

		double phi_1=0.5*(1-abs);
		double phi_2=0.5*(1+abs);
		double J=fT_1*0.5;

		double t_temp=phi_2*fT_1;

		double t=fT-t_temp;
		double g, dg;
		CubicB1(g,dg,t,fBWidth);

		g1+=phi_1*g*w*J;
		g2+=phi_2*g*w*J;
	}

}



int GreenFunctionT::KronDelta(int i, int j)
{
	if (i==j)
		return 1;
	else
		return 0;
}

double GreenFunctionT::Heaviside(double t)
{
	if (t>=0)
		return 1.0;
	else
		return 0.0;
}


void GreenFunctionT::CubicB1(double& B1, double& DB1, double t, double T)
{
	t=t/T;
	double tp2=pow(t,2);
	double tp3=pow(t,3);

	double p1, p2, p3, p4;

	p1=32.0/3.0;
	B1=(Heaviside(t)-Heaviside(t-0.25))*p1*tp3;

	p1=-32;
	p2=32;
	p3=-8;
	p4=2.0/3.0;
	B1+=(Heaviside(t-0.25)-Heaviside(t-0.5))*(p1*tp3+p2*tp2+p3*t+p4);

	p1=32;
	p2=-64;
	p3=40;
	p4=-22.0/3.0;
	B1+=(Heaviside(t-0.5)-Heaviside(t-0.75))*(p1*tp3+p2*tp2+p3*t+p4);

	p1=-32.0/3.0;
	p2=32;
	p3=-32;
	p4=32.0/3.0;
	B1+=(Heaviside(t-0.75)-Heaviside(t-1))*(p1*tp3+p2*tp2+p3*t+p4);

	DB1 =(Heaviside(t)-Heaviside(t-0.25))*32*tp2;
	DB1+=(Heaviside(t-0.25)-Heaviside(t-0.5))*(-96*tp2+64*t-8);
	DB1+=(Heaviside(t-0.5)-Heaviside(t-0.75))*(96*tp2-128*t+40);
	DB1+=(Heaviside(t-0.75)-Heaviside(t-1))*(-32*tp2+64*t-32);
	DB1=DB1/T;

}

double GreenFunctionT::CubicB1_int(double r, double C1, double C2, double t, double T)
{
	double B=0;

	//compute common factors
	double t1=r/C1;
	double t2=r/C2;

	double tp5=pow(t,5);
	double tp4=pow(t,4);
	double tp3=pow(t,3);
	double tp2=pow(t,2);

	double Tp5=pow(T,5);
	double Tp4=pow(T,4);
	double Tp3=pow(T,3);
	double Tp2=pow(T,2);

	double t1p5=pow(t1,5);
	double t1p4=pow(t1,4);
	double t1p3=pow(t1,3);
	double t1p2=pow(t1,2);

	double t2p5=pow(t2,5);
	double t2p4=pow(t2,4);
	double t2p3=pow(t2,3);
	double t2p2=pow(t2,2);

	double J1_cubic=4*t1p5-15*t1p4*t+20*t1p3*tp2-10*t1p2*tp3;
	double J2_cubic=4*t2p5-15*t2p4*t+20*t2p3*tp2-10*t2p2*tp3;

	double J1_quad=-3*t1p4+8*t1p3*t-6*t1p2*tp2;
	double J2_quad=-3*t2p4+8*t2p3*t-6*t2p2*tp2;

	double J1_lin=2*t1p3-3*t1p2*t;
	double J2_lin=2*t2p3-3*t2p2*t;

	double J1_const=-t1p2;
	double J2_const=-t2p2;

	//Integral of 1st segment of B-Spline function
	double T1=0;
	double T2=0.25*T;

	double K1_cubic, K2_cubic, K1_quad, K2_quad, K1_lin, K2_lin, K1_const, K2_const;
	K1_cubic=tp5-5*t*pow(T1,4)+4*pow(T1,5);
	K2_cubic=tp5-5*t*pow(T2,4)+4*pow(T2,5);

	double p_cubic, p_quad, p_lin, p_const;
	p_cubic=32.0/(3.0*Tp3);

	double I_cubic, I_quad, I_lin, I_const;

	I_cubic = Heaviside(t-t1-T1)*(J1_cubic+K1_cubic)-Heaviside(t-t2-T1)*(J2_cubic+K1_cubic);
	I_cubic-=(Heaviside(t-t1-T2)*(J1_cubic+K2_cubic)-Heaviside(t-t2-T2)*(J2_cubic+K2_cubic));
	I_cubic=I_cubic/(20.0*pow(r,2));

	B+=p_cubic*I_cubic;

	//*************************************************
	//Integral of 2nd segment of B-Spline function
	T1=0.25*T;
	T2=0.5*T;

	K1_cubic=tp5-5*t*pow(T1,4)+4*pow(T1,5);
	K2_cubic=tp5-5*t*pow(T2,4)+4*pow(T2,5);

	K1_quad=tp4-4*t*pow(T1,3)+3*pow(T1,4);
	K2_quad=tp4-4*t*pow(T2,3)+3*pow(T2,4);

	K1_lin=pow(t-T1, 2)*(t+2*T1);
	K2_lin=pow(t-T2, 2)*(t+2*T2);

	K1_const=pow(t-T1,2);
	K2_const=pow(t-T2,2);

	I_cubic = Heaviside(t-t1-T1)*(J1_cubic+K1_cubic)-Heaviside(t-t2-T1)*(J2_cubic+K1_cubic);
	I_cubic-=(Heaviside(t-t1-T2)*(J1_cubic+K2_cubic)-Heaviside(t-t2-T2)*(J2_cubic+K2_cubic));
	I_cubic=I_cubic/(20.0*pow(r,2));

	I_quad = Heaviside(t-t1-T1)*(J1_quad+K1_quad)-Heaviside(t-t2-T1)*(J2_quad+K1_quad);
	I_quad-=(Heaviside(t-t1-T2)*(J1_quad+K2_quad)-Heaviside(t-t2-T2)*(J2_quad+K2_quad));
	I_quad =I_quad/(12.0*pow(r,2));

	I_lin = Heaviside(t-t1-T1)*(J1_lin+K1_lin)-Heaviside(t-t2-T1)*(J2_lin+K1_lin);
	I_lin-=(Heaviside(t-t1-T2)*(J1_lin+K2_lin)-Heaviside(t-t2-T2)*(J2_lin+K2_lin));
	I_lin=I_lin/(6.0*pow(r,2));

	I_const = Heaviside(t-t1-T1)*(J1_const+K1_const)-Heaviside(t-t2-T1)*(J2_const+K1_const);
	I_const-=(Heaviside(t-t1-T2)*(J1_const+K2_const)-Heaviside(t-t2-T2)*(J2_const+K2_const));
	I_const=I_const/(2.0*pow(r,2));

	p_cubic=-32.0/Tp3;
	p_quad=32.0/Tp2;
	p_lin=-8.0/T;
	p_const=2.0/3.0;

	B+=p_cubic*I_cubic+p_quad*I_quad+p_lin*I_lin+p_const*I_const;

	//*************************************************
	//Integral of 3rd segment of B-Spline function
	T1=0.5*T;
	T2=0.75*T;

	K1_cubic=tp5-5*t*pow(T1,4)+4*pow(T1,5);
	K2_cubic=tp5-5*t*pow(T2,4)+4*pow(T2,5);

	K1_quad=tp4-4*t*pow(T1,3)+3*pow(T1,4);
	K2_quad=tp4-4*t*pow(T2,3)+3*pow(T2,4);

	K1_lin=pow(t-T1, 2)*(t+2*T1);
	K2_lin=pow(t-T2, 2)*(t+2*T2);

	K1_const=pow(t-T1,2);
	K2_const=pow(t-T2,2);

	I_cubic = Heaviside(t-t1-T1)*(J1_cubic+K1_cubic)-Heaviside(t-t2-T1)*(J2_cubic+K1_cubic);
	I_cubic-=(Heaviside(t-t1-T2)*(J1_cubic+K2_cubic)-Heaviside(t-t2-T2)*(J2_cubic+K2_cubic));
	I_cubic=I_cubic/(20.0*pow(r,2));

	I_quad = Heaviside(t-t1-T1)*(J1_quad+K1_quad)-Heaviside(t-t2-T1)*(J2_quad+K1_quad);
	I_quad-=(Heaviside(t-t1-T2)*(J1_quad+K2_quad)-Heaviside(t-t2-T2)*(J2_quad+K2_quad));
	I_quad =I_quad/(12.0*pow(r,2));

	I_lin = Heaviside(t-t1-T1)*(J1_lin+K1_lin)-Heaviside(t-t2-T1)*(J2_lin+K1_lin);
	I_lin-=(Heaviside(t-t1-T2)*(J1_lin+K2_lin)-Heaviside(t-t2-T2)*(J2_lin+K2_lin));
	I_lin=I_lin/(6.0*pow(r,2));

	I_const = Heaviside(t-t1-T1)*(J1_const+K1_const)-Heaviside(t-t2-T1)*(J2_const+K1_const);
	I_const-=(Heaviside(t-t1-T2)*(J1_const+K2_const)-Heaviside(t-t2-T2)*(J2_const+K2_const));
	I_const=I_const/(2.0*pow(r,2));

	p_cubic=32.0/Tp3;
	p_quad=-64.0/Tp2;
	p_lin=40/T;
	p_const=-22.0/3.0;

	B+=p_cubic*I_cubic+p_quad*I_quad+p_lin*I_lin+p_const*I_const;


	//*************************************************
	//Integral of 4th segment of B-Spline function
	T1=0.75*T;
	T2=T;

	K1_cubic=tp5-5*t*pow(T1,4)+4*pow(T1,5);
	K2_cubic=tp5-5*t*pow(T2,4)+4*pow(T2,5);

	K1_quad=tp4-4*t*pow(T1,3)+3*pow(T1,4);
	K2_quad=tp4-4*t*pow(T2,3)+3*pow(T2,4);

	K1_lin=pow(t-T1, 2)*(t+2*T1);
	K2_lin=pow(t-T2, 2)*(t+2*T2);

	K1_const=pow(t-T1,2);
	K2_const=pow(t-T2,2);

	I_cubic = Heaviside(t-t1-T1)*(J1_cubic+K1_cubic)-Heaviside(t-t2-T1)*(J2_cubic+K1_cubic);
	I_cubic-=(Heaviside(t-t1-T2)*(J1_cubic+K2_cubic)-Heaviside(t-t2-T2)*(J2_cubic+K2_cubic));
	I_cubic=I_cubic/(20.0*pow(r,2));

	I_quad = Heaviside(t-t1-T1)*(J1_quad+K1_quad)-Heaviside(t-t2-T1)*(J2_quad+K1_quad);
	I_quad-=(Heaviside(t-t1-T2)*(J1_quad+K2_quad)-Heaviside(t-t2-T2)*(J2_quad+K2_quad));
	I_quad =I_quad/(12.0*pow(r,2));

	I_lin = Heaviside(t-t1-T1)*(J1_lin+K1_lin)-Heaviside(t-t2-T1)*(J2_lin+K1_lin);
	I_lin-=(Heaviside(t-t1-T2)*(J1_lin+K2_lin)-Heaviside(t-t2-T2)*(J2_lin+K2_lin));
	I_lin=I_lin/(6.0*pow(r,2));

	I_const = Heaviside(t-t1-T1)*(J1_const+K1_const)-Heaviside(t-t2-T1)*(J2_const+K1_const);
	I_const-=(Heaviside(t-t1-T2)*(J1_const+K2_const)-Heaviside(t-t2-T2)*(J2_const+K2_const));
	I_const=I_const/(2.0*pow(r,2));

	p_cubic=-32.0/(3.0*Tp3);
	p_quad=32/Tp2;
	p_lin=-32.0/T;
	p_const=32.0/3.0;

	B+=p_cubic*I_cubic+p_quad*I_quad+p_lin*I_lin+p_const*I_const;

	return B;
}



const double* GreenFunctionT::GetKelvinTraction()
{
	return fKelvinTraction;
}


const double* GreenFunctionT::GetKelvinConstant()
{
	return fDT_Kelvin_1;
}

const double* GreenFunctionT::GetKelvinLinear()
{
	return fDT_Kelvin_1;
}

void GreenFunctionT::GetKelvinRampedLinear(double** DT_1, double** DT_2)
{
	*DT_1=fDT_Kelvin_1;
	*DT_2=fDT_Kelvin_2;
}

void GreenFunctionT::GetStokesConstant(double** DU_1, double** DT_1)
{
	*DU_1=fDU_Stokes_1;

	*DT_1=fDT_Stokes_1;
}

void GreenFunctionT::GetStokesLinear(double** DU_1, double** DU_2, double** DT_1, double** DT_2)
{
	*DU_1=fDU_Stokes_1;
	*DU_2=fDU_Stokes_2;

	*DT_1=fDT_Stokes_1;
	*DT_2=fDT_Stokes_2;
}



double GreenFunctionT::Gaussian_Abs(int num_q, int q_i)
{
	switch (num_q){

		case 6:
		{
		      double shape_coords[6]={
		         		0.6612093864662645, -0.6612093864662645, -0.2386191860831969,
		                0.2386191860831969, -0.9324695142031521,  0.9324695142031521};
		      return shape_coords[q_i];

		}
		      break;


	    case 10:
	    {
	           double shape_coords[10]={
	        		    -0.148874339, 0.148874339,-0.433395394, 0.433395394,-0.679409568,
	                    0.679409568,-0.865063367, 0.865063367,-0.973906529, 0.973906529};

	           return shape_coords[q_i];
	    }

	           break;

	    case 20:
	    {
	    		double shape_coords[20]={
	    				-0.076526521, 0.076526521,-0.227785851, 0.227785851,
	    				-0.373706089, 0.373706089,-0.510867002, 0.510867002,
	    				-0.636053681, 0.636053681,-0.746331906, 0.746331906,
	    				-0.839116972, 0.839116972,-0.912234428, 0.912234428,
	    				-0.963971927, 0.963971927,-0.993128599, 0.993128599};

	    		return shape_coords[q_i];
	    }
	    		break;


	    case 40:
	    {
	    	double shape_coords[40]={
	    			-0.0387724175060508, 0.0387724175060508, -0.1160840706752552, 0.1160840706752552,
	    			-0.1926975807013711, 0.1926975807013711, -0.2681521850072537, 0.2681521850072537,
	    			-0.3419940908257585, 0.3419940908257585, -0.4137792043716050, 0.4137792043716050,
	    			-0.4830758016861787, 0.4830758016861787, -0.5494671250951282, 0.5494671250951282,
	    			-0.6125538896679802, 0.6125538896679802, -0.6719566846141796, 0.6719566846141796,
	    			-0.7273182551899271, 0.7273182551899271, -0.7783056514265194, 0.7783056514265194,
	    			-0.8246122308333117, 0.8246122308333117, -0.8659595032122595, 0.8659595032122595,
	    			-0.9020988069688743, 0.9020988069688743, -0.9328128082786765, 0.9328128082786765,
	    			-0.9579168192137917, 0.9579168192137917, -0.9772599499837743, 0.9772599499837743,
	    			-0.9907262386994570, 0.9907262386994570, -0.9982377097105593, 0.9982377097105593};

	    		return shape_coords[q_i];
	    }

	    		break;
	    default:
	    	throw "ShapeFunctionT::SetNQ, Unsupported number of Gaussian points";
		}
}

double GreenFunctionT::Gaussian_Weight(int num_q, int q_i)
{
	switch (num_q){

		case 6:
		{
		      double weights[6]={
		    		  	0.3607615730481386,0.3607615730481386,0.4679139345726910,
		                0.4679139345726910,0.1713244923791704,0.1713244923791704};

		      return weights[q_i];
		}
		      break;


	    case 10:
	    {
	           double weights[10]={
	        		   0.295524225,0.295524225,0.269266719,0.269266719,0.219086363,
	        		   0.219086363,0.149451349,0.149451349,0.066671344,0.066671344};

	           return weights[q_i];
	    }

	           break;

	    case 20:
	    {
	    		double weights[20]={
	    				0.152753387,  0.152753387, 0.149172986, 0.149172986,
	    				0.142096109,  0.142096109, 0.131688638, 0.131688638,
	    				0.118194532,  0.118194532, 0.10193012,  0.10193012,
	    				0.083276742,  0.083276742, 0.062672048, 0.062672048,
	    				0.04060143,   0.04060143,  0.017614007, 0.017614007};

	    		return weights[q_i];
	    }
	    		break;
	    case 40:
	    {
	    	double weights[40]={
	    			0.0775059479784248, 0.0775059479784248, 0.0770398181642480, 0.0770398181642480,
	    			0.0761103619006262,	0.0761103619006262, 0.0747231690579683, 0.0747231690579683,
	    			0.0728865823958041, 0.0728865823958041, 0.0706116473912868, 0.0706116473912868,
	    			0.0679120458152339, 0.0679120458152339, 0.0648040134566010, 0.0648040134566010,
	    			0.0613062424929289,	0.0613062424929289, 0.0574397690993916,	0.0574397690993916,
	    			0.0532278469839368,	0.0532278469839368,	0.0486958076350722, 0.0486958076350722,
	    			0.0438709081856733, 0.0438709081856733,	0.0387821679744720,	0.0387821679744720,
	    			0.0334601952825478,	0.0334601952825478,	0.0279370069800234, 0.0279370069800234,
	    			0.0222458491941670,	0.0222458491941670, 0.0164210583819079,	0.0164210583819079,
	    			0.0104982845311528, 0.0104982845311528,	0.0045212770985332,	0.0045212770985332};

	    	return weights[q_i];
	    }
	    	break;

	    default:
	    	throw "ShapeFunctionT::SetNQ, Unsupported number of Gaussian points";
		}
}











