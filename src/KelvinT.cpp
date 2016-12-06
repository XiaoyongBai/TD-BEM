//
//  KelvinT.cpp
//  GreenFunction
//
//  Created by Xiaoyong Bai on 5/13/16.
//
//

#include "KelvinT.h"
#include "MathOperationT.h"
#include "cmath"

using namespace GreenFunction;


KelvinT::KelvinT()
{
    //allocate members
    fR_v=new double[3];
    fN=new double[3];
    
    fDis=new double[9];
    fTrac=new double[9];
    fStress_X=new double[9];
    fStress_Y=new double[9];
    fStress_Z=new double[9];
}


KelvinT::~KelvinT()
{
    delete [] fR_v;
    delete [] fN;
    
    delete [] fDis;
    delete [] fTrac;
    
    delete [] fStress_X;
    delete [] fStress_Y;
    delete [] fStress_Z;
}


void KelvinT::SetMaterial(double E, double nu)
{
    fE=E;
    fNu=nu;
    
    fMu=fE/(2*(1+nu));
}


void KelvinT::SetGeometry(double* src, double* dest, double* normal)
{
    
    for (int i=0; i<3; i++) {
        fR_v[i]=dest[i]-src[i];
        fN[i]=normal[i];
    }
    
    fR=MathOperationT::VecNormal(3, fR_v);
    double nn=MathOperationT::VecNormal(3, fN);
    
    MathOperationT::VecScale(3, fR_v, 1/fR);
    MathOperationT::VecScale(3, fN, 1/nn);

}


void KelvinT::Compute()
{
    ComputeDisplacement();
    ComputeStress();
    ComputeTraction();
}



void KelvinT::ComputeDisplacement()
{
    
    for (int k=0; k<3; k++) {
        for (int i=0; i<3; i++) {
            int dik=KronDelta(i, k);
            fDis[i*3+k]=(3-4*fNu)*dik+fR_v[i]*fR_v[k];
        }
    }
    
    double scale= 1/(16*M_PI*fMu*(1-fNu)*fR);
    MathOperationT::VecScale(9, fDis, scale);
}


void KelvinT::ComputeStress()
{
    
    double factor=-1/(8*M_PI*(1-fNu)*fR*fR);
    
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            int k=1;
            fStress_X[i*3+j]=(1-2*fNu)*(KronDelta(i,k)*fR_v[j]+KronDelta(j,k)*fR_v[i]-KronDelta(i,j)*fR_v[k])
                                +3*fR_v[i]*fR_v[j]*fR_v[k];
            
            k=2;
            fStress_Y[i*3+j]=(1-2*fNu)*(KronDelta(i,k)*fR_v[j]+KronDelta(j,k)*fR_v[i]-KronDelta(i,j)*fR_v[k])
            +3*fR_v[i]*fR_v[j]*fR_v[k];
            
            k=3;
            fStress_Z[i*3+j]=(1-2*fNu)*(KronDelta(i,k)*fR_v[j]+KronDelta(j,k)*fR_v[i]-KronDelta(i,j)*fR_v[k])
            +3*fR_v[i]*fR_v[j]*fR_v[k];
            
            fStress_X[i*3+j]*=factor;
            fStress_Y[i*3+j]*=factor;
            fStress_Z[i*3+j]*=factor;
        }//end loop j
    }//end loop i
    
}

void KelvinT::ComputeTraction()
{
    double* temp=fTrac;
    
    //Compute Kelvin's traction
    
    double rdn=MathOperationT::VecDot(3, fR_v, fN);
    
    double divider=-(8*M_PI)*(1-fNu)*fR*fR;
    for ( int i=0; i<3; i++) {
        for( int k=0; k<3; k++) {
            int dik=KronDelta(i, k);
            double rirk=fR_v[i]*fR_v[k];
            double rink=fR_v[i]*fN[k];
            double rkni=fR_v[k]*fN[i];
            
            *temp = ((1-2*fNu)*dik+3*rirk)*rdn-(1-2*fNu)*(rkni-rink);
            *temp /= divider;
            
            temp++;
        }//end loop k
    }//end loop i
}


int KelvinT::KronDelta(int i, int j)
{
    if (i==j) {
        return 1;
    }
    else{
        return 0;
    }
}













