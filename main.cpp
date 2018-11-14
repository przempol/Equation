#include <iostream>
#include <cmath>
//#include <string>

using namespace std;

const double pi = 3.14159265359;


const unsigned int N = 100;
const double dx = 1./N;
const double dt = 0.0001;


const double kappa=0;
const double omega=0;


double calculateStationary(unsigned int n, unsigned int k){
	return sqrt(2.)*sin(n*pi*k*dx);
}

double calculateHamiltonian(double psi_, double psi0, double psi1, double tau,unsigned int k){
	return (-psi1-psi_+2*psi0)/(2*dx*dx)+kappa*(k*dx-0.5)*psi0*sin(omega*tau);
}


int main(int argc, char* argv[]){
	double* psiR=new double[N+1];
	double* psiI=new double[N+1];

	double* HRtmp=new double[N+1];
	double* HItmp=new double[N+1];
	HRtmp[0]=0.;
	HRtmp[N]=0.;
	HItmp[0]=0.;
	HItmp[N]=0.;
	

	double tau=0;
	
	for(unsigned int ii=0;ii<N+1;ii++){
		psiR[ii]=calculateStationary(1,ii);
		psiI[ii]=0;
	}
	
	for(unsigned int ii=1;ii<N;ii++){
		HRtmp[ii]=calculateHamiltonian(psiR[ii-1], psiR[ii], psiR[ii+1],tau,ii);
		HItmp[ii]=calculateHamiltonian(psiI[ii-1], psiI[ii], psiI[ii+1],tau,ii);
	}
	
	
	
	
	for(unsigned int sim=0;sim<167;sim++){
		//leap frog ?
		
		
		for(unsigned int ii=0;ii<N+1;ii++)	psiR[ii]=psiR[ii]+HItmp[ii]*dt/2.;
		
		for(unsigned int jj=1;jj<N;jj++)	HItmp[jj]=calculateHamiltonian(psiI[jj-1], psiI[jj], psiI[jj+1],tau+dt/2.,jj);
		
		for(unsigned int ii=0;ii<N+1;ii++)	psiI[ii]=psiI[ii]-HRtmp[ii]*dt;	
		
		for(unsigned int jj=1;jj<N;jj++)	HRtmp[jj]=calculateHamiltonian(psiR[jj-1], psiR[jj], psiR[jj+1],tau+dt,jj);
			
		for(unsigned int ii=0;ii<N+1;ii++)	psiR[ii]=psiR[ii]+HItmp[ii]*dt/2.;
		
		
		
		for(unsigned int ii=0;ii<N+1;ii++)	cout<<psiR[ii]<<endl;
		cout<<endl<<endl;
		//leapfrog end
		
	}
	
	
	//for(unsigned int ii=0;ii<N+1;ii++) cout<<psiR[ii]<<endl;
	
	
	delete psiR;
	delete psiI;
	delete HRtmp;
	delete HItmp;
	return  0;
}
