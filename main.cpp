#include <iostream>
#include <cmath>
#include <fstream>

#include <string>

using namespace std;

const double pi = 3.14159265359;


const unsigned int N = 100;
const double dx = 1./N;
const double dt = 1.00025e-4;
const double simTime = 20.;



unsigned int stat=1;

const double kappa=0.;
const double omega=3.*pi*pi/2.;


double calculateStationary(unsigned int n, unsigned int k){
	return sqrt(2.)*sin(n*pi*k*dx);
}

double calculateHamiltonian(double psi_, double psi0, double psi1, double tau,unsigned int k){
	return (-psi1-psi_+2*psi0)/(2*dx*dx)+kappa*(k*dx-0.5)*psi0*sin(omega*tau);
}

double calculateNorm(double* psiI,double* psiR){
	double ret=1.;
	for(unsigned int ii=0;ii<N+1;ii++)	ret+=(psiR[ii]*psiR[ii]+psiI[ii]*psiI[ii]);
	ret*=dx;
	return ret;
}

double calculateAvPos(double* psiI,double* psiR){ // calculate average position
	double ret=1.;
	for(unsigned int ii=0;ii<N+1;ii++)	ret+=(ii*dx*(psiR[ii]*psiR[ii]+psiI[ii]*psiI[ii]));
	ret*=dx;
	return ret;
}

double calculateAvEnergy(double* psiI,double* psiR,double* HI,double* HR){ // calculate average energy
	double ret=1.;
	for(unsigned int ii=0;ii<N+1;ii++)	ret+=(psiR[ii]*HR[ii]+psiI[ii]*HI[ii]);
	ret*=dx;
	return ret;
}



void saveData(string name, double* data, unsigned number){
	const char *namechar = name.c_str();
	ofstream data_file;
	data_file.open(namechar);
	data_file<<endl;
	data_file.close();
	
	data_file.open(namechar, ios::app);
	//for(unsigned int ii=0;ii<number;ii++)	data_file<<data[ii]<<endl;
	data_file<<"t"<<"\t"<<"V"<<"\t"<<"Ek"<<"\t"<<"Ec"<<"\t"<<"T"<<"\t"<<"P"<<endl;
	data_file.close();
}



int main(int argc, char* argv[]){
	ofstream density_file;
	ofstream energy_file;
	ofstream norm_file;
	ofstream pos_file;
	density_file.open("anim/density.dat");
	energy_file.open("anim/energy.dat");
	norm_file.open("anim/norm.dat");
	pos_file.open("anim/pos.dat");
	
	double* psiR=new double[N+1];
	double* psiI=new double[N+1];

	double* HamR=new double[N+1];
	double* HamI=new double[N+1];
	HamR[0]=0.;
	HamR[N]=0.;
	HamI[0]=0.;
	HamI[N]=0.;
	
	
	double* energyArray=new double[int(( simTime / dt )/100)+1];
	
	
	double sumEnergy=0;
	double varEnergy=0;
	double tau=0;
	
	for(unsigned int ii=0;ii<N+1;ii++){
		psiR[ii]=calculateStationary(stat,ii);
		psiI[ii]=0;
	}
	
	for(unsigned int ii=1;ii<N;ii++){
		HamR[ii]=calculateHamiltonian(psiR[ii-1], psiR[ii], psiR[ii+1],tau,ii);
		HamI[ii]=calculateHamiltonian(psiI[ii-1], psiI[ii], psiI[ii+1],tau,ii);
	}
	
	
	for(unsigned int sim=0;sim<simTime/dt;sim++){		//leap frog 		
		for(unsigned int ii=0;ii<N+1;ii++)	psiR[ii]=psiR[ii]+HamI[ii]*dt/2.;
		
		for(unsigned int jj=1;jj<N;jj++)	HamR[jj]=calculateHamiltonian(psiR[jj-1], psiR[jj], psiR[jj+1],tau+dt/2.,jj);
		
		for(unsigned int ii=0;ii<N+1;ii++)	psiI[ii]=psiI[ii]-HamR[ii]*dt;	
		
		for(unsigned int jj=1;jj<N;jj++)	HamI[jj]=calculateHamiltonian(psiI[jj-1], psiI[jj], psiI[jj+1],tau+dt,jj);
			
		for(unsigned int ii=0;ii<N+1;ii++)	psiR[ii]=psiR[ii]+HamI[ii]*dt/2.; 
		tau+=dt;
		//leapfrog end
		
		if(sim%100==0)	{
			for(unsigned int ii=0;ii<N+1;ii++) density_file<<psiR[ii]*psiR[ii]+psiI[ii]*psiI[ii]<<endl;
			density_file<<endl<<endl;
			
			energyArray[sim/100]=calculateAvEnergy(psiI,psiR,HamI,HamR);
			
			energy_file<<tau-dt<<"\t"<<energyArray[sim/100]<<endl;
			norm_file<<tau-dt<<"\t"<<calculateNorm(psiI,psiR)<<endl;
			pos_file<<tau-dt<<"\t"<<calculateAvPos(psiI,psiR)<<endl;
			
			sumEnergy+=energyArray[sim/100];
		}
	}	
	density_file.close();
	energy_file.close();
	norm_file.close();
	pos_file.close();
	
	
	sumEnergy/=(( simTime / dt )/100);
	
	for(unsigned int ii=0;ii<( simTime / dt )/100;ii++)	varEnergy+=(sumEnergy-energyArray[ii])*(sumEnergy-energyArray[ii]);
	varEnergy/=(( simTime / dt )/100 - 1);
	varEnergy=sqrt(varEnergy);
	cout<<varEnergy<<endl;

	
	delete[] psiR;
	delete[] psiI;
	delete[] HamR;
	delete[] HamI;
	delete[] energyArray;
	return  0;
}
