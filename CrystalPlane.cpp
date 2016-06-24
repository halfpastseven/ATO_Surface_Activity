#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#define SMALL 0.0001
using namespace std;

struct atom{
	atom():X(0),Y(0),Z(0), Flags(0){};
	//atom(const atom& A): X(A.X), Y(A.Y),Z(A.Z),Flags(A.Flags){}
	atom(float x, float y, float z):X(x),Y(y),Z(z), Flags(0){}
	float X;
	float Y;
	float Z;
	int Flags;
	inline void displace (float x, float y, float z) {X+=x; Y+=y; Z+=z;}
	inline void displace (const atom& Vect) {X+=Vect.X; Y+=Vect.Y; Z+=Vect.Z;}
};

struct plane {
	plane (int h, int k, int l, const atom& base):H(h),K(k),L(l),Base(base){}
	
	inline float assess (const atom& A) const{
		return int(((A.X-Base.X) * H + (A.Y-Base.Y)*K + (A.Z-Base.Z)*L +SMALL) * 100 );
	}
	
	int H;
	int K;
	int L;
	atom Base; 
};

enum atomType{
	Dop	   = 0x0001,
	In_110 = 0x0002, In_101 = 0x0004, In_200 = 0x0008, In_211 = 0x0010,
	Un_110 = 0x0020, Un_101 = 0x0040, Un_200 = 0x0080, Un_211 = 0x0100,
	Act110 = 0x0200, Act101 = 0x0400, Act200 = 0x0800, Act211 = 0x1000
};

int Xlo = 0;
int Xhi = 50;
int Ylo = 0;
int Yhi = 50;
int Zlo = 0;
int Zhi = 100;

inline int search(float X, float Y, float Z){
	if ((X+SMALL)<Xlo||(X-SMALL)>Xhi||
		(Y+SMALL)<Xlo||(Y-SMALL)>Yhi||
		(Z+SMALL)<Xlo||(Z-SMALL)>Zhi) return -1;
	if (int (X*2 + SMALL) % 2) {// centre atom
		return (Xhi-Xlo+1)*(Yhi-Ylo+1)*(Zhi-Zlo+1) 
				+ (int(X)-Xlo)*(Yhi-Ylo)*(Zhi-Zlo) 
				+ (int(Y)-Ylo)*(Zhi-Zlo)
				+ (int(Z)-Zlo);
	}
	else{
		return (int(X+SMALL)-Xlo)*(Yhi-Ylo+1)*(Zhi-Zlo+1) 
				+ (int(Y+SMALL)-Ylo)*(Zhi-Zlo+1)
				+ (int(Z+SMALL)-Zlo);
	}
}

inline int search (const atom& A) {return search (A.X, A.Y, A.Z);}




/*
Debugger

*/

inline string showBit (int i){
	string c(16,'.');
	for (int j=15; j>=0&&i; j--) {
		if (i%2) {c[j] = '*';}
		i>>=1;
	}
	c.insert(15,"|");
	c.insert(11,"|");
	c.insert(7,"|");
	c.insert(3,"|");
	return c;
}

inline void CheckAtom (const atom& A){
	printf ("Atom %i (%.1f, %.1f, %.1f): \t%s\n",
			search(A),A.X, A.Y, A.Z,
			showBit(A.Flags).c_str());
}

/**** End of Debugger***/
vector<double> MonteCarloSim (float PercentSb){ 
	vector <atom> A(0);
	vector <size_t> Sb(0);
	vector <size_t> inP110(0);
	vector <size_t> inP101(0);
	vector <size_t> inP200(0);
	vector <size_t> inP211(0);
	
	A.reserve ((Xhi-Xlo+1)*(Yhi-Ylo+1)*(Zhi-Zlo+1)*2);
	Sb.reserve((Xhi-Xlo+1)*(Yhi-Ylo+1)*(Zhi-Zlo+1)*2);
	//Prepare the super lattice
	for (int iX = Xlo; iX<=Xhi; iX++){ // Corner atoms
		for (int iY = Ylo; iY<=Yhi; iY++){
			for (int iZ = Zlo; iZ<=Zhi; iZ++){
				A.push_back (atom(iX,iY,iZ));
				if ((rand()%1000000)/10000.0<PercentSb) {
					Sb.push_back(A.size()-1);// cout<<iX<<"\t"<<iY<<"\t"<<iZ<<endl;
					A[A.size()-1].Flags |= Dop;  //Random doping
				}
				//CheckAtom (A[A.size()-1]);
			}
		}
	}

	for (int iX = Xlo; iX<Xhi; iX++){ //Center atoms
		for (int iY = Ylo; iY<Yhi; iY++){
			for (int iZ = Zlo; iZ<Zhi; iZ++){
				A.push_back (atom(iX+0.5,iY+0.5,iZ+0.5));
				if ((rand()%1000000)/10000.0<PercentSb) {
					Sb.push_back(A.size()-1);// cout<<iX+0.5<<"\t"<<iY+0.5<<"\t"<<iZ+0.5<<endl;
					A[A.size()-1].Flags |= Dop;  //Random doping
				}
			}
		}
	}
	
	
	//Prepare the planes
	const atom BaseAtom ( A[search ((Xhi-Xlo)/2.0+Xlo,(Yhi-Ylo)/2.0+Ylo,(Zhi-Zlo)/2.0+Zlo)]);
	/*cout <<"Base : ("<<search ((Xhi-Xlo)/2.0+Xlo,(Yhi-Ylo)/2.0+Ylo,(Zhi-Zlo)/2.0+Zlo)<<") ("
		 <<BaseAtom.X<<" , "<<BaseAtom.Y<<" , "<<BaseAtom.Z<<")"<<endl;*/
	const plane P110 (1,1,0, BaseAtom);
	const plane P101 (1,0,1, BaseAtom);
	const plane P200 (2,0,0, BaseAtom);
	const plane P211 (2,1,1, BaseAtom);
	
	
	//Assess the atoms with respect to plane
	//cout<<"Size of the lattice "<<A.size()<<endl;
	for (int iA=0; iA<A.size(); iA++){
	
		int PositionWrtPlane110 = P110.assess(A[iA]);
		if (!PositionWrtPlane110) {A[iA].Flags|=In_110; inP110.push_back(iA);}
		else if (PositionWrtPlane110<0) A[iA].Flags|=Un_110;	

		int PositionWrtPlane101 = P101.assess(A[iA]);
		if (!PositionWrtPlane101) {A[iA].Flags|=In_101; inP101.push_back(iA);}
		else if (PositionWrtPlane101<0) A[iA].Flags|=Un_101;
		
		int PositionWrtPlane200 = P200.assess(A[iA]);
		if (!PositionWrtPlane200) {A[iA].Flags|=In_200; inP200.push_back(iA);}
		else if (PositionWrtPlane200<0) A[iA].Flags|=Un_200;	
		
		int PositionWrtPlane211 = P211.assess(A[iA]); //printf ("PositionWrtPlane211 : %i\n", PositionWrtPlane211);
		if (!PositionWrtPlane211) {A[iA].Flags|=In_211; inP211.push_back(iA);}
		else if (PositionWrtPlane211<0) A[iA].Flags|=Un_211;	
				
		/*cout<<"("<<A[iA].X<<" ,"<<A[iA].Y<<" ,"<<A[iA].Z<<") : "
				 <<"   (110) "<<PositionWrtPlane110
				 <<"   (101) "<<PositionWrtPlane101
				 <<"   (200) "<<PositionWrtPlane200
				 <<"   (211) "<<PositionWrtPlane211<<"\n";*/
	}
	
	//Assess the activity around Sb
	for (int iSb=0; iSb<Sb.size(); iSb++){
		
		atom Peripheral [10];
		for (int i=0; i<10; i++) Peripheral[i].displace(A[Sb[iSb]]);
		
		//Position all the peripheral atoms
		Peripheral[0].displace (0,0,1); Peripheral[1].displace (0,0,-1); 
		Peripheral[2].displace (-0.5, -0.5, -0.5);Peripheral[3].displace (-0.5, -0.5, +0.5);
		Peripheral[4].displace (-0.5, +0.5, -0.5);Peripheral[5].displace (-0.5, +0.5, +0.5);
		Peripheral[6].displace (+0.5, -0.5, -0.5);Peripheral[7].displace (+0.5, -0.5, +0.5);
		Peripheral[8].displace (+0.5, +0.5, -0.5);Peripheral[9].displace (+0.5, +0.5, +0.5);
		
		A[Sb[iSb]].Flags |= (Act110|Act101|Act200|Act211);  // Assuming the centre Sb is active;
				
		for (int i=0; i<10; i++) {
			int indexPeriph = search (Peripheral[i]);
			if (indexPeriph==-1) continue;
			
			if ((A[Sb[iSb]].Flags & (In_110|Un_110))&&(A[indexPeriph].Flags & (In_110|Un_110))) {
				//If this both are in/under plane
				if (!(A[indexPeriph].Flags & Dop))  A[indexPeriph].Flags |= Act110; //Activate peripheral atom
				else A[Sb [iSb]].Flags &= (-1-Act110);  //De-activate the centre atom
			}
			
			if ((A[Sb[iSb]].Flags & (In_101|Un_101))&&(A[indexPeriph].Flags & (In_101|Un_101))) {
				if (!(A[indexPeriph].Flags & Dop))  A[indexPeriph].Flags |= Act101; //Activate peripheral atom
				else A[Sb [iSb]].Flags &= (-1-Act101);
			}	
			
			if ((A[Sb[iSb]].Flags & (In_200|Un_200))&&(A[indexPeriph].Flags & (In_200|Un_200))) {
				if (!(A[indexPeriph].Flags & Dop))  A[indexPeriph].Flags |= Act200; //Activate peripheral atom
				else A[Sb [iSb]].Flags &= (-1-Act200);
			}	
			
			if ((A[Sb[iSb]].Flags & (In_211|Un_211))&&(A[indexPeriph].Flags & (In_211|Un_211))) {
				if (!(A[indexPeriph].Flags & Dop))  A[indexPeriph].Flags |= Act211; //Activate peripheral atom
				else A[Sb [iSb]].Flags &= (-1-Act211);
			}	
			
		}

	}	
	
	/** Assess the surface activity**/
	
	int Count110 = 0;
	int Count101 = 0;
	int Count200 = 0;
	int Count211 = 0;
	for (int i=0; i<inP110.size(); i++) Count110 += (A[inP110[i]].Flags & Act110) ? 1:0;
	for (int i=0; i<inP101.size(); i++) Count101 += (A[inP101[i]].Flags & Act101) ? 1:0;
	for (int i=0; i<inP200.size(); i++) Count200 += (A[inP200[i]].Flags & Act200) ? 1:0;
	for (int i=0; i<inP211.size(); i++) Count211 += (A[inP211[i]].Flags & Act211) ? 1:0;
	
	vector <double> results (5);
	results[0] = Sb.size()/ double (A.size());
	results[1] = Count110 / double (inP110.size());
	results[2] = Count101 / double (inP101.size());
	results[3] = Count200 / double (inP200.size());
	results[4] = Count211 / double (inP211.size()); //printf ("%i out of %i of 211\n", Count211, inP211.size());

	/*vector<size_t> vA (inP211);
	for (int i=0; i<vA.size(); i++) CheckAtom (A[vA[i]]);*/

		
	return results;
}

int main (int argc, char* argv[]){
	srand (time(NULL));
	time_t Start = time (NULL);
	const int Trials = 50;
	
	//vector <vector<double> > Result (99);
	ofstream OutputTrials("OutputTrials.txt");
	OutputTrials<<"Doping\tA110\tA101\tA200\tA211"<<endl;

	ofstream OutputResults("OutputResults.txt");
	OutputResults<<"Doping\tA110\tA101\tA200\tA211"<<endl;
	
	int TestNo=0;
	vector<vector<double> > Ravg(0); //For storing the average
	for (float Dop=0.1; Dop<100; Dop*=1.2){
		cout<<"Doping: "<<Dop<<"% "<<endl;
		
		vector<vector<double> > Rcheck(0);
		
		for (int iTrial =0; iTrial<Trials; iTrial++){
			vector <double> R = MonteCarloSim (Dop);
			Rcheck.push_back (R);
			printf("%.8f       %.8f  %.8f  %.8f  %.8f \n",R[0],R[1],R[2],R[3],R[4]);
			OutputTrials <<R[0]<<"\t"<<R[1]<<"\t"<<R[2]<<"\t"<<R[3]<<"\t"<<R[4]<<endl;
		}
		Ravg.push_back (vector<double> (5,0));
		
		for (int iR =0; iR<5; iR++) {
			for (int iTrial=0; iTrial<Trials;iTrial++) Ravg [TestNo][iR] += Rcheck[iTrial][iR]/Trials;
		}
		
		cout<<"-------Summary---------"<<endl;
		cout <<Ravg[TestNo][0]<<"\t"<<Ravg[TestNo][1]<<"\t"<<Ravg[TestNo][2]
				<<"\t"<<Ravg[TestNo][3]<<"\t"<<Ravg[TestNo][4]<<endl;
		OutputResults <<Ravg[TestNo][0]<<"\t"<<Ravg[TestNo][1]<<"\t"<<Ravg[TestNo][2]
				<<"\t"<<Ravg[TestNo][3]<<"\t"<<Ravg[TestNo][4]<<endl;
		
		int TellTime = time (NULL)-Start;
		
		printf ("\nTime elapsed: %i s\n", TellTime);
		TestNo++;
	}
	return 0;
}
