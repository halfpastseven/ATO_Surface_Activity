#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <string>
#include <cmath>
//#define SMALL 0.0001
using namespace std;

/* Symbol Conventions:
** X, Y, Z: True coordinates (including lower case)
** A, B, C: Number of lattice parameters in each direction (including lower case) (lo, hi = low, high)
** H, K, L: Miller indices
** newX, newY, newZ: Transformed true coordinates
** a, c: lattice parameters*/

//#define OUTPUT_ATOM_COORDINATES //Enable this to output the coordinates of the atoms in a plane 

/***Lattice***/
const float a(0.4738);
const float c(0.3187);
const float O_coord (0.30478);

int Alo = -50; //Lattice Size, default
int Ahi = 50;
int Blo = -50;
int Bhi = 50;
int Clo = -100;
int Chi = 100;

inline int Round_Int (float x){
	return (x>0) ? (x+0.5) : (x-0.5);
}

struct atom{
	atom():A(0),B(0),C(0), Flags(0){};
	//atom(const atom& Atoms): A(Atoms.A), B(Atoms.B),C(Atoms.C),Flags(Atoms.Flags){}
	atom(float A_, float B_, float C_):A(A_),B(B_),C(C_), Flags(0){}
	atom(const atom& At):A(At.A),B(At.B),C(At.C), Flags(At.Flags){}
	float A;
	float B;
	float C;
	int Flags;
	inline void displace (float A_, float B_, float C_) {A+=A_; B+=B_; C+=C_;}
	inline void displace (const atom& Vect) {A+=Vect.A; B+=Vect.B; C+=Vect.C;}
	inline float Sqr () const {return (A*A+B*B+C*C);}
	inline float Abs () const {return sqrt(A*A+B*B+C*C);}
};

struct plane {
	plane (int h, int k, int l, const atom& base):H(h),K(k),L(l),Base(base){
		HKL_Square = H*H+K*K+L*L;
	}
	
	inline float assessPlane (const atom& At) const{
		return Round_Int(((At.A-Base.A) * H + (At.B-Base.B)*K + (At.C-Base.C)*L) * 100 );
	}
	
	inline float assessSquareDistanceToPlane (const atom& At) const{
		float SqrDis = (At.A-Base.A)*H + (At.B-Base.B)*K + (At.C-Base.C)*L;
		return SqrDis*SqrDis / HKL_Square;
	}

	int H;
	int K;
	int L;
	atom Base; 
	int HKL_Square;
	inline float Sqr () const {return (HKL_Square);}
	inline float Abs () const {return sqrt(HKL_Square);}
};

enum atomType{
	Dop	   = 0x01,
	In_Plane = 0x02, Under_Plane = 0x04,
	Active = 0x08
};

/***Containers***/
vector <atom> Atoms;
vector <size_t> Sb;
vector <atom> O;
atom BaseAtom;

/***Methods***/
int PrepareLattice(int Alo_, int Ahi_, 
		int Blo_, int Bhi_,	int Clo_, int Chi_,
		float PercentSb=0.5);
//Responsible for setting the atoms in lattice

float addSb (float Level){//Responsible for randomly replacing Sn with Sb
	
	const unsigned Resolution = 0x1000000;
	for (int iA = 0; iA<Atoms.size(); iA++){
		if ((rand()%Resolution)/double(Resolution) < Level) {
			Atoms[iA].Flags |= Dop; 
			Sb.push_back(iA);
		}
	}
	return Sb.size()/float(Atoms.size());
}

vector <float> assess(const plane& P);
//Responsible for assessing the activity of planes

bool fSame (float f1, float f2) {return !(Round_Int ((f1-f2)*1000));}
bool isZero (float f1) {return !(Round_Int (f1*1000));}

inline int search(float A, float B, float C){
	if (Round_Int(A)<Alo||Round_Int(A)>Ahi||
		Round_Int(B)<Alo||Round_Int(B)>Bhi||
		Round_Int(C)<Alo||Round_Int(C)>Chi) return -1;
		
	if (Round_Int (A*2) % 2) {// The atom is a centre atom "half"
		return (Ahi-Alo+1)*(Bhi-Blo+1)*(Chi-Clo+1) 
				+ (Round_Int(A-0.5)-Alo)*(Bhi-Blo)*(Chi-Clo) 
				+ (Round_Int(B-0.5)-Blo)*(Chi-Clo)
				+ (Round_Int(C-0.5)-Clo); 
	}
	else{// The atom is a corner atom "integral"
		return (Round_Int(A)-Alo)*(Bhi-Blo+1)*(Chi-Clo+1) 
				+ (Round_Int(B)-Blo)*(Chi-Clo+1)
				+ (Round_Int(C)-Clo);
	}
}

inline int search (const atom& At) {
	return search (At.A, At.B, At.C);
}

inline string showBit (int i){
	const int lengthOfBit=4;
	string c(lengthOfBit,'.');
	for (int j=lengthOfBit-1; j>=0 && i; j--, i>>=1) if (i%2) c[j] = '*';
	return c;
}

inline void CheckAtom (const atom& At){
	printf ("Atom %i (%.5f, %.5f, %.5f): \t%s\n",
			search(At),At.A, At.B, At.C,
			showBit(At.Flags).c_str());
}

void resetAtomFlags (bool resetSb=0){
	for (int iA=0; iA<Atoms.size(); iA++) {
		Atoms[iA].Flags &= (-1-In_Plane); //clear a flag without affecting other flags
		Atoms[iA].Flags &= (-1-Under_Plane);
		Atoms[iA].Flags &= (-1-Active);
		if (resetSb) Atoms[iA].Flags &= (-1-Dop);
	}
	if (resetSb) Sb.resize(0);
}

int  PrepareLattice (
		int Alo_, int Ahi_, 
		int Blo_, int Bhi_,
		int Clo_, int Chi_, float PercentSb){ 

	Alo=Alo_; Ahi=Ahi_; Blo=Blo_; Bhi=Bhi_; Clo=Clo_; Chi=Chi_; //set pars
	
	Atoms.reserve ((Ahi-Alo+1)*(Bhi-Blo+1)*(Chi-Clo+1)*2); //Guaranteed enough
	Sb.reserve((Ahi-Alo+1)*(Bhi-Blo+1)*(Chi-Clo+1)*2);
	O.reserve ((Ahi-Alo+1)*(Bhi-Blo+1)*(Chi-Clo+1)*4); 

	//Prepare the super lattice
	for (int iA = Alo; iA<=Ahi; iA++){ // Corner Sn atoms
		for (int iB = Blo; iB<=Bhi; iB++){
			for (int iC = Clo; iC<=Chi; iC++){
				Atoms.push_back (atom(iA,iB,iC));
			}
		}
	}

	for (int iA = Alo; iA<Ahi; iA++){ //Center atoms
		for (int iB = Blo; iB<Bhi; iB++){
			for (int iC = Clo; iC<Chi; iC++){
				Atoms.push_back (atom(iA+0.5,iB+0.5,iC+0.5));
			}
		}
	}
	
	addSb (PercentSb);	
	
	for (int iA = Alo; iA<Ahi; iA++){ 
		for (int iB = Blo; iB<Bhi; iB++){
			for (int iC = Clo; iC<=Chi; iC++){
				
				O.push_back (atom(iA    +O_coord,	iB    +O_coord,	iC)); //(b,b,0)
				O.push_back (atom(iA +1 -O_coord,	iB +1 -O_coord,	iC)); //(1-b,1-b,0)

				if (iA<Ahi && iB<Bhi && iC<Chi){ //Make sure atom does not go beyond the boundary
					O.push_back (atom(iA    +O_coord,	iB +1 -O_coord,	iC+0.5)); //(b, 1-b, 0.5)
					O.push_back (atom(iA +1 -O_coord,	iB    +O_coord, iC+0.5)); //(1-b, b, 0.5)
					O.push_back (atom(iA    +O_coord,	iB    +O_coord,	iC+1)); //(b,b,1)
					O.push_back (atom(iA +1 -O_coord,	iB +1 -O_coord,	iC+1)); //(1-b,1-b,1)
				}
				//CheckAtom (Atoms[Atoms.size()-1]);
			}
		}
	}	
	//Prepare the planes
	return Atoms.size();
}

vector<float> assess(const plane& P){ //Assess the atom flags with respect to plane

	vector <size_t> atomsInPlane;
	
	for (int iA=0; iA<Atoms.size(); iA++){
		int PositionToPlane = P.assessPlane(Atoms[iA]);
		if (!PositionToPlane) {
			Atoms[iA].Flags|=In_Plane; 
			atomsInPlane.push_back(iA);
		}
		else if (PositionToPlane<0) Atoms[iA].Flags|=Under_Plane;
	}

	//Assess the activity around Sb.
	int Cnt = 0;
	for (int iSb=0; iSb<Sb.size(); iSb++){
		
		const float Threshold = 1.1;
		/* Add an assessment of atom with respect to the distance of the plane */
		float SquareDistance = P.assessSquareDistanceToPlane (Atoms[Sb[iSb]]);
		
		if (SquareDistance > Threshold*Threshold ||
				(!(Atoms[Sb[iSb]].Flags & (In_Plane|Under_Plane)))) continue;

		const size_t NofPerihperals(10); //Number of peripheral Sn/Sb atoms
		
		atom Peripheral [NofPerihperals];
		for (int i=0; i<NofPerihperals; i++) Peripheral[i].displace(Atoms[Sb[iSb]]);
		
		//Position all the peripheral atoms
		Peripheral[0].displace (0,0,1); Peripheral[1].displace (0,0,-1); 
		Peripheral[2].displace (-0.5, -0.5, -0.5);Peripheral[3].displace (-0.5, -0.5, +0.5);
		Peripheral[4].displace (-0.5, +0.5, -0.5);Peripheral[5].displace (-0.5, +0.5, +0.5);
		Peripheral[6].displace (+0.5, -0.5, -0.5);Peripheral[7].displace (+0.5, -0.5, +0.5);
		Peripheral[8].displace (+0.5, +0.5, -0.5);Peripheral[9].displace (+0.5, +0.5, +0.5);
				
		Atoms[Sb[iSb]].Flags |= Active;  // Assuming the centre Sb is active;
				
		for (int i=0; i<NofPerihperals; i++) {
			int indexPeriph = search (Peripheral[i]);
			if (indexPeriph==-1) continue;  //Searching out of boundary of the supercell

			if ((Atoms[indexPeriph].Flags & (In_Plane|Under_Plane))) {
				//If this both are in/under plane
				if (!(Atoms[indexPeriph].Flags & Dop)) {
					Atoms[indexPeriph].Flags |= Active; //Activate peripheral atom
				}
				else {
					Atoms[Sb [iSb]].Flags &= (-1-Active);   //De-activate the centre atom
				}
			}
			//printf ("\n");
		}
	}	
	
	/** Assess the surface activity**/
	
	int CountActive = 0;
	for (int i=0; i<atomsInPlane.size(); i++) {
		CountActive += (Atoms[atomsInPlane[i]].Flags & Active)? 1:0; //If active, then add 1
	}
	
	vector <float> results (2);
	results[0] = Sb.size()/ double (Atoms.size());
	results[1] = CountActive / double (atomsInPlane.size());
	
#ifdef OUTPUT_ATOM_COORDINATES
	
	/*Real coordinates X, Y, Z will be involved here*/
	
	string OutputFilename ("SurfaceAtoms_");
	OutputFilename.push_back ('0'+P.H);
	OutputFilename.push_back ('0'+P.K);
	OutputFilename.push_back ('0'+P.L);
	OutputFilename = OutputFilename + ".txt";
	ofstream ExportPlane (OutputFilename.c_str());
	
	atom newX;
	atom newY;
	atom newZ (P.H/a, P.K/a, P.L/c);
	float Normalizer = newZ.Abs();
	newZ.A/=Normalizer; newZ.B/=Normalizer; newZ.C/=Normalizer;	

	float SmallestDistance=10;
	for (int iA=-1; iA<=1; iA++){
		for (int iB=-1; iB<=1; iB++){
			for (int iC=-1; iC<=1; iC++){
				if (isZero(iA*P.H+iB*P.K+iC*P.L) && (iA||iB||iC)) {
					float X = iA*a;
					float Y = iB*a;
					float Z = iC*c;
					
					float Distance = sqrt(X*X + Y*Y + Z*Z);
					if (Distance <SmallestDistance) {
						newX = atom(X, Y, Z);
						SmallestDistance = Distance;
					}
				}
			}
		}
	}
	
	for (int iA=-1; iA<=0; iA++){
		for (int iB=-1; iB<=0; iB++){
			for (int iC=-1; iC<=0; iC++){
				if (isZero((iA+0.5)*P.H+(iB+0.5)*P.K+(iC+0.5)*P.L)) {
					float X = (iA+0.5)*a;
					float Y = (iB+0.5)*a;
					float Z = (iC+0.5)*c;
					
					float Distance = sqrt(X*X + Y*Y + Z*Z);
					if (Distance <SmallestDistance) {
						newX = atom(X, Y, Z);
						SmallestDistance = Distance;
					}
				}
			}
		}
	}
	
	//newX = atom (0,1,0); //Enable manually overriding the newX setting

	Normalizer = newX.Abs();
	newX.A /= Normalizer;
	newX.B /= Normalizer;
	newX.C /= Normalizer;

	newY = atom(newZ.B * newX.C - newZ.C * newX.B,
				newZ.C * newX.A - newZ.A * newX.C,
				newZ.A * newX.B - newZ.B * newX.A);
	
	CheckAtom (newX);
	CheckAtom (newY);
	CheckAtom (newZ);
	CheckAtom (BaseAtom);
	printf ("%.5f %.5f %.5f\n", newX.Abs(), newY.Abs(), newZ.Abs());
	printf ("%.5f %.5f %.5f\n", newX.A*newY.A + newX.B * newY.B + newX.C * newY.C,
								newY.A*newZ.A + newY.B * newZ.B + newY.C * newZ.C,
								newX.A*newZ.A + newX.B * newZ.B + newX.C * newZ.C);
	vector <atom> Transf (atomsInPlane.size());
	//Transform to new XYZ
	
	/*Transform Matrix 
	 * [newX.A	newX.B	newX.C ][ X ]
	 * [newY.A	newY.B	newY.C ][ Y ]
	 * [newZ.A	newZ.B	newZ.C ][ Z ]
	 */

	for (size_t iSurfAtom =0; iSurfAtom<atomsInPlane.size(); iSurfAtom++){
		atom AtomTC = Atoms[atomsInPlane[iSurfAtom]];
		AtomTC.displace (-BaseAtom.A, -BaseAtom.B, -BaseAtom.C);
		AtomTC.A *=a;
		AtomTC.B *=a;
		AtomTC.C *=c;		//TC is the relative coord to BaseAtom in true XYZ
		
		Transf [iSurfAtom] = atom (
			newX.A * AtomTC.A + newX.B * AtomTC.B + newX.C * AtomTC.C,
			newY.A * AtomTC.A + newY.B * AtomTC.B + newY.C * AtomTC.C,
			newZ.A * AtomTC.A + newZ.B * AtomTC.B + newZ.C * AtomTC.C
		);
		ExportPlane	<<((Atoms[atomsInPlane[iSurfAtom]].Flags&Dop) ? 'B' : 'N')<<"\t"
					<<(Atoms[atomsInPlane[iSurfAtom]].A-BaseAtom.A)*a<<"\t"
					<<(Atoms[atomsInPlane[iSurfAtom]].B-BaseAtom.B)*a<<"\t"
					<<(Atoms[atomsInPlane[iSurfAtom]].C-BaseAtom.C)*c<<"\t"
					//<<AtomTC.Abs()<<"\t"
					<<Transf [iSurfAtom].A<<"\t"
					<<Transf [iSurfAtom].B<<"\t"
					<<Transf [iSurfAtom].C<<"\t"
					//<<Transf [iSurfAtom].Abs()<<"\t"
					<<endl;
		//CheckAtom (Atoms[atomsInPlane[iSurfAtom]]);
	}
	
	vector <atom> OTransf (O.size());
	for (size_t iSurfO =0; iSurfO<O.size(); iSurfO++){
		atom AtomTC = O[iSurfO];
		AtomTC.displace (-BaseAtom.A, -BaseAtom.B, -BaseAtom.C);
		AtomTC.A *=a;
		AtomTC.B *=a;
		AtomTC.C *=c;		//TC is the relative coord to BaseAtom in true XYZ
		
		OTransf [iSurfO] = atom (
			newX.A * AtomTC.A + newX.B * AtomTC.B + newX.C * AtomTC.C,
			newY.A * AtomTC.A + newY.B * AtomTC.B + newY.C * AtomTC.C,
			newZ.A * AtomTC.A + newZ.B * AtomTC.B + newZ.C * AtomTC.C
		);
		if ((OTransf [iSurfO].C * OTransf [iSurfO].C) < 0.001)
			ExportPlane	<<'O'<<"\t"
					<<(AtomTC.A-BaseAtom.A)*a<<"\t"
					<<(AtomTC.B-BaseAtom.B)*a<<"\t"
					<<(AtomTC.C-BaseAtom.C)*c<<"\t"

					<<OTransf [iSurfO].A<<"\t"
					<<OTransf [iSurfO].B<<"\t"
					<<OTransf [iSurfO].C<<"\t"
					<<endl;
	}
	
#endif
	return results;
}

int main (int argc, char* argv[]){
	
	if (argc<2) {
		printf("Insufficient arguments");
		return 1;
	}
	srand (time(NULL));
	time_t Start = time (NULL);//Time count
	
	stringstream InputHKL (string(argv[1],3));
	
	int H = argv[1][0]-'0';
	int K = argv[1][1]-'0';
	int L = argv[1][2]-'0';
	
	float SbLevel = 0;
	
	BaseAtom = atom (int ((Ahi+Alo)/2), int ((Bhi+Blo)/2), int ((Chi+Clo)/2)); 
	//!!!Careful: Choice of a base atom can affect the arrangement of atoms

	plane TestPlane (H, K, L, BaseAtom);
	PrepareLattice ( Alo, Ahi, Blo, Bhi, Clo, Chi, 0);
	
	/***Prepare the output files***/
	
	string Output ("Output");
	ofstream OutputTrials((Output+argv[1]+".txt").c_str());
	ofstream OutputAverage((Output+argv[1]+"_Average.txt").c_str());
	
	OutputTrials<<"Doping\tActivity"<<endl;

	vector <vector<float> > Result ; //Rows = Sb level, Cols = Trials
	vector<float> Ravg; //For storing the average activity
	vector<float> Rstd; //For storing the standard deviation;
	const size_t Trials = 100;
	
	OutputAverage<<"SbLevel\t"<<"Activity_"<<argv[1]<<"\tStd_"<<argv[1]<<endl;
	for (float SbLevel=0.001; SbLevel<1; SbLevel*=1.1){
		
		
		Result.push_back (vector<float> (Trials));
		
		Ravg.push_back (0);
		Rstd.push_back (0);
		OutputTrials<<SbLevel<<'\t';
		for (size_t iTrial =0; iTrial<Trials; iTrial++){
			resetAtomFlags(true);
			addSb (SbLevel);
			vector <float> R = assess (TestPlane);
			Result[Result.size()-1][iTrial] = R[1];
			Ravg [Ravg.size()-1] += R[1] / Trials;
			Rstd [Rstd.size()-1] += R[1] * R[1];
			OutputTrials<<R[1]<<"\t";
		}
		Rstd [Rstd.size()-1] = sqrt(
				(Rstd [Rstd.size()-1] - Ravg [Rstd.size()-1]*Ravg [Rstd.size()-1]*Trials)
				/(Trials-1));
				
		OutputTrials << endl;
		OutputAverage<< SbLevel << "\t" <<Ravg[Ravg.size()-1] << "\t" <<Rstd[Rstd.size()-1]<<endl;

		printf ("Sb Level : %.5f,  Activity : %.5f +/- %.5f", SbLevel, Ravg[Ravg.size()-1], Rstd[Rstd.size()-1]);

		int TellTime = time (NULL)-Start;
		Start = time (NULL);
		printf ("\nTime elapsed: %i s\n", TellTime);
	}

	return 0;
}
