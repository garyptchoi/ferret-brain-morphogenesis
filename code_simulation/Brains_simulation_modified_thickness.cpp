// Modifying the thickness H (0.5x) in the brain growth simulation
#include <omp.h>
#include "vema.h"
#include "eig3.h"
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
using namespace std;

class Tetra{
public:
	int n1, n2, n3, n4;
	Matrix G; // Growth tensor
	Tetra(): n1(0), n2(0), n3(0), n4(0) {}
	Tetra(int n1_, int n2_, int n3_, int n4_): n1(n1_), n2(n2_), n3(n3_), n4(n4_) {}
};
class Face{
public:
	int n1, n2, n3;
	Face(): n1(0), n2(0), n3(0) {}
	Face(int n1_, int n2_, int n3_): n1(n1_), n2(n2_), n3(n3_) {}
};

void createNNLtriangle(vector<int>*, Vector*, vector<Face>&, int*, int, int, double, double);
Vector closestPointTriangle(Vector&, Vector&, Vector&, Vector&, double&, double&, double&);
void dist2surf(Vector*, int*, int, int, int*, double*);
void writetxt(Vector*, vector<Face>&, int*, int*, int, int);
void writeoff(Vector*, vector<Face>&, int*, int*, int, int);

double ub, vb, wb; // Barycentric coordinates

int main(int argc, char* argv[]) {

	char inpa[10], inpb[10], inpc[10], inpd[10], inpe[10];
	ifstream filu;
	
	// note: should rescale the mesh so that the major radius = 1
	filu.open("P8-tet.mesh"); 
	const double a = 0.005; // Mesh spacing - set manually based on the average spacing in the mesh or smaller
	const int di = 500; // Output data once every di steps 
	// Read nodes
	filu >> inpa;
	int nn = atoi(inpa);
	Vector* Ut = new Vector[nn]();
	Vector* Ut0 = new Vector[nn]();
	for (int i = 0; i < nn; i++) {
		filu >> inpa >> inpb >> inpc;
		Ut0[i] = Vector(atof(inpa), atof(inpb), atof(inpc));
		Ut[i] = Ut0[i];
	}
	cout << "Number of nodes = nn = " << nn << "\n";
	if(nn == 0){
		cout << "Error in reading the mesh file!!\n";
		return 1;
	}
	// Read elements
	filu >> inpa;
	int ne = atoi(inpa);
	Tetra* tets = new Tetra[ne]();
	for (int i = 0; i < ne; i++) {
		filu >> inpa >> inpb >> inpc >> inpd >> inpe;
		tets[i] = Tetra(atoi(inpb)-1, atoi(inpc)-1, atoi(inpd)-1, atoi(inpe)-1); // ****** note orientation
	}
	cout << "Number of tets = ne = " << ne << "\n";
	// Read faces
	filu >> inpa;
	int nf = atoi(inpa);
	vector<Face> faces;
	for (int i = 0; i < nf; i++) {
		filu >> inpa >> inpb >> inpc >> inpd;
		faces.push_back(Face(atoi(inpb)-1, atoi(inpc)-1, atoi(inpd)-1));
	}
	cout << "Number of faces = nf = " << nf << "\n";
	filu.close();
		
	// Determine surface nodes (SN = global indices of the surface nodes, SNb = surface indices of global nodes)
	int nsn = 0; // Number of nodes at the surface
	int* SNb = new int[nn]; // Nodal index map from full mesh to surface
	for (int i = 0; i < nn; i++) SNb[i] = 0;
	for (int i = 0; i < nf; i++) { SNb[faces[i].n1] = 1; SNb[faces[i].n2] = 1; SNb[faces[i].n3] = 1; }
	for (int i = 0; i < nn; i++) if (SNb[i] == 1) nsn++;
	cout << "Number of surface nodes = nsn = " << nsn << endl;
	int* SN = new int[nsn]; // Nodal index map from surface to full mesh
	int p = 0; // Iterator
	for (int i = 0; i < nn; i++) if (SNb[i] == 1) { SN[p] = i; SNb[i] = p; p++; }
	
	cout << "Define parameters..." << "\n";
	double H = 0.1; // Thickness of growing layer
	double Hcp; // Cortical plate thickness for visualization
	const double mu = 1.0; // Shear modulus
	const double K = 5.0; // Bulk modulus 
	const double rho = 0.025; // Mass density - adjust to run the simulation faster or slower 
	const double gamma = 0.5; // Damping constant
	const double dt = 0.05*sqrt(rho*a*a/K); // Time step
	
	const double mpy = 0; // Midplane position, for separating the two sides (y-coord)
	
	double at; // Relative growth
	
	const double alfagt = 1.9;// Tangential growth rate of gray matter
	const double alfagn = 0;
	const double alfaw = 0;
	const double alfar = 0;
	
	int step = 0; // Current timestep
	double t = 0.0; // Current time
	Vector* Vt = new Vector[nn](); // Velocities
	Vector* Ft = new Vector[nn](); // Forces
	double* Vn0 = new double[nn]; // Nodal volumes in reference state
	double* Vn = new double[nn]; // Deformed nodal volumes
	Matrix I; // Unit matrix
	double Ue, Area; // Elastic energy
	
	double zoom = 1.0; // Zoom variable for visualization
	vector<int>* NNLt = new vector<int>[nsn];
	
	double maxDist; // Maximum displacement since the last update of the proximity list
	Vector* Utold = new Vector[nsn]; // Stores positions when proximity list is updated
	const double hs = 0.6*a; // Thickness of contact proximity skin // 0.7*a
	const double hc = 0.2*a; // Thickness of repulsive skin // 0.3*a
	double ub, vb, wb; // Barycentric coordinates
	int tp; // Loop iterator
	
		
	int ninv; // count number of inverted elements
		
	// Find closest surface nodes (csn) and distances to them (d2s)
	int* csn = new int[nn];	// Nearest surface nodes
	double* d2s = new double[nn];	// Distances to nearest surface nodes
	dist2surf(Ut, SN, nn, nsn, csn, d2s);
	Vector* N0 = new Vector[nsn];// Normals in reference state
	Vector Ntmp;
	for (int i = 0; i < nf; i++) {// Find normals
		Ntmp = (Ut0[faces[i].n2] - Ut0[faces[i].n1]).cross(Ut0[faces[i].n3] - Ut0[faces[i].n1]);
		N0[SNb[faces[i].n1]] += Ntmp;
		N0[SNb[faces[i].n2]] += Ntmp;
		N0[SNb[faces[i].n3]] += Ntmp;
	}
	for (int i = 0; i < nsn; i++) N0[i].normalize();
	

	// Mark non-growing areas
	double* gr = new double[nn];
	int count_nongrow = 0;
	for (int i = 0; i < nn; i++) {
		Vector qp = Ut0[i];
//		double rqp = Vector(0, qp.y, qp.z).length(); 
		double rqp = Vector(0, qp.y, qp.z+0.2).length(); 
		gr[i] = 1.0;
		
		// depending on the brain structure, enforce non-growing constraints
		if (qp.x >0 && rqp < 0.5){
			// for fixing the middle part
			gr[i] = gr[i] - 5.0*(0.5-rqp);
			count_nongrow++;
		}
		if (qp.z < -0.4) {
			// fixing the bottom part
			gr[i] = gr[i] + 5.0*(0.4+qp.z);	
			count_nongrow++;
		}
		gr[i] = max(gr[i], 0.0);
		
		
	}
	cout << "# non-growing node = " << count_nongrow << endl; 
	
		
	// Check minimum and maximum edge lengths at the surface
	double mine = 1e9, maxe = 0.0, ave = 0.0;
	for (int i = 0; i < nf; i++) {
		mine = min((Ut[faces[i].n2] - Ut[faces[i].n1]).length(), mine);
		mine = min((Ut[faces[i].n3] - Ut[faces[i].n1]).length(), mine);
		mine = min((Ut[faces[i].n3] - Ut[faces[i].n2]).length(), mine);
		maxe = max((Ut[faces[i].n2] - Ut[faces[i].n1]).length(), maxe);
		maxe = max((Ut[faces[i].n3] - Ut[faces[i].n1]).length(), maxe);
		maxe = max((Ut[faces[i].n3] - Ut[faces[i].n2]).length(), maxe);
		ave += (Ut[faces[i].n3] - Ut[faces[i].n2]).length() + (Ut[faces[i].n3] - Ut[faces[i].n1]).length() + (Ut[faces[i].n2] - Ut[faces[i].n1]).length();
	}
	ave /= 3.0*nf;
	cout << "Edge length information:" << endl;
	cout << "min edge =" << mine  << ", average edge = " << ave << ", max edge = " << maxe << ",	user input = " << a << endl;

	// Find center of mass and dimension of the mesh
	double maxx = -1e9, minx = 1e9, maxy = -1e9, miny = 1e9, maxz = -1e9, minz = 1e9;
	Vector cog;
	for (int i = 0; i < nn; i++) {
		maxx = max(maxx, Ut[i].x); minx = min(minx, Ut[i].x);
		maxy = max(maxy, Ut[i].y); miny = min(miny, Ut[i].y);
		maxz = max(maxz, Ut[i].z); minz = min(minz, Ut[i].z);
		cog += Ut[i];
	}
	cog /= nn;
		
	double maxd = max(max(max(abs(maxx-cog.x), abs(minx-cog.x)), max(abs(maxy-cog.y), abs(miny-cog.y))), max(abs(maxz-cog.z), abs(minz-cog.z)));
	
	cout << "center of mass = (" << cog.x << ", " << cog.y << ", " << cog.z << ") " << endl;
	cout << "Min and Max coordinates:" << endl;
	cout << "x-coord: " << minx << " " << maxx << endl;
	cout << "y-coord: " << miny << " " << maxy << endl;
	cout << "z-coord: " << minz << " " << maxz << endl;
	cout << "Maximum diameter = " << maxd << endl;

	cout << ub << " " << vb << " " << wb << endl;
	
	cout.precision(3);
	ofstream dat;
	dat.precision(5);
	dat.open("Brains_log.dat");
	
	cout << "Main loop started.\n";	
	cout << "       step          t         Ue       Area     Volume       ninv" << endl;		
	while (t < 1.0) { // Main loop
			
		// Folding simulation
		H = 0.5*(0.1 - 0.05*t); // reduce the thickness
		
		// Contact processing for self avoidance
		maxDist = 0.0;
		#pragma omp parallel for reduction(max:maxDist)
		for (int i = 0; i < nsn; i++ ) {
			maxDist = max(maxDist, (Ut[SN[i]] - Utold[i]).length()); // find the maximum displacement of surface nodes
		}		
		
		if (maxDist > 0.5*(hs-hc)) { // Update proximity lists
		
			cout << step << endl;
			createNNLtriangle(NNLt, Ut, faces, SN, nsn, nf, hs, 8.0*a/5);
			
			for (int i = 0; i < nsn; i++) Utold[i] = Ut[SN[i]];
		}
		#pragma omp parallel for private(tp)
		for (int i = 0; i < nsn; i++) { // Determine the closest triangles to points
			for (tp = 0; tp < NNLt[i].size(); tp++) {// Loop trough triangle-proximity list
				int pt = SN[i];
				int tri = NNLt[i][tp];
				// find the triangle closest to the point
				Vector cc = closestPointTriangle(Ut[pt], Ut[faces[tri].n1], Ut[faces[tri].n2], Ut[faces[tri].n3], ub, vb, wb) - Ut[pt];
				double rc = cc.length(); // obtain the spacing
				if (rc < hc){ // Calculate contact force if within the contact range 
					cc.normalize();
					
					Vector Ntri = (Ut[faces[tri].n2] - Ut[faces[tri].n1]).cross(Ut[faces[tri].n3] - Ut[faces[tri].n1]); 
					Ntri.normalize(); 
					
					
					Vector fn = cc*(rc-hc)/hc*K*a*a;
					
					if (fn.dot(Ntri) < 0.0) fn -= Ntri*fn.dot(Ntri)*2.0; // restore the node to the right side of the triangle
					
					Ft[faces[tri].n1] -= fn*ub;
					Ft[faces[tri].n2] -= fn*vb;
					Ft[faces[tri].n3] -= fn*wb;
					Ft[pt] += fn;
				}
			}
		}
		
		
		// Nodal volumes for nodal pressure calculation
		#pragma omp parallel for
		for (int i = 0; i < nn; i++) { Vn0[i] = 0.0; Vn[i] = 0.0; }
		#pragma omp parallel for
		for (int i = 0; i < ne; i++) {
			int n1 = tets[i].n1;
			int n2 = tets[i].n2;
			int n3 = tets[i].n3;
			int n4 = tets[i].n4;
			
			// Undeformed
			Vector xr1 = Ut0[n2] - Ut0[n1];
			Vector xr2 = Ut0[n3] - Ut0[n1];
			Vector xr3 = Ut0[n4] - Ut0[n1];
			Matrix Ar = Matrix(xr1, xr2, xr3);
			Ar = tets[i].G.prod(Ar);
			
			double vol0 = Ar.det()/6.0;
			Vn0[n1] += vol0/4.0;
			Vn0[n2] += vol0/4.0;
			Vn0[n3] += vol0/4.0;
			Vn0[n4] += vol0/4.0;
			
			// Deformed
			Vector x1 = Ut[n2] - Ut[n1];
			Vector x2 = Ut[n3] - Ut[n1];
			Vector x3 = Ut[n4] - Ut[n1];
			Matrix A = Matrix(x1, x2, x3);
			double vol = A.det()/6.0;
			Vn[n1] += vol/4.0;
			Vn[n2] += vol/4.0;
			Vn[n3] += vol/4.0;
			Vn[n4] += vol/4.0;
		}
		
		// Deformations
		Ue = 0.0; 
		ninv = 0;	
		#pragma omp parallel for reduction(+:Ue)
		for (int i = 0; i < ne; i++) {
			int n1 = tets[i].n1;
			int n2 = tets[i].n2;
			int n3 = tets[i].n3;
			int n4 = tets[i].n4;
			
			Vector xr1 = Ut0[n2] - Ut0[n1];
			Vector xr2 = Ut0[n3] - Ut0[n1];
			Vector xr3 = Ut0[n4] - Ut0[n1];
			Matrix Ar = Matrix(xr1, xr2, xr3); // Undeformed state
			Ar = tets[i].G.prod(Ar);
			
			Vector x1 = Ut[n2] - Ut[n1];
			Vector x2 = Ut[n3] - Ut[n1];
			Vector x3 = Ut[n4] - Ut[n1];

			Vector N1 = x3.cross(x1);
			Vector N2 = x2.cross(x3);
			Vector N3 = x1.cross(x2);
			Vector N4 = (x2 - x3).cross(x1 - x3);

			Matrix A = Matrix(x1, x2, x3); // Deformed state
			double vol = A.det()/6.0;
			Matrix F = A.prod(Ar.inv()); // Deformation gradient
			Matrix B = F.prod(F.trans());
			double J = F.det();
			double powJ23 = 1.0 + 2.0/3.0*(J - 1.0) - 1.0/9.0*(J-1.0)*(J-1.0);
			Matrix Ts = (B - I*B.trace()/3.0)*mu/(J*powJ23); // Deviatoric stress
			
			double J1 = Vn[n1]/Vn0[n1];
			double J2 = Vn[n2]/Vn0[n2];
			double J3 = Vn[n3]/Vn0[n3];
			double J4 = Vn[n4]/Vn0[n4];
			double Ja = (J1 + J2 + J3 + J4)/4.0;
			Matrix Tv = I*K*(Ja-1.0); // Isotropic stress
			
			Matrix T = Ts + Tv;
			
			double us = 0.5*mu*(B.trace()/powJ23 - 3.0);
			double uv = 0.5*K*( (J1-1.0)*(J1-1.0) + (J2-1.0)*(J2-1.0) + (J3-1.0)*(J3-1.0) + (J4-1.0)*(J4-1.0) )*0.25;
			Ue += (us + uv)*vol/J;
						
						
			if (F.det() < 0.0) {
					ninv++;
			}
				
			Ft[n1] += T.prod(N1 + N2 + N3)/6.0;
			Ft[n2] += T.prod(N1 + N3 + N4)/6.0;
			Ft[n3] += T.prod(N2 + N3 + N4)/6.0;
			Ft[n4] += T.prod(N1 + N2 + N4)/6.0;					
			
			
			// Update growth tensor G = g(y) I + (1-g(y)) n X n
			double gm = 1.0/(1.0 + exp(10.0*(0.25*(d2s[n1]+d2s[n2]+d2s[n3]+d2s[n4])/H - 1.0)))* 0.25*(gr[n1]+gr[n2]+gr[n3]+gr[n4]); 
			double vm = 1.0 - gm;
			
			
			Vector Ns = (N0[csn[n1]] + N0[csn[n2]] + N0[csn[n3]] + N0[csn[n4]]);
			Ns.normalize();
			
			// tangential 
			Matrix Lggt = (I - Matrix(Ns.x*Ns.x, Ns.x*Ns.y, Ns.x*Ns.z, Ns.x*Ns.y, Ns.y*Ns.y, Ns.y*Ns.z, Ns.x*Ns.z, Ns.y*Ns.z, Ns.z*Ns.z))*alfagt;
			// normal
			Matrix Lggn = Matrix(Ns.x*Ns.x, Ns.x*Ns.y, Ns.x*Ns.z, Ns.x*Ns.y, Ns.y*Ns.y, Ns.y*Ns.z, Ns.x*Ns.z, Ns.y*Ns.z, Ns.z*Ns.z)*alfagn; 
			
			Matrix Lgg = Lggt + Lggn; 
			Matrix Lgw = I*alfaw; 
			Matrix Lg = Lgg*gm + Lgw*vm; 
			tets[i].G += Lg.prod(tets[i].G)*dt;
		}
		
		// Midplane
		#pragma omp parallel for 
		for (int i = 0; i < nsn; i++) {
			int pt = SN[i];
			if ( Ut0[pt].y < mpy - 0.005 && Ut[pt].y > mpy - 0.005 ) {
				Ft[pt].y = 0;
			}
			if ( Ut0[pt].y > mpy + 0.005 && Ut[pt].y < mpy + 0.005 ) {
				Ft[pt].y = 0;
			}
		}	
						
		// Output
			Area = 0.0;
			for (int i = 0; i < nf; i++) {
				Vector N = (Ut[faces[i].n2]-Ut[faces[i].n1]).cross(Ut[faces[i].n3]-Ut[faces[i].n1]);
				Area += 0.5*N.length()*(gr[faces[i].n1] + gr[faces[i].n2] + gr[faces[i].n3])/3.0; // Only growing area
			}				
			double Volume = 0.0; // Volume
			for (int i = 0; i < nn; i++) Volume += Vn[i];
			cout  << setw(11) << step << setw(11) << t << setw(11) << Ue << setw(11) << Area << setw(11) << Volume << setw(11) << ninv << endl;
			dat << setw(13) << step << setw(13) << t << setw(13) << Ue << setw(13) << Area << setw(13) << Volume << setw(11) << ninv << endl;
		if (step%di == 0) {

//			writetxt(Ut, faces, SN, SNb, nsn, step);
			writeoff(Ut, faces, SN, SNb, nsn, step);

		}
		
		if (Ue > 100000 || Ue < -100000){// just to escape from error
			cout << "Error! Program terminated.\n";
			return 1;
		}
		#pragma omp parallel for 
		for (int i = 0; i < nn; i++) {

			Ft[i] -= Vt[i]*gamma*Vn0[i];  
			Vt[i] += Ft[i]/(Vn0[i]*rho)*dt; 
			Ut[i] += Vt[i]*dt;
			
			

			Ft[i].clear();
		}

		t += dt;
		step++;
	}		
	dat.close();
	cout << "Completed." << endl;
	cout  << setw(11) << step << setw(11) << t << setw(11) << Ue << setw(11) << Area << setw(11) << ninv << endl;
	return 0;
}

// Generates point-triangle proximity lists using the linked cell algorithm
void createNNLtriangle(vector<int>* NNLt, Vector* Ut, vector<Face>& faces, int* SN, int nsn, int nf, double hs, double mw) {
	int mx = max(1, (int)(5.0/mw)); 
	vector<int> head(mx*mx*mx, -1); 
	vector<int> list(nf);
	int xa, ya, za, xi, yi, zi;
	int pt, tri;
	Vector cog; 
	
	for (int i = 0; i < nf; i++) { // Divide triangles into cells
		Vector cog = (Ut[faces[i].n1] + Ut[faces[i].n2] + Ut[faces[i].n3])/3.0; // center of each triangle on the deformed surface
		int xa = (int)((cog.x + 2.5)/5.0*mx); 
		int ya = (int)((cog.y + 2.5)/5.0*mx); 
		int za = (int)((cog.z + 2.5)/5.0*mx); 
		if(mx*mx*za + mx*ya + xa < mx*mx*mx && mx*mx*za + mx*ya + xa >=0){ // like coordinates (xa, ya, za) 
		list[i]=head[mx*mx*za + mx*ya + xa];
		head[mx*mx*za + mx*ya + xa] = i;
		}
		else{
			cout << "i = " << i << ", mx*mx*za + mx*ya + xa = " << mx*mx*za + mx*ya + xa << endl;
			cout << xa << "  " << ya << "  " <<  za << "  " << mx << "  " <<endl;
			return;
		}

	}
	
	#pragma omp parallel for
	for (int i = 0; i < nsn; i++) { // Search cells around each point and build proximity list
		int pt = SN[i]; 
		NNLt[i].clear(); 
		int xa = (int)((Ut[pt].x + 2.5)/5.0*mx); 
		int ya = (int)((Ut[pt].y + 2.5)/5.0*mx);
		int za = (int)((Ut[pt].z + 2.5)/5.0*mx);
		for (int xi = max(0, xa-1); xi <= min(mx-1, xa+1); xi++)
		for (int yi = max(0, ya-1); yi <= min(mx-1, ya+1); yi++)
		for (int zi = max(0, za-1); zi <= min(mx-1, za+1); zi++) {
			int tri = head[mx*mx*zi + mx*yi + xi];
			while (tri != -1) {
				if ( pt != faces[tri].n1 && pt != faces[tri].n2 && pt != faces[tri].n3 ) {				
					if ( (closestPointTriangle(Ut[pt], Ut[faces[tri].n1], Ut[faces[tri].n2], Ut[faces[tri].n3], ub, vb, wb) - Ut[pt]).length() < hs) {
						NNLt[i].push_back(tri); // append tri no. to the list 
					}
				}
				tri = list[tri]; // recursion, find next
			}
		}
	}
}
Vector closestPointTriangle(Vector& p, Vector& a, Vector& b, Vector& c, double& u, double& v, double& w) {
	
	Vector ab = b - a;
	Vector ac = c - a;
	Vector ap = p - a;
	double d1 = ab.dot(ap);
	double d2 = ac.dot(ap);
	if (d1 <= 0.0 && d2 <= 0.0) {
		return a;
	}
	Vector bp = p - b;
	double d3 = ab.dot(bp);
	double d4 = ac.dot(bp);
	if (d3 >= 0.0 && d4 <= d3) {
		return b;
	}
	double vc = d1*d4 - d3*d2;
	if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
		return a + ab * v;
	}
	Vector cp = p - c;
	double d5 = ab.dot(cp);
	double d6 = ac.dot(cp);
	if (d6 >= 0.0 && d5 <= d6) {
		return c;
	}
	double vb = d5*d2 - d1*d6;
	if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
		return a + ac * w;
	}
	double va = d3*d6 - d5*d4;
	if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
		return b + (c - b) * w;
	}
	double denom = 1.0 / (va + vb + vc);
	return a + ab * v + ac * w;
}

void dist2surf(Vector* Ut, int* SN, int nn, int nsn, int* csn, double* d2s) {
	int p, j;
	#pragma omp parallel for private(p, j)
	for (int i = 0; i < nn; i++) {
		double d2min = 1e9;
		for (j = 0; j < nsn; j++) {
			double d2 = (Ut[SN[j]] - Ut[i]).dot(Ut[SN[j]] - Ut[i]);
			if (d2 < d2min) {
				d2min = d2;
				p = j;
			}
		}
		csn[i] = p;
		d2s[i] = sqrt(d2min);
	}
}

void writetxt(Vector* Ut, vector<Face>& faces, int* SN, int* SNb, int nsn, int step) {

	char povname[50];
	sprintf(povname, "B%d.txt", step);
	ofstream pov(povname);
	pov.setf(ios::fixed);
	pov.precision(5);

	// Normals
	Vector* N = new Vector[nsn];
	Vector Ntmp;
	for (int i = 0; i < faces.size(); i++) {
		Ntmp = (Ut[faces[i].n2] - Ut[faces[i].n1]).cross(Ut[faces[i].n3] - Ut[faces[i].n1]);
		N[SNb[faces[i].n1]] += Ntmp;
		N[SNb[faces[i].n2]] += Ntmp;
		N[SNb[faces[i].n3]] += Ntmp;
	}
	for (int i = 0; i < nsn; i++) N[i].normalize();
	
	pov << nsn << " \n";
	for (int i = 0; i < nsn; i++) {
		pov << Ut[SN[i]].x << " " << Ut[SN[i]].y << " " << Ut[SN[i]].z << " \n";
	}
	pov << nsn << " \n";
	for (int i = 0; i < nsn; i++) {
		pov << N[i].x << " " << N[i].y << " " << N[i].z << " \n";
	}
	pov << faces.size() << " \n";
	for (int i = 0; i < faces.size(); i++) {
		pov << SNb[faces[i].n1] << " " << SNb[faces[i].n2] << " " << SNb[faces[i].n3] << " \n"; 
	}
	pov.close();
	delete [] N;
}


void writeoff(Vector* Ut, vector<Face>& faces, int* SN, int* SNb, int nsn, int step) {

	char povname[50];
	sprintf(povname, "B%d.off", step);
	ofstream pov(povname);
	pov.setf(ios::fixed);
	pov.precision(5);

	pov << "OFF\n";
	pov << nsn << " " << faces.size() << " 0\n";
	for (int i = 0; i < nsn; i++) {
		pov << Ut[SN[i]].x << " " << Ut[SN[i]].y << " " << Ut[SN[i]].z << " \n";
	}
	for (int i = 0; i < faces.size(); i++) {
		pov << "3 "<< SNb[faces[i].n1] << " " << SNb[faces[i].n2] << " " << SNb[faces[i].n3] << " \n"; 
	}
	pov.close();
}
