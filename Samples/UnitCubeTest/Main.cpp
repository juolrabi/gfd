/**
 * A program for testing convergence of solution of u'' + curl curl u = j.
 * First order problem stands as
 *    u' = -curl v + j,
 *    v' = curl u.
 * 
 * Analytic solution is given as 
 *    u = t² (x-x²) (y-y²) [y, -x]^T,
 *    v = t³ xy (5x + 5y - 4 - 6xy) / 3,
 *    j = t a + t³ b,
 *    a = 2 (x-x²) (y-y²) [y  -x]^T,
 *    b = [x (5x + 10y - 12xy - 4), 
		  -y (5y + 10x - 12xy - 4)]^T / 3.
 * 
 * Author: Jukka Räbinä, University of Jyväskylä, 2020.
 */

#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Mesh/PartMesh.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include "../../GFD/Types/Types.hpp"
#include "../../GFD/Types/MpiEasy.hpp"
#include "../../GFD/Discrete/Dec.hpp"
#include <iostream>
#include <chrono>

using namespace std;
using namespace gfd;

// modify these constants:
const uint space_steps = 8; // number of space elements per unit length
const uint time_steps = 12; // number of time elements per unit time
const uint grid_type = 0; // 2d grid type: 0 = squares, 1 = triangles, 2 = snubsquare

double funcU(const Buffer<double> &q) { // analytic function for u
	const Vector2 v(q[0], q[1]);
	const Vector2 p(q[2], q[3]);
	const double common = p.x * (1.0 - p.x) * p.y * (1.0 - p.y);
	return common * (p.y * v.x - p.x * v.y);
}
double funcV(const Buffer<double> &q) { // analytic function for v
	const double v(q[0]);
	const Vector2 p(q[1], q[2]);
	return p.x * p.y * (5.0 * p.x + 5.0 * p.y - 6.0 * p.x * p.y - 4.0) * v / 3.0;
}
double funcA(const Buffer<double> &q) { // analytic function for a
	const Vector2 v(q[0], q[1]);
	const Vector2 p(q[2], q[3]);
	const double common = 2.0 * p.x * (1.0 - p.x) * p.y * (1.0 - p.y);
	return common * (p.y * v.x - p.x * v.y);
}
double funcB(const Buffer<double> &q) { // analytic function for b
	const Vector2 v(q[0], q[1]);
	const Vector2 p(q[2], q[3]);
	return (p.x * (5.0 * p.x + 10.0 * p.y - 12.0 * p.x * p.y - 4.0) * v.x 
		- p.y * (5.0 * p.y + 10.0 * p.x - 12.0 * p.x * p.y - 4.0) * v.y) / 3.0;
}
bool createMesh(PartMesh &mesh) {
	if(mesh.getNumberOfParts() > 1) {
		if(mesh.getPart() == 0) cout << "MPI mesh generation is not implemented." << endl;
		return false;
	}
	BuilderMesh bmesh(mesh.getDimension());
	const double hx = 1.0 / double(space_steps);
	if(grid_type == 0) bmesh.createGrid(Vector4(0,0,0,0), Vector4(1,1,0,0), hx);
	else if(grid_type == 1) bmesh.createTriangleGrid(Vector2(0,0), Vector2(1,1), hx);
	else if(grid_type == 2) bmesh.createSnubSquareGrid(Vector2(0,0), Vector2(1,1), hx);
	bmesh.fillBoundaryFlags(1);
	
	if(mesh.getPart() == 0) mesh.swap(bmesh);
	return true;
}
void drawSolution(const PartMesh &mesh, const Column<double> &u, const Vector2 &minp, const Vector2 &maxp, const double scale, const string &path) {
	MeshDrawer drawer;
	const Vector3 vo(0.5 * (minp + maxp), 0.0);
	const Vector3 vp(vo + Vector3(0, -0, 10.0 * (maxp - minp).len()));
	const Vector3 vx = TwoVector3(Vector3(0, 1, 0), vp - vo).dual().unit() / (maxp.x - minp.x);
	const Vector3 vy = TwoVector3(vp - vo, vx).dual().unit() / (maxp.x - minp.x);
	drawer.initPosition(Vector4(vp, 0), Vector4(vo, 0), Vector4(vx, 0), Vector4(vy, 0));
	Picture pic(1000, 1000 * (maxp.y - minp.y) / (maxp.x - minp.x));
	drawer.initPicture(&pic);

	Buffer<Vector3> col(mesh.getFaceLocals(), Vector3(0, 0, 0));
	Buffer<double> full_u;
	Dec dec(mesh, 0, mesh.getDimension());
	dec.getFullBuffer(fg_prim1, u, full_u);
	for (uint i=0; i<col.size(); i++) {
		const Buffer<uint> &e = mesh.getFaceEdges(i);
		SymMatrix2 A(0,0,0);
		Vector2 b(0,0);
		for(uint j=0; j<e.size(); j++) {
			const Vector2 vi = mesh.getEdgeVector(e[j]).toVector2();
			A += vi.outerProduct();
			b += vi * full_u[e[j]];
		}
		b = A.inverse() * b;
		col[i] = scale * Vector3(b.x, b.y, 0.0);
	}
	drawer.drawBoundaryFaces(mesh, col);
	for(uint j=0; j<pic.getHeight(); j++) {
		for(uint i=0; i<pic.getWidth(); i++) {
			const Vector4 color = pic.getColor(i, j);
			double c[3] = {color.x, color.y, color.z};
			sumMPI(c, 3);
			pic.setColor(i, j, Vector4(c[0], c[1], c[2], 1));
		}
	}
	if(getMPIrank() == 0) pic.save(path);
}
int main() {
	initMPI();
	auto starttime = chrono::system_clock::now();

	// create mesh
	PartMesh mesh(0, 1, 2);
	if(!createMesh(mesh)) return 0;
	if(getMPIrank() == 0) {
		cout << "Mesh generated in " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;
		starttime = chrono::system_clock::now();
	}

	// initialize matrices and vectors
	Dec dec(mesh, 0,  mesh.getDimension());
	Derivative d1;
	dec.integrateDerivative(fg_prim1, d1);
	Diagonal<double> h2;
	dec.integrateHodge(HodgeUnit1, 0, fg_prim2, h2);
	Diagonal<double> h1i;
	dec.integrateHodge(HodgeUnit2, 0, fg_dual1, UintSet(0), h1i);
	Sparse<double> delta1;
	const double ht = 1.0 / double(time_steps);
	delta1.setScale(-ht * ht, h1i * transpose(d1) * h2);
	Column<double> u(0.0);
	dec.integrateZeroForm(fg_prim1, u);
	Column<double> v(0.0);
	dec.integrateZeroForm(fg_prim2, v);
	Column<double> a(0.0);
	dec.integrateForm(funcA, 20, fg_prim1, a);
	Column<double> b(0.0);
	dec.integrateForm(funcB, 20, fg_prim1, UintSet(0), b);

	// integrate forward in time
	for(uint i=0; i<time_steps; i++) {
		const double time0 = i * ht;
		const double time1 = (i + 1) * ht;
		const double sq0 = time0 * time0;
		const double sq1 = time1 * time1;
		const double integral_t = 0.5 * (sq1 - sq0); 
		const double integral_t3 = 0.25 * (sq1 * sq1 - sq0 * sq0); 
		v += d1 * u;
		u += delta1 * v + scaled(integral_t, a) + scaled(integral_t3, b);
	}
	v += scaled(0.5, d1 * u); // even solution for v
	v.scale(ht); // rescale v to match funcV
	if(getMPIrank() == 0) {
		cout << "Time-integration processed in " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;
		starttime = chrono::system_clock::now();
	}

	// compare to analytic solution
	// here we compare discrete forms
	// if you want to compare continuous fields, choose interpolation strategy and then compare interpolated fields
	Column<double> u0(0.0);
	dec.integrateForm(funcU, 20, fg_prim1, u0);
	Column<double> v0(0.0);
	dec.integrateForm(funcV, 20, fg_prim2, v0);
	Diagonal<double> h1;
	h1i.trimSparse();
	h1.setInverse(h1i);
	const double u_error_norm = sqrt((u - u0).getDot(h1 * (u - u0)));
	const double u0_norm = sqrt(u0.getDot(h1 * u0));
	const double v_error_norm = sqrt((v - v0).getDot(h2 * (v - v0)));
	const double v0_norm = sqrt(v0.getDot(h2 * v0));
	if(getMPIrank() == 0) {
		cout << "   u_error_norm = " << u_error_norm << endl;
		cout << "   u0_norm = " << u0_norm << endl;
		cout << "   u_relative_norm = " << u_error_norm / u0_norm << endl;
		cout << "   v_error_norm = " << v_error_norm << endl;
		cout << "   v0_norm = " << v0_norm << endl;
		cout << "   v_relative_norm = " << v_error_norm / v0_norm << endl;
	}
	if(getMPIrank() == 0) {
		cout << "Error norms printed in " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;
		starttime = chrono::system_clock::now();
	}

	// draw solution u and analytic solution u0
	drawSolution(mesh, u, Vector2(0,0), Vector2(1,1), 20.0, "kuva.bmp");
	drawSolution(mesh, u0, Vector2(0,0), Vector2(1,1), 20.0, "kuva0.bmp");
	if(getMPIrank() == 0) cout << "Solution pictures saved in " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;
	finalizeMPI();

    return 0;
}
