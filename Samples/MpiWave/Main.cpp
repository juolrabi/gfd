#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Mesh/PartMesh.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include "../../GFD/Types/Types.hpp"
#include "../../GFD/Types/Random.hpp"
#include "../../GFD/Types/Buffer.hpp"
#include "../../GFD/Types/MpiEasy.hpp"
#include "../../GFD/Discrete/Dec.hpp"
#include <iostream>

using namespace std;
using namespace gfd;

void get2Form(const Vector4 &p, double *result) {
	const double psq = 20.0 * (p + Vector4(0.2,0.2,0.0,0.0)).lensq();
	result[0] = (psq < PI ? 1.0 + cos(psq) : 0.0);
}

int main() {
	initMPI();
	uint i, j;

	// create or load mesh
	BuilderMesh bmesh(2);
	bmesh.createSnubSquareGrid(Vector2(-1,-1), Vector2(1,1), 0.2);
	PartMesh mesh(getMPIrank(), getMPIranks(), 2);
	mesh.createPartFromFlags(bmesh);

	// initialize wave
	Dec dec(mesh);
	Form<double> f2;
	dec.integrateForm(fg_prim2, get2Form, 2, f2);

	cout << "rank " << getMPIrank() << " " << f2.m_height << endl;
/*
	Form<double> f2(mesh.getFaceLocals(), 0.0);
	for(i=0; i<f2.m_height; i++) {
		const double psq = 20.0 * (mesh.getFaceAverage2(i) + Vector2(0.2,0.2)).lensq();
		f2.m_val[i] = (psq < PI ? 1.0 + cos(psq) : 0.0) * mesh.getFaceVector2(i).xy;
	}
*/
	Form<double> f1(mesh.getEdgeLocals(), 0.0);

	Buffer< Buffer< pair<uint,sign> > > buf(f2.m_height);
	for(i=0; i<f2.m_height; i++) {
		const Buffer<uint> &ele = mesh.getFaceEdges(i);
		buf[i].resize(ele.size());
		for(j=0; j<ele.size(); j++) buf[i][j] = pair<uint,sign>(ele[j], mesh.getFaceIncidence(i, ele[j]));
	}
	Derivative d1;
	d1.setFull(f1.m_height, buf, mesh.getExternalEdges());
	buf.clear();
	Hodge<double> h2(f2.m_height, 0.0);
	for(i=0; i<f2.m_height; i++) h2.m_val[i] = mesh.getFaceHodge(i);
	Hodge<double> h1(f1.m_height, 0.0);
	for(i=0; i<f1.m_height; i++) h1.m_val[i] = 1.0 / mesh.getEdgeHodge(i);
	Sparse<double> delta1;
	delta1.setTimes(delta1.setTimes(h1, delta1.setTranspose(d1)), h2).scaleLeft(-0.0005);

	for(uint iter=0; iter<=200; iter++) {
		f2 += d1 * f1;
		f1 += delta1 * f2;

		if((iter % 10) != 0) continue;

		// draw the mesh
		Buffer<Vector3> col(mesh.getFaceSize(), Vector3(0,0,0));
		for(i=0; i<f2.m_height; i++) col[i] = f2.m_val[i] / mesh.getFaceVector2(i).xy * Vector3(1,1,1);

		Picture pic(500,500);
		MeshDrawer drawer;
		const Vector3 vo(0,0,0);
		const Vector3 vp(0,0,10);
		const Vector3 vx = TwoVector3(Vector3(0,1,0), vp-vo).dual().unit() / 2.5;
		const Vector3 vy = TwoVector3(vp-vo, vx).dual().unit() / 2.5;
		drawer.initPosition(Vector4(vp,0), Vector4(vo,0), Vector4(vx,0), Vector4(vy,0));
		drawer.initPicture(&pic);
		drawer.drawBoundaryFaces(mesh, col);
		Text name;
		name << "mesh" << iter/10 << "_" << getMPIrank() << ".bmp";
		pic.save(name.str(), true);

	}

/*
	// initialize forms and integrate initial solution
	Form<double> f0;
	bm.integrateForm(fg_prim0, get0Form, 2, f0);

	// initialize derivatives and integrate Hodge operators
	Derivative d0;
*/
/*

	Buffer< Buffer< pair<uint,Sign> > > inc(mesh.getEdgeLocals());
	for(uint k=0; k<mesh.getEdgeLocals(); k++) {
		const Buffer<uint> n = mesh.getEdgeNodes(k);
		for(uint i=0; i<n.size(); i++) inc[k].push_back(pair<uint,Sign>(n[i], mesh.getEdgeIncidence(k, n[i])));
	}
	const uint nlocs = mesh.getNodeLocals();
	Buffer< pair<uint,uint> > ext(mesh.getNodeSize() - nlocs);
	for(uint i=0; i<ext.size(); i++) ext[i] = pair<uint,uint>(mesh.getNodePart(nlocs+i), mesh.getNodeLink(nlocs+i));
	//d0.setFull(nlocs, inc, ext);

*/
/*	bm.integrateDerivative(fg_prim0, d0);

	Derivative d0T;
	d0T.setTranspose(d0);
	Hodge<double> h1;
	bm.integrateHodge(fg_prim1, get1Hodge, 2, h1);
	Hodge<double> dh0;
	bm.integrateHodge(fg_dual0, get0Hodge, 2, dh0);
	Sparse<double> delta0 = dh0 * d0T * h1;
	Form<double> f1;
	f1.setFullOfZeros(d0.m_height);

	uint i, j;
	Buffer<uint> c0(f0.m_height, 0);
	Buffer<uint> c1(f1.m_height, 0);
	uint steps0 = 1;
{
	ullong sum = 0;
	uint num = 0;
	Sparse<double> test;
	test.setTimes(d0, delta0);
	for(i=0; i<test.m_height; i++) {
		c1[i] = uint(sqrt(0.5 * test.getValue(i, i)) + 1.0);
		sum += c1[i];
		num++;
	}
	test.setTimes(delta0, d0);
	for(i=0; i<test.m_height; i++) {
		c0[i] = uint(sqrt(0.5 * test.getValue(i, i)) + 1.0);
		sum += c0[i];
		num++;
	}
	sumMPI(&sum, 1);
	sumMPI(&num, 1);
	cout << "average " << double(sum) / double(num) << endl;
	steps0 = uint(sum / num);

	for(i=0; i<c0.size(); i++) {
		uint steps = steps0;
		for(j=0; steps<c0[i]; j++) steps *= 2;
		c0[i] = j;
	}
	for(i=0; i<c1.size(); i++) {
		uint steps = steps0;
		for(j=0; steps<c1[i]; j++) steps *= 2;
		c1[i] = j;
	}
//	c0[i] /= uint(sqrt(mindiag));
//	for(i=0; i<c1.size(); i++) c1[i] /= uint(sqrt(mindiag));
//	delta0.scaleRight(-2.0 / mindiag);
}


	SparseSplit<double> sd0(c1, d0, c0);
	SparseSplit<double> sdelta0(delta0, sd0.getCategories());
	cout << sdelta0.m_term.size() << " " << sd0.m_term.size() << endl;
	uint maxstep = 1;
	for(i=0; i<sdelta0.m_term.size() || i<sd0.m_term.size(); i++) {
		if(i<sd0.m_term.size()) sd0.m_term[i].scaleRight(1.0 / double(steps0 * maxstep));
		maxstep *= 2;
		if(i<sdelta0.m_term.size()) sdelta0.m_term[i].scaleRight(-1.0 / double(steps0 * maxstep));
	}
	delta0.scaleRight(-1.0 / double(steps0 * maxstep));
*/
/*	const uint maxpow = 2;
	uint maxstep = 1;
	for(i=1; i<maxpow; i++) maxstep *= 3;

	Random rnd(5345345 * getMPIrank());
	Buffer<uint> c0(f0.m_height, 0);
	for(i=0; i<c0.size(); i++) c0[i] = rnd.getUint() % maxpow;
	Buffer<uint> c1(f1.m_height, 0);
	for(i=0; i<c1.size(); i++) c1[i] = rnd.getUint() % maxpow;
	SparseSplit<double> sd0(c1, d0);
	SparseSplit<double> sdelta0(c0, delta0);

	double fac = 1.0;
	for(i=1; i<sd0.m_term.size(); i++) {
		fac /= 3.0;
		sd0.m_term[i].scaleRight(fac);
	}
	fac = 1.0;
	for(i=1; i<sdelta0.m_term.size(); i++) {
		fac /= 3.0;
		sdelta0.m_term[i].scaleRight(fac);
	}
*/
/*
	// iterate
	for(uint k=0; k<=10*steps0; k++)
	{
		for(i=0; i<maxstep; i++) {
			uint step = maxstep;
			for(j=0; j<sd0.m_term.size(); j++) {
				if((i % step) == (step / 2)) {
					f1 += sd0.m_term[j] * f0;
					if(k == 0) cout << "odd " << j << endl;
				}
				step /= 2;
			}
			step = maxstep;
			for(j=0; j<sdelta0.m_term.size(); j++) {
				step /= 2;
				if((i % step) == (step / 2)) {
					f0 += sdelta0.m_term[j] * f1;
					if(k == 0) cout << "even " << j + 1 << endl;
				}
			}

		if((k % steps0) != 0) continue;

		Picture pic0(500,500);
		Buffer<double> val;
		for(j=0; j<pic0.getHeight(); j++)
		{
			for(i=0; i<pic0.getWidth(); i++)
			{
				const Vector4 p(0.006 * i - 1.0, 0.006 * j - 1.0,0,0);
				bm.interpolateForm(fg_prim0, f0, p, val);
				sumMPI(&val[0], val.size());
				const Vector4 col((val.size() > 0 ? val[0] : 0.0), (val.size() > 1 ? val[1] : 0.0), (val.size() > 2 ? val[2] : 0.0), 1.0);
				pic0.setColor(i, j, col);
			}
		}
		if(getMPIrank() == 0)
		{
			Text path;
			path << "kuva" << k / steps0 << ".bmp" << endl;
			pic0.save(path.str(), true);
		}
	}


	// draw f0 at the end
	double check = 0.0;
	Picture pic0(500,500);
	Buffer<double> val(3, 0.0);
	for(uint j=0; j<pic0.getHeight(); j++) {
		for(uint i=0; i<pic0.getWidth(); i++) {
			const Vector4 p(0.006 * i - 1.0, 0.006 * j - 1.0,0,0);
			bm.interpolateForm(fg_prim0, f0, p, val);
			sumMPI(&val[0], val.size());
			const Vector4 col((val.size() > 0 ? val[0] : 0.0), (val.size() > 1 ? val[1] : 0.0), (val.size() > 2 ? val[2] : 0.0), 1.0);
			check += col.toVector3().len();
			pic0.setColor(i, j, col);
		}
	}
	if(getMPIrank() == 0) {
		cout << "check = " << check << endl;
		pic0.save("kuva.bmp", true);
	}
*/

	finalizeMPI();

    return 0;
}
