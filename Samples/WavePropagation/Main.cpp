#include "../../GFD/BlockDec/BlockMesh.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include "../../GFD/Types/Types.hpp"
#include "../../GFD/Types/Random.hpp"
#include "../../GFD/Types/Buffer.hpp"
#include "../../GFD/Types/MpiEasy.hpp"
#include "../../GFD/Discrete/Split.hpp"
#include "../../GFD/Discrete/Dec.hpp"
#include <iostream>

using namespace std;
using namespace gfd;

void timesScalar(double *result, const double *right)
{
	result[0] *= right[0];
}
void plusScalarTimesScalar(double *result, const double *left, const double *right)
{
	result[0] += left[0] * right[0];
}

void get0Form(const Vector4 &p, double *result)
{
	const double plen = 3.0 * p.len();
	if(plen < 1.0) result[0] = 3.0 * (1.0 + cos(PI * plen));
	else result[0] = 0.0;
}

void get1(const Vector4 &p, double *result) {
	result[0] = cos(3*p.lensq());
}
void get2(const Vector4 &p, double *result) {
	result[0] = cos(3*p.lensq());
	result[1] = sin(1*p.lensq());
}
void get3(const Vector4 &p, double *result) {
	result[0] = cos(3*p.lensq());
	result[1] = sin(1*p.lensq());
	result[2] = p.lensq();
}

void get1Form(const Vector4 &p, double *result)
{
/*	result[0] = 1.0;
	result[1] = 0.0;
	result[2] = 0.0;
*/
//	result[0] = cos(3*p.lensq());
	result[0] = -0.6 * p.x * sin(3*p.lensq());
	result[1] = -0.6 * p.y * sin(3*p.lensq());
	result[2] = -0.6 * p.z * sin(3*p.lensq());


/*	result[0] = result[1] = result[2] = 0.0;
	if(p.lensq() < 1.0) result[0] = 1.0;
	if((Vector4(1,1,0,0) - p).lensq() < 1.0) result[2] = 1.0;
//	result[1] = sin(10.0 * p.x);
*/}

void get0Hodge(const Vector4 &p, double *result)
{
//	if((Vector4(1,1,0,0) - p).lensq() < 1.0) result[0] = 0.5;
//	else result[0] = 1.0;
	result[0] = 1.0;
}
void get1Hodge(const Vector4 &p, double *result)
{
	result[0] = 1.0;
	result[1] = 0.0;
	result[2] = 1.0;
	//result[3] = 0.0;
	//result[4] = 0.0;
	//result[5] = 1.0;
}

void get0Hodge2d(const Vector4 &p, double *result)
{
	result[0] = 1.0;
//	if((Vector4(0.5,0.5,0,0) - p).lensq() < 2.0) result[0] = 1.0;
//	else result[0] = 0.0;
}
void get1Hodge2d(const Vector4 &p, double *result)
{
/*	result[0] = 1.0;
	result[1] = 0.0;
	result[2] = 1.0;
*/	if((Vector4(0.5,0.5,0,0) - p).lensq() < 2.0) {
		result[0] = 1.0;
		result[1] = 0.0;
		result[2] = 1.0;
	}
	else {
		result[0] = 0.0;
		result[1] = 0.0;
		result[2] = 0.0;
	}
}


int main()
{
	initMPI();

	// 2d
	const uint w = 21;
	const double h = (w>1 ? 3.0 / double(w - 1) : 3.0);
	Buffer<uchar> data(w * w, 3);
//	data[4] = 11;
	for(uint i=0; i<w; i++)
	{
		const uint ii = (i-10) * (i-10);
		for(uint j=0; j<w; j++)
		{
			const uint jj = (j-10) * (j-10);
			if(ii + jj < 7) data[j * w + i] = 11;
			else if(ii + jj < 25) data[j * w + i] = 4;
		}
	}
	BlockMesh bm;
	bm.init(Vector2(-1,-1), Vector2(h,h), w, w, &data[0]);

	// create a PartMesh from BlockMesh
	PartMesh mesh(getMPIrank(), getMPIranks(), 2);
	bm.toMesh(mesh);

{ 
	FormGrade grade = fg_dual1;
	FormGrade hgrade = FormGradeDual(grade);
	const uint num = 2;
	Form<double> f0;
	bm.integrateForm(grade, get2, num, f0);
	Hodge<double> h0;
	bm.integrateHodge(grade, get1Hodge, num, h0);
	Form<double> f1;
	f1.setTimes(h0, f0);

	//Form<double> decf1;
	//bm.integrateForm(hgrade, get1, num, decf1);

	Dec dec(mesh);
	Form<double> decf0;
	dec.integrateForm(grade, get2, num, decf0);
	Hodge<double> dech0;
	dec.integrateHodge(grade, get1Hodge, num, dech0);
	Form<double> decf1;
	decf1.setTimes(dech0, decf0);

	double sum = 0.0;
	double summ = 0.0;
	for(uint j=0; j<f1.m_height; j++) {
		const double val = f1.getValue(j);
		const double dif = val - decf1.getValue(j);
		sum += dif * dif;
		summ += val * val;
		//cout << mesh.getNodeHodge(j) << " " << dech0.m_val[j] << endl;
	}
	sumMPI(&sum, 1);
	sumMPI(&summ, 1);

	// draw f0 at the end
	double check = 0.0;
	Picture pic0(500,500);
	Buffer<double> val(3, 0.0);
	for(uint j=0; j<pic0.getHeight(); j++) {
		for(uint i=0; i<pic0.getWidth(); i++) {
			const Vector4 p(0.006 * i - 1.0, 0.006 * j - 1.0,0,0);
			bm.interpolateForm(hgrade, f1, p, val);
			sumMPI(&val[0], val.size());
			const Vector4 col((val.size() > 0 ? val[0] : 0.0), (val.size() > 1 ? val[1] : 0.0), (val.size() > 2 ? val[2] : 0.0), 1.0);
			check += col.toVector3().len();
			pic0.setColor(i, j, col);
		}
	}
	if(getMPIrank() == 0) {
		cout << "sum = " << sum << endl;
		cout << "summ = " << summ << endl;
		cout << "check = " << check << endl;
		pic0.save("kuva.bmp", true);
	}
}

	finalizeMPI();

    return 0;


	// draw the mesh
	MeshDrawer drawer;
	const Vector3 vo(0,0,0);
	const Vector3 vp(0,0,10);
	const Vector3 vx = TwoVector3(Vector3(0,1,0), vp-vo).dual().unit() / 4.5;
	const Vector3 vy = TwoVector3(vp-vo, vx).dual().unit() / 4.5;
	drawer.initPosition(Vector4(vp,0), Vector4(vo,0), Vector4(vx,0), Vector4(vy,0));
	drawer.initSvg(500, 500);
	drawer.drawBoundaryFaces(mesh, Vector3(0,0.5,1));
	Text name;
	name << "mesh" << getMPIrank() << ".svg";
	drawer.saveSvg(name.str());

	Text text;
	mesh.writeStatistics(text);
	text.save("stat.txt");

	// initialize forms and integrate initial solution
	Form<double> f0;
	bm.integrateForm(fg_prim0, get0Form, 2, f0);

	// initialize derivatives and integrate Hodge operators
	Derivative d0;

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
	bm.integrateDerivative(fg_prim0, d0);
/*	{
	Form<double> fff = d0 * f0;
	double sum = 0.0;
	double maxdiag = 0.0;
	for(uint i=0; i<fff.m_height; i++) {
		const double diag = fff.getValue(i);
		if(diag > maxdiag) maxdiag = diag;
		sum += diag;
	}
	maxMPI(&maxdiag, 1);
	sumMPI(&sum, 1);
	cout << maxdiag << " is max and sum is " << sum << endl;
	}
*/
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


/*	Random rnd(5345345 * getMPIrank());
	Buffer<uint> c0(f0.m_height, 0);
	for(i=0; i<c0.size(); i++) c0[i] = rnd.getUint() % 2;
	Buffer<uint> c1(f1.m_height, 0);
	for(i=0; i<c1.size(); i++) c1[i] = rnd.getUint() % 2;
*/

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

//			f0 += delta0 * f1;
//			if(k == 0) cout << "even" << endl;
		}

/*		for(i=0; i<2*maxstep; i++) {
			uint step = 2*maxstep;
			for(j=0; j<maxpow; j++) {
				if((i % step) == 0) {
					if(j < sdelta0.m_term.size()) f0 += sdelta0.m_term[j] * f1;
					if(k == 0) cout << "even " << j << endl;
				}
				else if(((i + maxstep) % step) == 0) {
					if(j < sd0.m_term.size()) f1 += sd0.m_term[j] * f0;
					if(k == 0) cout << "odd " << j << endl;
				}
				step /= 3;
			}
		}
*/

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


	finalizeMPI();

    return 0;
}
