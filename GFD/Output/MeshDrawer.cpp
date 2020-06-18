#include "MeshDrawer.hpp"
#include <fstream>
#include <iostream>

using namespace gfd;

void MeshDrawer::drawPrimalNodes(const Mesh &mesh, const Vector3 &col, const UintSet &flag)
{
	uint i;
//	const double h = 0.1;
	for(i=0; i<mesh.getNodeSize(); i++)
	{
		if(!flag.includes(mesh.getNodeFlag(i))) continue;
		const Vector4 p = mesh.getNodePosition(i);
		drawPoint(p, col);
/*		drawLine(p-Vector4(0,h,0,0), p+Vector4(0,h,0,0), col);
		drawLine(p-Vector4(0,0,h,0), p+Vector4(0,0,h,0), col);
		drawLine(p-Vector4(0,0,0,h), p+Vector4(0,0,0,h), col);
*/	}
}

void MeshDrawer::drawPrimalEdges(const Mesh &mesh, const Vector3 &col, const UintSet &flag)
{
	uint i;
	for(i=0; i<mesh.getEdgeSize(); i++)
	{
		if(!flag.includes(mesh.getEdgeFlag(i))) continue;
		const Buffer<uint> &en = mesh.getEdgeNodes(i);
		const Vector4 p0 = mesh.getNodePosition(en[0]);
		const Vector4 p1 = mesh.getNodePosition(en[1]);
		drawLine(p0, p1, col);
	}
}

void MeshDrawer::drawPrimalEdges(const Mesh& mesh, const Buffer<Vector3>& col, const UintSet& flag)
{
	uint i;
	for (i = 0; i < mesh.getEdgeSize(); i++)
	{
		if (!flag.includes(mesh.getEdgeFlag(i))) continue;
		const Buffer<uint>& en = mesh.getEdgeNodes(i);
		const Vector4 p0 = mesh.getNodePosition(en[0]);
		const Vector4 p1 = mesh.getNodePosition(en[1]);
		drawLine(p0, p1, col[i]);
	}
}

void MeshDrawer::drawDualEdges(const Mesh &mesh, const Vector3 &col, const UintSet &flag)
{
	uint i;
	if(mesh.getFaceSize() == 0)
	{
		for(i=0; i<mesh.getNodeSize(); i++)
		{
			if(!flag.includes(mesh.getNodeFlag(i))) continue;
			const Buffer<uint> &dn = mesh.getNodeEdges(i);
			if(dn.size() == 1) drawLine(mesh.getEdgePosition(dn[0]), mesh.getNodePosition(i), col);
			else if(dn.size() == 2) drawLine(mesh.getEdgePosition(dn[0]), mesh.getEdgePosition(dn[1]), col);
		}
		return;
	}
	if(mesh.getBodySize() == 0)
	{
		for(i=0; i<mesh.getEdgeSize(); i++)
		{
			if(!flag.includes(mesh.getEdgeFlag(i))) continue;
			const Buffer<uint> &dn = mesh.getEdgeFaces(i);
			if(dn.size() == 1) drawLine(mesh.getFacePosition(dn[0]), mesh.getEdgePosition(i), col);
			else if(dn.size() == 2) drawLine(mesh.getFacePosition(dn[0]), mesh.getFacePosition(dn[1]), col);
		}
		return;
	}
	if(mesh.getQuadSize() == 0)
	{
		for(i=0; i<mesh.getFaceSize(); i++)
		{
			if(!flag.includes(mesh.getFaceFlag(i))) continue;
			const Buffer<uint> &dn = mesh.getFaceBodies(i);
			if(dn.size() == 1) drawLine(mesh.getBodyPosition(dn[0]), mesh.getFacePosition(i), col);
			else if(dn.size() == 2) drawLine(mesh.getBodyPosition(dn[0]), mesh.getBodyPosition(dn[1]), col);
		}
		return;
	}
	for(i=0; i<mesh.getBodySize(); i++)
	{
		if(!flag.includes(mesh.getBodyFlag(i))) continue;
		const Buffer<uint> &dn = mesh.getBodyQuads(i);
		if(dn.size() == 1) drawLine(mesh.getQuadPosition(dn[0]), mesh.getBodyPosition(i), col);
		else if(dn.size() == 2) drawLine(mesh.getQuadPosition(dn[0]), mesh.getQuadPosition(dn[1]), col);
	}
}

void MeshDrawer::drawBoundaryFaces(const Mesh &mesh, const Vector3 &col, const UintSet &flag)
{
	uint i, j;
	uint fs = 0;
	Buffer<uint> f(mesh.getFaceSize());
	Buffer<double> sq(m_isvg ? mesh.getFaceSize() : 0);
	for(i=0; i<mesh.getFaceSize(); i++)
	{
		if(!flag.includes(mesh.getFaceFlag(i))) continue;
		const Buffer<uint> &b = mesh.getFaceBodies(i);
		if(b.size() >= 2) continue;
//std::cout << "a " << mesh.getFaceNodes(i).size();
//const Buffer<uint> &e = mesh.getFaceEdges(i);
//for(j=0; j<e.size(); j++) std::cout << " " << e[j] << " {" << mesh.getEdgeNodes(e[j])[0] << " " << mesh.getEdgeNodes(e[j])[1] << "}";
//std::cout << std::endl;
		const Vector4 fp = mesh.getFaceAverage(i);
//std::cout << "b" << std::endl;
		const Vector4 vcam = getPosition() - fp;
		if(b.size() == 1 && vcam.dot(mesh.getFaceDeviation(i, mesh.getBodyAverage(b[0]))) >= 0.0) continue;
		if(m_isvg) sq[fs] = vcam.lensq();
		f[fs++] = i;
	}

	// order faces
	if(m_isvg)
	{
		for(i=1; i<fs; i++)
		{
			const uint fi = f[i];
			const double sqi = sq[i];
			for(j=i; j>0; j--)
			{
				if(sq[j-1] > sqi) break;
				sq[j] = sq[j-1];
				f[j] = f[j-1];
			}
			sq[j] = sqi;
			f[j] = fi;
		}
	}

	for(i=0; i<fs; i++)
	{
		const Buffer<uint> n = mesh.getFaceNodes(f[i]);
		Buffer<Vector4> np(n.size());
		Buffer<Vector3> nc(n.size());

		for(j=0; j<n.size(); j++)
		{
			np[j] = mesh.getNodePosition(n[j]);
			nc[j] = col; // sininen
		}
		drawPolygon(np, nc);
	}
}

void MeshDrawer::drawBoundaryFaces(const Mesh &mesh, const Buffer<Vector3> &col, const UintSet &flag)
{
	uint i, j;
	uint fs = 0;
	Buffer<uint> f(mesh.getFaceSize());
	Buffer<double> sq(m_isvg ? mesh.getFaceSize() : 0);
	for(i=0; i<mesh.getFaceSize(); i++)
	{
		if(!flag.includes(mesh.getFaceFlag(i))) continue;
		const Buffer<uint> &b = mesh.getFaceBodies(i);
		if(b.size() >= 2) continue;
		const Vector4 fp = mesh.getFaceAverage(i);
		const Vector4 vcam = getPosition() - fp;
		if(b.size() == 1 && vcam.dot(mesh.getFaceDeviation(i, mesh.getBodyAverage(b[0]))) >= 0.0) continue;
		if(m_isvg) sq[fs] = vcam.lensq();
		f[fs++] = i;
	}


	// order faces
	if(m_isvg)
	{
		for(i=1; i<fs; i++)
		{
			const uint fi = f[i];
			const double sqi = sq[i];
			for(j=i; j>0; j--)
			{
				if(sq[j-1] > sqi) break;
				sq[j] = sq[j-1];
				f[j] = f[j-1];
			}
			sq[j] = sqi;
			f[j] = fi;
		}
	}

	for(i=0; i<fs; i++)
	{
		const Buffer<uint> n = mesh.getFaceNodes(f[i]);
		Buffer<Vector4> np(n.size());
		Buffer<Vector3> nc(n.size());

		for(j=0; j<n.size(); j++)
		{
			np[j] = mesh.getNodePosition(n[j]);
			nc[j] = col[f[i]];
		}
		drawPolygon(np, nc);
	}
}

void MeshDrawer::drawBoundaryFacesByNodeColor(const Mesh &mesh, const Buffer<Vector3> &col, const UintSet &flag)
{
	uint i, j;
	uint fs = 0;
	Buffer<uint> f(mesh.getFaceSize());
	Buffer<double> sq(m_isvg ? mesh.getFaceSize() : 0);
	for(i=0; i<mesh.getFaceSize(); i++)
	{
		if(!flag.includes(mesh.getFaceFlag(i))) continue;
		const Buffer<uint> &b = mesh.getFaceBodies(i);
		if(b.size() >= 2) continue;
		const Vector4 fp = mesh.getFaceAverage(i);
		const Vector4 vcam = getPosition() - fp;
		if(b.size() == 1 && vcam.dot(mesh.getFaceDeviation(i, mesh.getBodyAverage(b[0]))) >= 0.0) continue;
		if(m_isvg) sq[fs] = vcam.lensq();
		f[fs++] = i;
	}


	// order faces
	if(m_isvg)
	{
		for(i=1; i<fs; i++)
		{
			const uint fi = f[i];
			const double sqi = sq[i];
			for(j=i; j>0; j--)
			{
				if(sq[j-1] > sqi) break;
				sq[j] = sq[j-1];
				f[j] = f[j-1];
			}
			sq[j] = sqi;
			f[j] = fi;
		}
	}

	for(i=0; i<fs; i++)
	{
		const Buffer<uint> n = mesh.getFaceNodes(f[i]);
		Buffer<Vector4> np(n.size());
		Buffer<Vector3> nc(n.size());

		for(j=0; j<n.size(); j++)
		{
			np[j] = mesh.getNodePosition(n[j]);
			nc[j] = col[n[j]];
		}
		drawPolygon(np, nc);
	}
}

