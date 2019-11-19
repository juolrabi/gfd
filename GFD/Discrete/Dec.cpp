#include "Dec.hpp"

using namespace gfd;

Sparse<sign> &Dec::integrateDerivative(const FormGrade grade, Sparse<sign> &d) const
{
	uint i, j;
	if(grade == fg_prim0 || grade == fg_dual1) {
		Buffer< Buffer< pair<uint,sign> > > inc(m_mesh.getEdgeLocals());
		for(i=0; i<inc.size(); i++) {
            const Buffer<uint> &ele = m_mesh.getEdgeNodes(i);
            inc[i].resize(ele.size());
            for(j=0; j<ele.size(); j++) inc[i][j] = pair<uint,sign>(ele[j], m_mesh.getEdgeIncidence(i, ele[j]));
        }
		d.setFull(m_mesh.getNodeLocals(), inc, m_mesh.getExternalNodes());
		if(grade == fg_dual1) d.setTranspose(d);
        else d.setTranspose(d.setTranspose(d));
		return d; 
	}
	if(grade == fg_prim1 || grade == fg_dual2) {
		Buffer< Buffer< pair<uint,sign> > > inc(m_mesh.getFaceLocals());
		for(i=0; i<inc.size(); i++) {
            const Buffer<uint> &ele = m_mesh.getFaceEdges(i);
            inc[i].resize(ele.size());
            for(j=0; j<ele.size(); j++) inc[i][j] = pair<uint,sign>(ele[j], m_mesh.getFaceIncidence(i, ele[j]));
        }
		d.setFull(m_mesh.getEdgeLocals(), inc, m_mesh.getExternalEdges());
		if(grade == fg_dual2) d.setTranspose(d);
		return d;
	}
	if(grade == fg_prim2 || grade == fg_dual3) {
		Buffer< Buffer< pair<uint,sign> > > inc(m_mesh.getBodyLocals());
		for(i=0; i<inc.size(); i++) {
            const Buffer<uint> &ele = m_mesh.getBodyFaces(i);
            inc[i].resize(ele.size());
            for(j=0; j<ele.size(); j++) inc[i][j] = pair<uint,sign>(ele[j], m_mesh.getBodyIncidence(i, ele[j]));
        }
		d.setFull(m_mesh.getFaceLocals(), inc, m_mesh.getExternalFaces());
		if(grade == fg_dual3) d.setTranspose(d);
		return d;
	}
	if(grade == fg_prim3 || grade == fg_dual4)
	{
		Buffer< Buffer< pair<uint,sign> > > inc(m_mesh.getQuadLocals());
		for(i=0; i<inc.size(); i++) {
            const Buffer<uint> &ele = m_mesh.getQuadBodies(i);
            inc[i].resize(ele.size());
            for(j=0; j<ele.size(); j++) inc[i][j] = pair<uint,sign>(ele[j], m_mesh.getQuadIncidence(i, ele[j]));
        }
		d.setFull(m_mesh.getBodyLocals(), inc, m_mesh.getExternalBodies());
		if(grade == fg_dual4) d.setTranspose(d);
		return d;
	}
	cout << "BlockMesh::integrate(Derivative &) -> Cannot integrate derivative due to unsuitable grade." << endl;
	return d;
}

Buffer< pair<uint,uint> > Dec::getMyExternals(const Buffer< pair<uint,uint> > &ext) const {
	uint i, j, k;
    const uint rank = getMPIrank();
    const uint ranks = getMPIranks();
	for(j=0; j<ranks; j++) {
        if(j == rank) continue;
		for(i=0,k=0; i<ext.size(); i++) {
			if(ext[i].first == j) k++;
		}
		sendMPI(&k, sizeof(uint), j, 0);
		if(k == 0) continue;
		Buffer<uint> link(k);
		for(i=0,k=0; i<ext.size(); i++) {
			if(ext[i].first == j) link[k++] = ext[i].second;
		}
		sendMPI(&link[0], k * sizeof(uint), j, 1);
	}
	Buffer< pair<uint,uint> > mext;
	for(j=0; j<ranks; j++) {
        if(j == rank) continue;
		recvMPI(&k, sizeof(uint), j, 0);
        if(k == 0) continue;
        const uint size = mext.size();
        mext.resize(size + k);
		Buffer<uint> link(k);
		recvMPI(&link[0], k * sizeof(uint), j, 1);
		for(i=0; i<k; i++) mext[size + i] = pair<uint,uint>(j, link[i]);
	}
	return mext;
}

