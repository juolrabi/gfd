#include "Dec.hpp"

using namespace gfd;

Dec::Dec(const PartMesh &mesh) 
: m_mesh(mesh) {
	m_lowdim = 0;
	if(m_mesh.getEdgeSize() == 0) m_highdim = 0;
	else if(m_mesh.getFaceSize() == 0) m_highdim = 1;
	else if(m_mesh.getBodySize() == 0) m_highdim = 2;
	else if(m_mesh.getQuadSize() == 0) m_highdim = 3;
	else m_highdim = 4;
	maxMPI(&m_highdim, 1);
 }

Dec::Dec(const PartMesh &mesh, const uint lowdim, const uint highdim) 
: m_mesh(mesh) {
	setLowDimension(lowdim);
	setHighDimension(highdim);
}

void Dec::setLowDimension(const uint lowdim) { 
	m_lowdim = lowdim; 
	if(m_lowdim > 4) m_lowdim = 4;
}
void Dec::setHighDimension(const uint highdim) { 
	m_highdim = highdim; 
	if(m_highdim > 4) m_highdim = 4;
}

void Dec::combineExternals(const Buffer< pair<uint,uint> > &ext, Buffer<Quadrature> &q) const {
	Quadrature *extq = &q[q.size() - ext.size()];
	const Buffer< pair<uint,uint> > mext = getMyExternals(ext);
	const uint rank = getMPIrank();
	for(uint i=0; i<ext.size(); i++) {
		Quadrature &qi = extq[i];
		if(ext[i].first == rank) q[ext[i].second].combine(qi);
		else {
			uint size = qi.size();
			sendMPI(&size, sizeof(uint), ext[i].first, 0);
			if(size == 0) continue;
			sendMPI(&qi[0], size * sizeof(double), ext[i].first, 1);
		}
		qi.clear();
	}
	for(uint i=0; i<mext.size(); i++) {
		uint size;
		recvMPI(&size, sizeof(uint), mext[i].first, 0);
		if(size == 0) continue;
		Buffer<double> buf(size);
		recvMPI(&buf[0], size * sizeof(double), mext[i].first, 1);
		q[mext[i].second].combine(buf);
	}
}

Diagonal<bool> &Dec::integrateFlags(const FormGrade grade, const UintSet &flag, Diagonal<bool> &result) const {
	switch(FormGradeDimension(grade)) {
	case 0: {
		initResult(m_mesh.getNodeLocals(), result);
		for(uint i=0; i<result.m_val.size(); i++) result.m_val[i] = flag.includes(m_mesh.getNodeFlag(i));
		return result;
	}
	case 1: {
		initResult(m_mesh.getEdgeLocals(), result);
		for(uint i=0; i<result.m_val.size(); i++) result.m_val[i] = flag.includes(m_mesh.getEdgeFlag(i));
		return result;
	}
	case 2: {
		initResult(m_mesh.getFaceLocals(), result);
		for(uint i=0; i<result.m_val.size(); i++) result.m_val[i] = flag.includes(m_mesh.getFaceFlag(i));
		return result;
	}
	case 3: {
		initResult(m_mesh.getBodyLocals(), result);
		for(uint i=0; i<result.m_val.size(); i++) result.m_val[i] = flag.includes(m_mesh.getBodyFlag(i));
		return result;
	}
	default: {
		initResult(m_mesh.getQuadLocals(), result);
		for(uint i=0; i<result.m_val.size(); i++) result.m_val[i] = flag.includes(m_mesh.getQuadFlag(i));
		return result;
	}
	}
}

Sparse<sign> &Dec::integrateDerivative(const FormGrade grade, Sparse<sign> &result) const
{
	uint i, j;
	switch(grade) {
	case fg_prim0:
	case fg_dual1: {
		Buffer< Buffer< pair<uint,sign> > > inc(m_mesh.getEdgeLocals());
		for(i=0; i<inc.size(); i++) {
            const Buffer<uint> &ele = m_mesh.getEdgeNodes(i);
            inc[i].resize(ele.size());
            for(j=0; j<ele.size(); j++) inc[i][j] = pair<uint,sign>(ele[j], m_mesh.getEdgeIncidence(i, ele[j]));
        }
		result.setFull(m_mesh.getNodeLocals(), inc, m_mesh.getExternalNodes());
		if(grade == fg_dual1) result.setTranspose(result);
		return result; 
	}
	case fg_prim1:
	case fg_dual2: {
		Buffer< Buffer< pair<uint,sign> > > inc(m_mesh.getFaceLocals());
		for(i=0; i<inc.size(); i++) {
            const Buffer<uint> &ele = m_mesh.getFaceEdges(i);
            inc[i].resize(ele.size());
            for(j=0; j<ele.size(); j++) inc[i][j] = pair<uint,sign>(ele[j], m_mesh.getFaceIncidence(i, ele[j]));
        }
		result.setFull(m_mesh.getEdgeLocals(), inc, m_mesh.getExternalEdges());
		if(grade == fg_dual2) result.setTranspose(result);
		return result;
	}
	case fg_prim2:
	case fg_dual3: {
		Buffer< Buffer< pair<uint,sign> > > inc(m_mesh.getBodyLocals());
		for(i=0; i<inc.size(); i++) {
            const Buffer<uint> &ele = m_mesh.getBodyFaces(i);
            inc[i].resize(ele.size());
            for(j=0; j<ele.size(); j++) inc[i][j] = pair<uint,sign>(ele[j], m_mesh.getBodyIncidence(i, ele[j]));
        }
		result.setFull(m_mesh.getFaceLocals(), inc, m_mesh.getExternalFaces());
		if(grade == fg_dual3) result.setTranspose(result);
		return result;
	}
	case fg_prim3:
	case fg_dual4: {
		Buffer< Buffer< pair<uint,sign> > > inc(m_mesh.getQuadLocals());
		for(i=0; i<inc.size(); i++) {
            const Buffer<uint> &ele = m_mesh.getQuadBodies(i);
            inc[i].resize(ele.size());
            for(j=0; j<ele.size(); j++) inc[i][j] = pair<uint,sign>(ele[j], m_mesh.getQuadIncidence(i, ele[j]));
        }
		result.setFull(m_mesh.getBodyLocals(), inc, m_mesh.getExternalBodies());
		if(grade == fg_dual4) result.setTranspose(result);
		return result;
	}
	default: {
		cout << "BlockMesh::integrateDerivative -> Cannot integrate derivative due to unsuitable grade." << endl;
		return result;
	}
	}
}

Sparse<sign> &Dec::integrateDerivative(const FormGrade grade, const UintSet &flag, Sparse<sign> &result) const {
	uint i, j;
	switch(grade) {
	case fg_prim0: {
		Buffer< Buffer< pair<uint,sign> > > inc(m_mesh.getEdgeLocals());
		for(i=0; i<m_mesh.getNodeSize(); i++) {
			if(!flag.includes(m_mesh.getNodeFlag(i))) continue;
            const Buffer<uint> &ele = m_mesh.getNodeEdges(i);
            for(j=0; j<ele.size(); j++) {
				if(ele[j] >= inc.size()) continue;
				inc[ele[j]].push_back(pair<uint,sign>(i, m_mesh.getEdgeIncidence(ele[j], i)));
			}
        }
		return result.setFull(m_mesh.getNodeLocals(), inc, m_mesh.getExternalNodes());
	}
	case fg_dual1: {
		Buffer< Buffer< pair<uint,sign> > > inc(m_mesh.getEdgeLocals());
		for(i=0; i<inc.size(); i++) {
			if(!flag.includes(m_mesh.getEdgeFlag(i))) continue;
            const Buffer<uint> &ele = m_mesh.getEdgeNodes(i);
            inc[i].resize(ele.size());
            for(j=0; j<ele.size(); j++) {
				inc[i][j] = pair<uint,sign>(ele[j], m_mesh.getEdgeIncidence(i, ele[j]));
			}
        }
		result.setFull(m_mesh.getNodeLocals(), inc, m_mesh.getExternalNodes());
		return result.setTranspose(result);
	}
	case fg_prim1: {
		Buffer< Buffer< pair<uint,sign> > > inc(m_mesh.getFaceLocals());
		for(i=0; i<m_mesh.getEdgeSize(); i++) {
			if(!flag.includes(m_mesh.getEdgeFlag(i))) continue;
            const Buffer<uint> &ele = m_mesh.getEdgeFaces(i);
            for(j=0; j<ele.size(); j++) {
				if(ele[j] >= inc.size()) continue;
				inc[ele[j]].push_back(pair<uint,sign>(i, m_mesh.getFaceIncidence(ele[j], i)));
			}
        }
		return result.setFull(m_mesh.getEdgeLocals(), inc, m_mesh.getExternalEdges());
	}
	case fg_dual2: {
		Buffer< Buffer< pair<uint,sign> > > inc(m_mesh.getFaceLocals());
		for(i=0; i<inc.size(); i++) {
			if(!flag.includes(m_mesh.getFaceFlag(i))) continue;
            const Buffer<uint> &ele = m_mesh.getFaceEdges(i);
            inc[i].resize(ele.size());
            for(j=0; j<ele.size(); j++) {
				inc[i][j] = pair<uint,sign>(ele[j], m_mesh.getFaceIncidence(i, ele[j]));
			}
        }
		result.setFull(m_mesh.getEdgeLocals(), inc, m_mesh.getExternalEdges());
		return result.setTranspose(result);
	}
	case fg_prim2: {
		Buffer< Buffer< pair<uint,sign> > > inc(m_mesh.getBodyLocals());
		for(i=0; i<m_mesh.getFaceSize(); i++) {
			if(!flag.includes(m_mesh.getFaceFlag(i))) continue;
            const Buffer<uint> &ele = m_mesh.getFaceBodies(i);
            for(j=0; j<ele.size(); j++) {
				if(ele[j] >= inc.size()) continue;
				inc[ele[j]].push_back(pair<uint,sign>(i, m_mesh.getBodyIncidence(ele[j], i)));
			}
        }
		return result.setFull(m_mesh.getFaceLocals(), inc, m_mesh.getExternalFaces());
	}
	case fg_dual3: {
		Buffer< Buffer< pair<uint,sign> > > inc(m_mesh.getBodyLocals());
		for(i=0; i<inc.size(); i++) {
			if(!flag.includes(m_mesh.getBodyFlag(i))) continue;
            const Buffer<uint> &ele = m_mesh.getBodyFaces(i);
            inc[i].resize(ele.size());
            for(j=0; j<ele.size(); j++) {
				inc[i][j] = pair<uint,sign>(ele[j], m_mesh.getBodyIncidence(i, ele[j]));
			}
        }
		result.setFull(m_mesh.getFaceLocals(), inc, m_mesh.getExternalFaces());
		return result.setTranspose(result);
	}
	case fg_prim3: {
		Buffer< Buffer< pair<uint,sign> > > inc(m_mesh.getQuadLocals());
		for(i=0; i<m_mesh.getBodySize(); i++) {
			if(!flag.includes(m_mesh.getBodyFlag(i))) continue;
            const Buffer<uint> &ele = m_mesh.getBodyQuads(i);
            for(j=0; j<ele.size(); j++) {
				if(ele[j] >= inc.size()) continue;
				inc[ele[j]].push_back(pair<uint,sign>(i, m_mesh.getQuadIncidence(ele[j], i)));
			}
        }
		return result.setFull(m_mesh.getBodyLocals(), inc, m_mesh.getExternalBodies());
	}
	case fg_dual4: {
		Buffer< Buffer< pair<uint,sign> > > inc(m_mesh.getQuadLocals());
		for(i=0; i<inc.size(); i++) {
			if(!flag.includes(m_mesh.getQuadFlag(i))) continue;
            const Buffer<uint> &ele = m_mesh.getQuadBodies(i);
            inc[i].resize(ele.size());
            for(j=0; j<ele.size(); j++) {
				inc[i][j] = pair<uint,sign>(ele[j], m_mesh.getQuadIncidence(i, ele[j]));
			}
        }
		result.setFull(m_mesh.getBodyLocals(), inc, m_mesh.getExternalBodies());
		return result.setTranspose(result);
	}
	default: {
		cout << "BlockMesh::integrateDerivative -> Cannot integrate derivative due to unsuitable grade." << endl;
		return result;
	}
	}
}

Sparse<double> &Dec::integrateCurvatureDerivative(SymMatrix4 curv(const Vector4 &), const FormGrade grade, const UintSet &flag, Sparse<double> &result) const {
	uint i, j;
	switch(grade) {
	case fg_prim0: {
		Buffer< Buffer< pair<uint,double> > > inc(m_mesh.getEdgeLocals());
		for(i=0; i<m_mesh.getNodeSize(); i++) {
			if(!flag.includes(m_mesh.getNodeFlag(i))) continue;
			const Vector4 p = m_mesh.getNodePosition(i);
            const Buffer<uint> &ele = m_mesh.getNodeEdges(i);
            for(j=0; j<ele.size(); j++) {
				if(ele[j] >= inc.size()) continue;
				const Vector4 d = m_mesh.getEdgePosition(ele[j]) - p;
				const double c = 1.0 + d.dot(curv(p + 0.5 * d) * d);
				inc[ele[j]].push_back(pair<uint,double>(i, m_mesh.getEdgeIncidence(ele[j], i) * c));
			}
        }
		return result.setFull(m_mesh.getNodeLocals(), inc, m_mesh.getExternalNodes());
	}
	case fg_dual1: {
		Buffer< Buffer< pair<uint,double> > > inc(m_mesh.getEdgeLocals());
		for(i=0; i<inc.size(); i++) {
			if(!flag.includes(m_mesh.getEdgeFlag(i))) continue;
			const Vector4 p = m_mesh.getEdgePosition(i);
            const Buffer<uint> &ele = m_mesh.getEdgeNodes(i);
            inc[i].resize(ele.size());
            for(j=0; j<ele.size(); j++) {
				const Vector4 d = m_mesh.getNodePosition(ele[j]) - p;
				const double c = 1.0 + d.dot(curv(p + 0.5 * d) * d);
				inc[i][j] = pair<uint,double>(ele[j], m_mesh.getEdgeIncidence(i, ele[j]) * c);
			}
        }
		result.setFull(m_mesh.getNodeLocals(), inc, m_mesh.getExternalNodes());
		return result.setTranspose(result);
	}
	case fg_prim1: {
		Buffer< Buffer< pair<uint,double> > > inc(m_mesh.getFaceLocals());
		for(i=0; i<m_mesh.getEdgeSize(); i++) {
			if(!flag.includes(m_mesh.getEdgeFlag(i))) continue;
			const Vector4 p = m_mesh.getEdgePosition(i);
            const Buffer<uint> &ele = m_mesh.getEdgeFaces(i);
            for(j=0; j<ele.size(); j++) {
				if(ele[j] >= inc.size()) continue;
				const Vector4 d = m_mesh.getFacePosition(ele[j]) - p;
				const double c = 1.0 + d.dot(curv(p + 0.5 * d) * d);
				inc[ele[j]].push_back(pair<uint,double>(i, m_mesh.getFaceIncidence(ele[j], i) * c));
			}
        }
		return result.setFull(m_mesh.getEdgeLocals(), inc, m_mesh.getExternalEdges());
	}
	case fg_dual2: {
		Buffer< Buffer< pair<uint,double> > > inc(m_mesh.getFaceLocals());
		for(i=0; i<inc.size(); i++) {
			if(!flag.includes(m_mesh.getFaceFlag(i))) continue;
			const Vector4 p = m_mesh.getFacePosition(i);
            const Buffer<uint> &ele = m_mesh.getFaceEdges(i);
            inc[i].resize(ele.size());
            for(j=0; j<ele.size(); j++) {
				const Vector4 d = m_mesh.getEdgePosition(ele[j]) - p;
				const double c = 1.0 + d.dot(curv(p + 0.5 * d) * d);
				inc[i][j] = pair<uint,double>(ele[j], m_mesh.getFaceIncidence(i, ele[j]) * c);
			}
        }
		result.setFull(m_mesh.getEdgeLocals(), inc, m_mesh.getExternalEdges());
		return result.setTranspose(result);
	}
	case fg_prim2: {
		Buffer< Buffer< pair<uint,double> > > inc(m_mesh.getBodyLocals());
		for(i=0; i<m_mesh.getFaceSize(); i++) {
			if(!flag.includes(m_mesh.getFaceFlag(i))) continue;
			const Vector4 p = m_mesh.getFacePosition(i);
            const Buffer<uint> &ele = m_mesh.getFaceBodies(i);
            for(j=0; j<ele.size(); j++) {
				if(ele[j] >= inc.size()) continue;
				const Vector4 d = m_mesh.getBodyPosition(ele[j]) - p;
				const double c = 1.0 + d.dot(curv(p + 0.5 * d) * d);
				inc[ele[j]].push_back(pair<uint,double>(i, m_mesh.getBodyIncidence(ele[j], i) * c));
			}
        }
		return result.setFull(m_mesh.getFaceLocals(), inc, m_mesh.getExternalFaces());
	}
	case fg_dual3: {
		Buffer< Buffer< pair<uint,double> > > inc(m_mesh.getBodyLocals());
		for(i=0; i<inc.size(); i++) {
			if(!flag.includes(m_mesh.getBodyFlag(i))) continue;
			const Vector4 p = m_mesh.getBodyPosition(i);
            const Buffer<uint> &ele = m_mesh.getBodyFaces(i);
            inc[i].resize(ele.size());
            for(j=0; j<ele.size(); j++) {
				const Vector4 d = m_mesh.getFacePosition(ele[j]) - p;
				const double c = 1.0 + d.dot(curv(p + 0.5 * d) * d);
				inc[i][j] = pair<uint,double>(ele[j], m_mesh.getBodyIncidence(i, ele[j]) * c);
			}
        }
		result.setFull(m_mesh.getFaceLocals(), inc, m_mesh.getExternalFaces());
		return result.setTranspose(result);
	}
	case fg_prim3: {
		Buffer< Buffer< pair<uint,double> > > inc(m_mesh.getQuadLocals());
		for(i=0; i<m_mesh.getBodySize(); i++) {
			if(!flag.includes(m_mesh.getBodyFlag(i))) continue;
			const Vector4 p = m_mesh.getBodyPosition(i);
            const Buffer<uint> &ele = m_mesh.getBodyQuads(i);
            for(j=0; j<ele.size(); j++) {
				if(ele[j] >= inc.size()) continue;
				const Vector4 d = m_mesh.getQuadPosition(ele[j]) - p;
				const double c = 1.0 + d.dot(curv(p + 0.5 * d) * d);
				inc[ele[j]].push_back(pair<uint,double>(i, m_mesh.getQuadIncidence(ele[j], i) * c));
			}
        }
		return result.setFull(m_mesh.getBodyLocals(), inc, m_mesh.getExternalBodies());
	}
	case fg_dual4: {
		Buffer< Buffer< pair<uint,double> > > inc(m_mesh.getQuadLocals());
		for(i=0; i<inc.size(); i++) {
			if(!flag.includes(m_mesh.getQuadFlag(i))) continue;
			const Vector4 p = m_mesh.getQuadPosition(i);
            const Buffer<uint> &ele = m_mesh.getQuadBodies(i);
            inc[i].resize(ele.size());
            for(j=0; j<ele.size(); j++) {
				const Vector4 d = m_mesh.getBodyPosition(ele[j]) - p;
				const double c = 1.0 + d.dot(curv(p + 0.5 * d) * d);
				inc[i][j] = pair<uint,double>(ele[j], m_mesh.getQuadIncidence(i, ele[j]) * c);
			}
        }
		result.setFull(m_mesh.getBodyLocals(), inc, m_mesh.getExternalBodies());
		return result.setTranspose(result);
	}
	default: {
		cout << "BlockMesh::integrateCurvatureDerivative -> Cannot integrate derivative due to unsuitable grade." << endl;
		return result;
	}
	}
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

