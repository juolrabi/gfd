#ifndef _MESHDRAWER_HPP_INCLUDED_
#define _MESHDRAWER_HPP_INCLUDED_

#include "../Types/UintSet.hpp"
#include "Camera.hpp"
#include "../Mesh/Mesh.hpp"

namespace gfd
{

class MeshDrawer : public Camera
{
public:
	MeshDrawer() { clear(); }
	virtual ~MeshDrawer() {clear(); }

	void drawPrimalNodes(const Mesh &mesh, const Vector3 &col = Vector3(1,1,1), const UintSet &flag = UINTSETALL);

	void drawPrimalEdges(const Mesh &mesh, const Vector3 &col = Vector3(1,1,1), const UintSet &flag = UINTSETALL);
	void drawPrimalEdges(const Mesh& mesh, const Buffer<Vector3>& col, const UintSet& flag = UINTSETALL);
	void drawDualEdges(const Mesh &mesh, const Vector3 &col = Vector3(1,1,1), const UintSet &flag = UINTSETALL);

	void drawBoundaryFaces(const Mesh &mesh, const Vector3 &col = Vector3(1,1,1), const UintSet &flag = UINTSETALL);
	void drawBoundaryFaces(const Mesh &mesh, const Buffer<Vector3> &col, const UintSet &flag = UINTSETALL);
	void drawBoundaryFacesByNodeColor(const Mesh &mesh, const Buffer<Vector3> &col, const UintSet &flag = UINTSETALL);

protected:

};

}

#endif //_MESHDRAWER_HPP_INCLUDED_
