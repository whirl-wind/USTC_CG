#include <Engine/MeshEdit/Paramaterize.h>

#include <Engine/MeshEdit/MinSurf.h>

#include <Engine/Primitive/TriMesh.h>

using namespace Ubpa;

using namespace std;
using namespace Eigen;
#define My_Pi 3.1415927
Ubpa::Paramaterize::Paramaterize()
{
	heMesh->Clear();
	triMesh = nullptr;
}

Paramaterize::Paramaterize(Ptr<TriMesh> triMesh) 
	: heMesh(make_shared<HEMesh<V>>()) {
	// TODO
	//cout << "Paramaterize::Paramaterize:" << endl
	//	<< "\t" << "not implemented" << endl;
	Init(triMesh);
}

void Paramaterize::Clear() {
	// TODO
	//cout << "Paramaterize::Clear:" << endl
	//	<< "\t" << "not implemented" << endl;
	heMesh->Clear();
	triMesh = nullptr;
}

bool Paramaterize::Init(Ptr<TriMesh> triMesh) {
	// TODO
	//cout << "Paramaterize::Init:" << endl
	//	<< "\t" << "not implemented" << endl;
	Clear();

	if (triMesh == nullptr)
		return true;

	if (triMesh->GetType() == TriMesh::INVALID) {
		printf("ERROR::MinSurf::Init:\n"
			"\t""trimesh is invalid\n");
		return false;
	}

	// init half-edge structure
	size_t nV = triMesh->GetPositions().size();
	vector<vector<size_t>> triangles;
	triangles.reserve(triMesh->GetTriangles().size());
	for (auto triangle : triMesh->GetTriangles())
		triangles.push_back({ triangle->idx[0], triangle->idx[1], triangle->idx[2] });
	heMesh->Reserve(nV);
	heMesh->Init(triangles);

	if (!heMesh->IsTriMesh() || !heMesh->HaveBoundary()) {
		printf("ERROR::MinSurf::Init:\n"
			"\t""trimesh is not a triangle mesh or hasn't a boundaries\n");
		heMesh->Clear();
		return false;
	}

	// triangle mesh's positions ->  half-edge structure's positions
	for (int i = 0; i < nV; i++) {
		auto v = heMesh->Vertices().at(i);
		v->pos = triMesh->GetPositions()[i].cast_to<vecf3>();
	}

	this->triMesh = triMesh;
	return true;
}

bool Paramaterize::Run() {
	// TODO
	//cout << "Paramaterize::Init:" << endl
	//	<< "\t" << "not implemented" << endl;
	if (heMesh->IsEmpty() || !triMesh) {
		printf("ERROR::MinSurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	Parameterize();

	// half-edge structure -> triangle mesh
	size_t nV = heMesh->NumVertices();
	size_t nF = heMesh->NumPolygons();
	vector<pointf3> positions;
	vector<unsigned> indice;
	positions.reserve(nV);
	indice.reserve(3 * nF);
	for (auto v : heMesh->Vertices())
		positions.push_back(v->pos.cast_to<pointf3>());
	for (auto f : heMesh->Polygons()) { // f is triangle
		for (auto v : f->BoundaryVertice()) // vertices of the triangle
			indice.push_back(static_cast<unsigned>(heMesh->Index(v)));
	}

	triMesh->Init(indice, positions);

	return true;
}

void Ubpa::Paramaterize::Parameterize()
{
	int count_b = 0;
	int row, col;
	row = col = heMesh->NumVertices();
	VectorXf bx(row);
	VectorXf by(row);
	bx.setZero();
	by.setZero();
	float R = 1.0;
	SparseMatrixType L;
	L.resize(row, col);
	vector<T> tripletList_L;
	for (auto v : heMesh->Vertices()) {
		int i = heMesh->Index(v);
		if (v->IsBoundary()) {
			count_b++;
			//cout << count_b << endl;
			tripletList_L.push_back(T(i, i, 1));
			continue;
		}
		int count = 0;
		for (auto vv : v->AdjVertices()) {
			int j = heMesh->Index(vv);
			tripletList_L.push_back(T(i, j, -1));
			count++;
		}
		tripletList_L.push_back(T(i, i, count));
	}
	L.setFromTriplets(tripletList_L.begin(), tripletList_L.end());
	L.makeCompressed();
	for (auto v : heMesh->Vertices()) {
		if (v->IsBoundary()) {
			int i = 0;

			auto he = v->HalfEdge();
			while(!he->IsBoundary()) he = he->RotateNext();
			auto he0 = he;
			int idx = heMesh->Index(he->Origin());
			//if (he->Origin()->IsBoundary()) cout << count_b - i << endl;;
			float xx = R* cos(-2 * My_Pi / (count_b * 1.0) * i);
			float yy = R* sin(-2 * My_Pi / (count_b * 1.0) * i);
			bx(idx) = xx;
			by(idx) = yy;
			he = he->Next();
			i++;
			while (he != he0) {
				idx = heMesh->Index(he->Origin());
				//if (he->Origin()->IsBoundary()) cout << count_b - i << endl;;
				xx = R* cos(-2 * My_Pi / (count_b * 1.0) * i);
				yy = R* sin(-2 * My_Pi / (count_b * 1.0) * i);
				bx(idx) = xx;
				by(idx) = yy;
				he = he->Next();
				i++;
			}
			break;
		}
	}
	
	VectorXf x(row);
	VectorXf y(row);
	Solve_LU solver;
	solver.compute(L);
	x = solver.solve(bx);
	y = solver.solve(by);
	/*
	cout << "L\n"<< L << endl;
	cout << "bx\n" << bx.transpose() << endl;
	cout << "by\n" << by.transpose() << endl;
	cout << "x\n" << x.transpose() << endl;
	cout << "y\n" << y.transpose() << endl;
	*/
	for (auto v : heMesh->Vertices()) {
		int i = heMesh->Index(v);
		v->pos = vecf3(x(i), y(i), 0);
	}
}

bool Paramaterize::EqualVector(vecf3 v1, vecf3 v2)
{
	if (abs(v1[0] - v2[0]) < 1e-6 && abs(v1[2] - v2[2]) < 1e-6 && abs(v1[1] - v2[1]) < 1e-6)
		return true;
	else
		return false;
}
