#pragma once

#include <Basic/HeapObj.h>
#include <UHEMesh/HEMesh.h>
#include <UGM/UGM>

#include <Eigen/Eigen>
#include <Eigen/SVD>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
//#include <Eigen/SPQRSupport>
#include <Eigen/SparseCholesky>

typedef Eigen::SparseMatrix<float> SparseMatrixType;
typedef Eigen::Triplet<float> T;
typedef Eigen::SparseLU<SparseMatrixType> Solve_LU;
typedef Eigen::SimplicialLLT<SparseMatrixType> Solve_LLT;
typedef Eigen::SparseQR <SparseMatrixType, Eigen::COLAMDOrdering<int> >  Solve_QR;

namespace Ubpa {
	class TriMesh;
	class MinSurf;

	// mesh boundary == 1
	class Paramaterize : public HeapObj {
	public:
		Paramaterize();
		Paramaterize(Ptr<TriMesh> triMesh);
	public:
		static const Ptr<Paramaterize> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<Paramaterize>(triMesh);
		}
	public:
		void Clear();
		bool Init(Ptr<TriMesh> triMesh);

		bool Run();

		static bool EqualVector(vecf3 v1, vecf3 v2);
	protected:
		// kernel part of the algorithm
		virtual void Parameterize();
	protected:
		class V;
		class E;
		class P;
		class V : public TVertex<V, E, P> {
		public:
			vecf3 pos;
		};
		class E : public TEdge<V, E, P> { };
		class P :public TPolygon<V, E, P> { };
	protected:
		Ptr<TriMesh> triMesh;
		const Ptr<HEMesh<V>> heMesh; // vertice order is same with triMesh
	};
}
