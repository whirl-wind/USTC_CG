#pragma once

#include <Basic/HeapObj.h>
//#include <Engine/Primitive/MassSpring.h>
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
	class Simulate : public HeapObj {
	public:
		Simulate(const std::vector<pointf3>& plist,
			const std::vector<unsigned>& elist) {
			edgelist = elist;
			this->positions.resize(plist.size());
			for (int i = 0; i < plist.size(); i++)
			{
				for (int j = 0; j < 3; j++)
				{
					this->positions[i][j] = plist[i][j];
				}
			}
		};
	public:
		static const Ptr<Simulate> New(const std::vector<pointf3>& plist,
			const std::vector<unsigned> &elist) {
			return Ubpa::New<Simulate>(plist, elist);
		}
	public:
		// clear cache data
		void Clear();

		// init cache data (eg. half-edge structure) for Run()
		bool Init();
		//bool Init();

		// call it after Init()
		bool Run();
		
		const std::vector<pointf3>& GetPositions() const { return positions; };

		const float GetStiff() { return stiff; };
		void SetStiff(float k) { stiff = k; Init();};
		const float GetTimeStep() { return h; };
		void SetTimeStep(float k) { h = k; Init();};
		std::vector<unsigned>& GetFix() { return this->fixed_id; };
		void SetFix(const std::vector<unsigned>& f) { this->fixed_id = f; Init();};
		const std::vector<pointf3>& GetVelocity() { return velocity; };
		//void SetVelocity(const std::vector<pointf3>& v) { velocity = v; };

		void SetLeftFix();
		void CheckBlow();
		void CheckBounce();
		void CheckFix();

	private:
		void Blow();
		void Bounce();
		// kernel part of the algorithm
		void SimulateOnce();
		void Fix(bool bv);
		void Optimize();
		void UpdateParameter();

	private:
		float h = 0.03f;  //步长
		float stiff;
		std::vector<unsigned> fixed_id;  //fixed point id

		//mesh data
		std::vector<unsigned> edgelist;

		//simulation data
		std::vector<pointf3> positions;
		std::vector<pointf3> velocity;
		std::vector<float> origin_length; //原长
		std::vector<float> k; //劲度系数
		std::vector<float> mass; //质量
		std::vector<size_t> fix;
		bool is_blow_;
		bool is_bounce_;
		bool is_fix_;

		//Caculate Matrix
		Solve_LU solver;
		SparseMatrixType K;
		SparseMatrixType M;
		SparseMatrixType M_inv;
		SparseMatrixType L;
		SparseMatrixType J;
		Eigen::VectorXf x;
		Eigen::VectorXf xf;
		Eigen::VectorXf y;
		Eigen::VectorXf v;
		Eigen::VectorXf d;
		Eigen::VectorXf b;
		Eigen::VectorXf f_int;
		Eigen::VectorXf f_ext;
	};
}
