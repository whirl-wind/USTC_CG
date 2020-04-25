#include <Engine/MeshEdit/Simulate.h>


#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;

#define My_g 9.8

void Simulate::Clear() {
	this->positions.clear();
	this->velocity.clear();
}

//周期性风
int Tt = 50;
int t = 0;
float df = 0.001;
float max_f = 0.1;
void Ubpa::Simulate::Blow()
{
	int m, s;
	m = positions.size();
	s = edgelist.size() / 2;
	bool s_t = true;
	for (int i = 0; i < m; i++)
	{
		if (f_ext(i * 3 + 0) <= 0.0 && t == Tt) {
			f_ext(i * 3 + 0) = max_f;
			f_ext(i * 3 + 2) = max_f;
			s_t = true;
		}
		else if (f_ext(i * 3 + 0) <= 0.0) {
			s_t = false;
		}
		else {
			f_ext(i * 3 + 0) -= df;
			f_ext(i * 3 + 2) -= df;
			s_t = true;
		}
	}
	if (s_t) t = 0;
	else t++;
}

void Ubpa::Simulate::Bounce()
{
	int m, s;
	m = positions.size();
	s = edgelist.size() / 2;
	for (int i = 0; i < m; i++)
	{
		if (x(i * 3 + 1) <= 0.0) {
			if (v(i * 3 + 1) < 0.0) {
				x(i * 3 + 1) = 0.0;
				v(i * 3 + 1) = -v(i * 3 + 1);
				f_ext(i * 3 + 1) = -f_int(i * 3 + 1);
			}
		}
		else {
			f_ext(i * 3 + 1) = -My_g * this->mass[i];
		}
	}
}

bool Simulate::Init() {
	is_fix_ = true;
	is_blow_ = false;
	is_bounce_ = false;
	//Clear();
	int m, s;
	m = positions.size();
	s = edgelist.size() / 2;
	this->origin_length.resize(s);
	this->k.resize(s);
	this->velocity.resize(m);
	this->mass.resize(m);
	
	M.resize(3 * m, 3 * m);
	M_inv.resize(3 * m, 3 * m);
	L.resize(3 * m, 3 * m);
	J.resize(3 * m, 3 * s);
	vector<T> tripletList_M;
	vector<T> tripletList_M_INV;
	vector<T> tripletList_L;
	vector<T> tripletList_J;
	VectorXf diag_L(m);
	diag_L.setZero();
	d.resize(3 * s);
	x.resize(3 * m);
	y.resize(3 * m);
	v.resize(3 * m);
	f_int.resize(3 * m);
	f_ext.resize(3 * m);
	f_int.setZero();
	f_ext.setZero();
	//springs
	for (int i = 0; i < s; i++) 
	{
		auto dd = this->positions[this->edgelist[2 * i]] - this->positions[this->edgelist[2 * i + 1]];
		this->origin_length[i] = dd.norm();
		dd.normalize_self();
		dd *= this->origin_length[i];
		this->k[i] = 50.0;
		
		for (int j = 0; j < 3; j++) {
			d(i * 3 + j) = dd[j];
			tripletList_L.push_back(T(this->edgelist[2 * i] * 3 + j, this->edgelist[2 * i + 1] * 3 + j, -this->k[i]));
			tripletList_L.push_back(T(this->edgelist[2 * i + 1] * 3 + j, this->edgelist[2 * i] * 3 + j, -this->k[i]));
			tripletList_J.push_back(T(this->edgelist[2 * i] * 3 + j, i * 3 + j, this->k[i]));
			tripletList_J.push_back(T(this->edgelist[2 * i + 1] * 3 + j, i * 3 + j, -this->k[i]));
		}
		diag_L(this->edgelist[2 * i]) += this->k[i];
		diag_L(this->edgelist[2 * i + 1]) += this->k[i];
	}

	//points
	for (int i = 0; i < m; i++)
	{
		this->mass[i] = 0.01;
		for (int j = 0; j < 3; j++)
		{
			this->velocity[i][j] = 0.0;
			v(i * 3 + j) = this->velocity[i][j];
			x(i * 3 + j) = this->positions[i][j];
			tripletList_L.push_back(T(i * 3 + j, i * 3 + j, diag_L(i)));
			tripletList_M.push_back(T(i * 3 + j, i * 3 + j, this->mass[i]));
			tripletList_M_INV.push_back(T(i * 3 + j, i * 3 + j, 1 / this->mass[i]));
		}
		f_ext(i * 3 + 1) = -My_g * this->mass[i];
		f_ext(i * 3 + 0) = 0;
		f_ext(i * 3 + 2) = 0;
		//v(i * 3 + 1) = 7.0;
	}

	L.setFromTriplets(tripletList_L.begin(), tripletList_L.end());
	L.makeCompressed();
	M.setFromTriplets(tripletList_M.begin(), tripletList_M.end());
	M.makeCompressed();
	M_inv.setFromTriplets(tripletList_M_INV.begin(), tripletList_M_INV.end());
	M_inv.makeCompressed();
	J.setFromTriplets(tripletList_J.begin(), tripletList_J.end());
	J.makeCompressed();

	Fix(is_fix_);
	SparseMatrixType G;
	G = K * (M + h * h * L) * K.transpose();
	solver.compute(G);

	return true;
}

bool Simulate::Run() {
	SimulateOnce();
	// half-edge structure -> triangle mesh

	return true;
}

void Ubpa::Simulate::SetLeftFix()
{
	//固定网格x坐标最小点
	fixed_id.clear();
	double x = 100000;
	for (int i = 0; i < positions.size(); i++)
	{
		if (positions[i][0] < x)
		{
			x = positions[i][0];
		}
	}

	for (int i = 0; i < positions.size(); i++)
	{
		if (abs(positions[i][0] - x) < 1e-5)
		{
			fixed_id.push_back(i);
		}
	}

	Init();
}

void Ubpa::Simulate::CheckBlow()
{
	is_blow_ = !is_blow_;
	if (is_blow_) {
		int m;
		m = positions.size();
		for (int i = 0; i < m; i++)
		{
			f_ext(i * 3 + 0) = max_f;
			f_ext(i * 3 + 2) = max_f;
		}
	}
	else {
		int m;
		m = positions.size();
		for (int i = 0; i < m; i++)
		{
			f_ext(i * 3 + 0) = 0;
			f_ext(i * 3 + 2) = 0;
		}
	}
}

void Ubpa::Simulate::CheckBounce()
{
	is_bounce_ = !is_bounce_;
	if (is_bounce_) {
		int m;
		m = positions.size();
		for (int i = 0; i < m; i++)
		{
			v(i * 3 + 1) = 7.0;
		}
	}
}

void Ubpa::Simulate::CheckFix()
{
	is_fix_ = !is_fix_;
	Fix(is_fix_);
}

void Simulate::SimulateOnce() {
	// TODO
	//cout << "WARNING::Simulate::SimulateOnce:" << endl;
	//		<< "\t" << "not implemented" << endl;
	Optimize();
	UpdateParameter();
	if(is_bounce_) Bounce();
	for (int i = 0; i < positions.size(); i++)
	{
		for (int j = 0; j < 3; j++) {
			this->velocity[i][j] = v(i * 3 + j);
			this->positions[i][j] = x(i * 3 + j);
		}
	}
	if (is_blow_) Blow();
}

void Ubpa::Simulate::Fix(bool bv)
{
	int m, n, s;
	n = positions.size();
	s = edgelist.size() / 2;
	fix.clear();
	if (bv) {
		float max_y = INT_MIN;
		for (int i = 0; i < n; i++) {
			if (this->positions[i][1] > max_y) max_y = this->positions[i][1];
		}
		for (int i = 0; i < n; i++) {
			if (this->positions[i][1] == max_y) fix.push_back(i);
		}
		m = n - fix.size();
		vector<T> tripletList_K;
		K.resize(3 * m, 3 * n);
		int j = 0;
		for (int i = 0; i < n; i++) {
			if (this->positions[i][1] != max_y) {
				tripletList_K.push_back(T(j * 3 + 0, i * 3 + 0, 1));
				tripletList_K.push_back(T(j * 3 + 1, i * 3 + 1, 1));
				tripletList_K.push_back(T(j * 3 + 2, i * 3 + 2, 1));
				j++;
			}
		}
		K.setFromTriplets(tripletList_K.begin(), tripletList_K.end());
		K.makeCompressed();
	}
	else {
		K.resize(3 * n, 3 * n);
		K.setIdentity();
	}
	b = x - K.transpose() * K * x;
	xf = K * x;
}

void Ubpa::Simulate::Optimize()
{
	int N = 10;
	float E = 1E-6;
	float x_norm = 0;
	y = x + h * v + h * h * M_inv * f_ext;
	//cout << y.transpose() << endl;
	VectorXf a = M * y - (M + h * h * L) * b;
	for (int i = 0; i < N; i++) {
		//Optimize x
		VectorXf bb;
		bb = K * (h * h * J * d + a);
		//cout << x.transpose() * L * x - 2 * x.transpose() * J * d + d.transpose() * d << endl;
		xf = solver.solve(bb);
		x = K.transpose() * xf + b;
		
		//Optimize d
		for (int i = 0; i < this->edgelist.size() / 2; i++)
		{
			vecf3 dd;
			for (int j = 0; j < 3; j++) {
				dd[j] = x(this->edgelist[2 * i] * 3 + j) - x(this->edgelist[2 * i + 1] * 3 + j);
			}
			dd.normalize_self();
			dd *= this->origin_length[i];
			for (int j = 0; j < 3; j++) {
				d(i * 3 + j) = dd[j];
			}
		}

		//Loss
		float x_n = x.norm();
		float loss = x_n - x_norm;
		if (loss < E) break;
		x_norm = x_n;
	}
}

void Ubpa::Simulate::UpdateParameter()
{
	//x,d updated through Optimization
	f_int.setZero();
	for (int i = 0; i < this->edgelist.size()/2; i++)
	{
		float dx = (x(this->edgelist[2 * i] * 3 + 0) - x(this->edgelist[2 * i + 1] * 3 + 0) - d(i * 3 + 0)) * (x(this->edgelist[2 * i] * 3 + 0) - x(this->edgelist[2 * i + 1] * 3 + 0) - d(i * 3 + 0)) + (x(this->edgelist[2 * i] * 3 + 1) - x(this->edgelist[2 * i + 1] * 3 + 1) - d(i * 3 + 1)) * (x(this->edgelist[2 * i] * 3 + 1) - x(this->edgelist[2 * i + 1] * 3 + 1) - d(i * 3 + 1)) + (x(this->edgelist[2 * i] * 3 + 2) - x(this->edgelist[2 * i + 1] * 3 + 2) - d(i * 3 + 2)) * (x(this->edgelist[2 * i] * 3 + 2) - x(this->edgelist[2 * i + 1] * 3 + 2) - d(i * 3 + 2));
		dx = sqrt(dx);
		for (int j = 0; j < 3; j++) {
			f_int(this->edgelist[2 * i] * 3 + j) += k[i] * dx * d(i * 3 + j);
			f_int(this->edgelist[2 * i + 1] * 3 + j) -= k[i] * dx * d(i * 3 + j);
		}
	}

	for (int i = 0; i < this->positions.size(); i++)
	{
		for (int j = 0; j < 3; j++) {
			v(i * 3 + j) = (x(i * 3 + j) - this->positions[i][j]) / h;
		}
	}
	//cout << v.transpose() << endl;
}
