#include "PathTracer.h"

#include <UBL/Image.h>
#include <omp.h>
#include <iostream>

#define USE_OPENMP_
#include <thread>


using namespace Ubpa;
using namespace std;

PathTracer::PathTracer(const Scene* scene, const SObj* cam_obj, Image* img)
	: scene{ scene },
	bvh{ const_cast<Scene*>(scene) },
	img{ img },
	cam{ cam_obj->Get<Cmpt::Camera>() },
	ccs{ cam->GenCoordinateSystem(cam_obj->Get<Cmpt::L2W>()->value) }
{
	scene->Each([this](const Cmpt::Light* light) ->bool {
		if (!vtable_is<EnvLight>(light->light.get()))
			return true; // continue

		env_light = static_cast<const EnvLight*>(light->light.get());
		return false; // stop
	});

	// TODO: preprocess env_light here
}

inline void UpdateProgress(string str, float progress)
{
	int barWidth = 70;

	std::cout << str << ": [";
	int pos = barWidth * progress;
	for (int i = 0; i < barWidth; ++i) {
		if (i < pos) std::cout << "=";
		else if (i == pos) std::cout << ">";
		else std::cout << " ";
	}
	std::cout << "] " << int(progress * 100.0) << " %\r";
	std::cout.flush();
};


void PathTracer::Run() {
	img->SetAll(0.f);

	const int spp = 2; // samples per pixel
	auto start = std::chrono::system_clock::now();
#ifdef USE_OPENMP_
	omp_set_num_threads(std::min(spp,32));
	for (size_t j = 0; j < img->height; j++) {
		for (size_t i = 0; i < img->width; i++) {
			float sumr, sumg, sumb;
			sumr = sumg = sumb = 0.0;
#pragma omp parallel for reduction(+: sumr, sumg, sumb)
			for (int k = 0; k < spp; k++) {
				Intersectors intersectors;
				float u = (i + rand01<float>() - 0.5f) / img->width;
				float v = (j + rand01<float>() - 0.5f) / img->height;
				rayf3 r = cam->GenRay(u, v, ccs);
				rgbf Lo = Shade(intersectors, intersectors.clostest.Visit(&bvh, r), -r.dir, true);
				sumr += Lo.at(0) / spp;
				sumg += Lo.at(1) / spp;
				sumb += Lo.at(2) / spp;
			}
			img->At<rgbf>(i, j) += rgbf(sumr, sumg, sumb);
		}
		float progress = (j + 1) / float(img->height);
		UpdateProgress("Rendering", progress);
	}
	UpdateProgress("Rendering", 1.f);
#else
#ifdef NDEBUG
	const size_t core_num = std::thread::hardware_concurrency();
	auto work = [this, core_num, spp](size_t id) {
		Intersectors intersectors;
		for (size_t j = id; j < img->height; j += core_num) {
			for (size_t i = 0; i < img->width; i++) {
				for (size_t k = 0; k < spp; k++) {
					float u = (i + rand01<float>() - 0.5f) / img->width;
					float v = (j + rand01<float>() - 0.5f) / img->height;
					rayf3 r = cam->GenRay(u, v, ccs);
					rgbf Lo;
					do { Lo = Shade(intersectors, intersectors.clostest.Visit(&bvh, r), -r.dir, true); }
					while (Lo.has_nan());
					img->At<rgbf>(i, j) += Lo / float(spp);
				}
			}
			float progress = (j + 1) / float(img->height);
			cout << progress << endl;
		}
	};
	vector<thread> workers;
	for (size_t i = 0; i < core_num; i++)
		workers.emplace_back(work, i);
	for (auto& worker : workers)
		worker.join();
#else
	Intersectors intersectors;
	for (size_t j = 0; j < img->height; j++) {
		for (size_t i = 0; i < img->width; i++) {
			for (size_t k = 0; k < spp; k++) {
				float u = (i + rand01<float>() - 0.5f) / img->width;
				float v = (j + rand01<float>() - 0.5f) / img->height;
				rayf3 r = cam->GenRay(u, v, ccs);
				rgbf Lo;
				do { Lo = Shade(intersectors, intersectors.clostest.Visit(&bvh, r), -r.dir, true); }
				while (Lo.has_nan());
				img->At<rgbf>(i, j) += Lo / spp;
			}
		}
		float progress = (j + 1) / float(img->height);
		UpdateProgress("Rendering", progress);
	}
	UpdateProgress("Rendering", 1.f);
#endif //NDEBUG
#endif //USE_OPENMP_
	auto stop = std::chrono::system_clock::now();
	std::cout << "Render complete: \n";
	std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::hours>(stop - start).count() << " hours\n";
	std::cout << "          : " << std::chrono::duration_cast<std::chrono::minutes>(stop - start).count() << " minutes\n";
	std::cout << "          : " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n";
}

rgbf PathTracer::Shade(const Intersectors& intersectors, const IntersectorClosest::Rst& intersection, const vecf3& wo, bool last_bounce_specular) {
	// TODO: HW9 - Trace
	// [ Tips ]
	// - EnvLight::Radiance(<direction>), <direction> is pointing to environment light
	// - AreaLight::Radiance(<uv>)
	// - rayf3: point, dir, tmin, **tmax**
	// - Intersectors::visibility.Visit(&bvh, <rayf3>)
	//   - tmin = EPSILON<float>
	//   - tmax = distance to light - EPSILON<float>
	// - Intersectors::cloest.Visit(&bvh, <rayf3>)
	//   - tmin as default (EPSILON<float>)
	//   - tmax as default (FLT_MAX)
	//
	// struct IntersectorClosest::Rst {
	//	 bool IsIntersected() const noexcept { return sobj != nullptr; }
	//	 const SObj* sobj{ nullptr }; // intersection sobj
	//	 pointf3 pos; // intersection point's position
	//	 pointf2 uv; // texcoord
	//   normalf n; // normal, normalized
	//	 vecf3 tangent; // perpendicular to normal, normalized
	// };

	constexpr rgbf error_color = rgbf{ 1.f,0.f,1.f };
	rgbf todo_color = rgbf{ 0.f,1.f,0.f };
	constexpr rgbf zero_color = rgbf{ 0.f,0.f,0.f };

	if (!intersection.IsIntersected()) {
		if (last_bounce_specular && env_light != nullptr) {
			// TODO: environment light
			todo_color = env_light->Radiance(-wo);
			return todo_color;
		}
		else
			return zero_color;
	}

	if (!intersection.sobj->Get<Cmpt::Material>()) {
		auto light = intersection.sobj->Get<Cmpt::Light>();
		if(!light) return error_color;

		if (last_bounce_specular) { // avoid double-count
			auto area_light = dynamic_cast<const AreaLight*>(light->light.get());
			if (!area_light) return error_color;
			// TODO: area light
			todo_color = area_light->Radiance(intersection.uv);
			return todo_color;
		}else
			return zero_color;
	}
	
	rgbf L_dir{ 0.f };
	rgbf L_indir{ 0.f };

	scene->Each([=, &intersectors, &L_dir](const Cmpt::Light* light, const Cmpt::L2W* l2w, const Cmpt::SObjPtr* ptr) {
		// TODO: L_dir += ...
		// - use PathTracer::BRDF to get BRDF value
		SampleLightResult sample_light_rst = SampleLight(intersection, wo, light, l2w, ptr);
		if (sample_light_rst.pd <= 0)
			return;
		if (sample_light_rst.is_infinity) {
			// TODO: L_dir of environment light
			// - only use SampleLightResult::L, n, pd
			// - SampleLightResult::x is useless		
			rayf3 r = rayf3(intersection.pos, sample_light_rst.n.cast_to<vecf3>());
			Intersectors inter_0;
			IntersectorClosest::Rst inter_0_rst = inter_0.clostest.Visit(&bvh, r);
			if (inter_0.visibility.Visit(&bvh, r))
				L_dir += env_light->Radiance(inter_0_rst.uv) * BRDF(intersection, wo, r.dir) * r.dir.dot(intersection.n.cast_to<vecf3>()) / sample_light_rst.pd;
		}
		else {
			// TODO: L_dir of area light
			rayf3 r = rayf3(intersection.pos, -(intersection.pos - sample_light_rst.x).normalize());
			Intersectors inter_0;
			IntersectorClosest::Rst inter_0_rst = inter_0.clostest.Visit(&bvh, r);
			if (inter_0.visibility.Visit(&bvh, rayf3(intersection.pos, -(intersection.pos - sample_light_rst.x).normalize(), EPSILON<float>, (intersection.pos - inter_0_rst.pos).norm() - EPSILON<float>))) {
				L_dir += sample_light_rst.L * BRDF(intersection, wo, r.dir) * r.dir.dot(intersection.n.cast_to<vecf3>()) * r.dir.dot(-sample_light_rst.n.cast_to<vecf3>()) / (intersection.pos - inter_0_rst.pos).norm2() / sample_light_rst.pd; 
			}
		}
	});

	// TODO: Russian Roulette
	// - rand01<float>() : random in [0, 1)
	if (rand01<float>() < Russian_Roulette_) {
		// TODO: recursion
		// - use PathTracer::SampleBRDF to get wi and pd (probability density)
		// - use PathTracer::BRDF to get BRDF value
		vecf3 wi;
		float pdf;
		tie(wi,pdf) = SampleBRDF(intersection, wo);

		rayf3 r = rayf3(intersection.pos, wi);
		Intersectors inter;
		IntersectorClosest::Rst inter_rst = inter.clostest.Visit(&bvh, r);

		if (inter_rst.IsIntersected() && intersection.sobj->Get<Cmpt::Material>()) {
			auto light = intersection.sobj->Get<Cmpt::Light>();
			if (!light) {
				L_indir += Shade(inter, inter_rst, -wi, true) * BRDF(intersection, wo, wi) * wi.dot(intersection.n.cast_to<vecf3>()) / pdf / Russian_Roulette_;
			}
			else {
				auto area_light = dynamic_cast<const AreaLight*>(light->light.get());
				if (!area_light) L_indir += Shade(inter, inter_rst, -wi, last_bounce_specular) * BRDF(intersection, wo, wi) * wi.dot(intersection.n.cast_to<vecf3>()) / pdf / Russian_Roulette_;
			}
		}
	}

	todo_color = L_dir;
	return todo_color; // TODO: combine L_dir and L_indir

}

PathTracer::SampleLightResult PathTracer::SampleLight(IntersectorClosest::Rst intersection, const vecf3& wo, const Cmpt::Light* light, const Cmpt::L2W* l2w, const Cmpt::SObjPtr* ptr) {
	PathTracer::SampleLightResult rst;
	if (vtable_is<AreaLight>(light->light.get())) {
		auto area_light = static_cast<const AreaLight*>(light->light.get());
		auto geo = ptr->value->Get<Cmpt::Geometry>();
		if (!geo) return rst; // invalid
		if (!vtable_is<Square>(geo->primitive.get())) return rst; // not support

		auto Xi = uniform_in_square<float>(); // [0, 1] x [0, 1]
		pointf3 pos_light_space{ 2 * Xi[0] - 1, 0, 2 * Xi[1] - 1 };
		scalef3 s = l2w->WorldScale();
		float area = s[0] * s[1] * Square::area;

		rst.L = area_light->Radiance(Xi.cast_to<pointf2>());
		rst.pd = 1.f / area;
		rst.x = l2w->value * pos_light_space;
		rst.is_infinity = false;
		rst.n = l2w->UpInWorld().cast_to<normalf>();
	}
	else if (vtable_is<EnvLight>(light->light.get())) {
		rst.is_infinity = true;

		auto mat = intersection.sobj->Get<Cmpt::Material>();
		if (!mat) return rst; // invalid
		auto brdf = dynamic_cast<const stdBRDF*>(mat->material.get());
		if (!brdf) return rst; // not support

		// multi-importance sampling, MIS

		auto env_light = static_cast<const EnvLight*>(light->light.get());

		float metalness = brdf->Metalness(intersection.uv);
		float roughness = brdf->Roughness(intersection.uv);
		float lambda = metalness * (1 - stdBRDF::Alpha(roughness)); // 0 - 1
		float p_mat = 1 / (2 - lambda); // 0.5 - 1
		float pd_mat, pd_env;
		vecf3 wi;
		rgbf Le;

		if (rand01<float>() < p_mat) {
			tie(wi, pd_mat) = SampleBRDF(intersection, wo);
			Le = env_light->Radiance(wi);
			pd_env = env_light->PDF(wi); // TODO: use your PDF
		}
		else {
			tie(Le, wi, pd_env) = env_light->Sample(intersection.n); // TODO: use your sampling method
			matf3 surface_to_world = svecf::TBN(intersection.n.cast_to<vecf3>(), intersection.tangent);
			matf3 world_to_surface = surface_to_world.inverse();
			svecf s_wo = (world_to_surface * wo).cast_to<svecf>();
			svecf s_wi = (world_to_surface * wi).cast_to<svecf>();
			rgbf albedo = brdf->Albedo(intersection.uv);
			pd_mat = brdf->PDF(albedo, metalness, roughness, s_wi, s_wo);
		}

		rst.L = Le;
		rst.n = -wi.cast_to<normalf>();
		rst.pd = p_mat * pd_mat + (1 - p_mat) * pd_env;
		rst.x = pointf3{ std::numeric_limits<float>::max() };
	}
	return rst;
}

std::tuple<vecf3, float> PathTracer::SampleBRDF(IntersectorClosest::Rst intersection, const vecf3& wo) {
	auto mat = intersection.sobj->Get<Cmpt::Material>();
	if (!mat) return { vecf3{0.f}, 0.f };
	auto brdf = dynamic_cast<const stdBRDF*>(mat->material.get());
	if (!brdf) return { vecf3{0.f}, 0.f };

	matf3 surface_to_world = svecf::TBN(intersection.n.cast_to<vecf3>(), intersection.tangent);
	matf3 world_to_surface = surface_to_world.inverse();
	svecf s_wo = (world_to_surface * wo).cast_to<svecf>();

	rgbf albedo = brdf->Albedo(intersection.uv);
	float metalness = brdf->Metalness(intersection.uv);
	float roughness = brdf->Roughness(intersection.uv);

	auto [s_wi, pdf] = brdf->Sample(albedo, metalness, roughness, s_wo);
	vecf3 wi = surface_to_world * s_wi;

	return { wi,pdf };
}

rgbf PathTracer::BRDF(IntersectorClosest::Rst intersection, const vecf3& wi, const vecf3& wo) {
	auto mat = intersection.sobj->Get<Cmpt::Material>();
	if (!mat) return rgbf{ 1.f,0.f,1.f };
	auto brdf = dynamic_cast<const stdBRDF*>(mat->material.get());
	if (!brdf) return rgbf{ 1.f,0.f,1.f };

	matf3 surface_to_world = svecf::TBN(intersection.n.cast_to<vecf3>(), intersection.tangent);
	matf3 world_to_surface = surface_to_world.inverse();
	svecf s_wi = (world_to_surface * wi).cast_to<svecf>();
	svecf s_wo = (world_to_surface * wo).cast_to<svecf>();

	rgbf albedo = brdf->Albedo(intersection.uv);
	float metalness = brdf->Metalness(intersection.uv);
	float roughness = brdf->Roughness(intersection.uv);

	return brdf->BRDF(albedo, metalness, roughness, s_wi, s_wo);
}