#pragma once

#include "Plane.h"
#include "Particle.h"
#include <glm/glm.hpp>
#include <random>
#include <vector>
#include <cmath>

#include "ctpl_stl.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace glm;

constexpr float GAS_CONST = 2000.0f;
constexpr float DT = 0.0008f;
constexpr float BOUND_DAMPING = -0.5f;
constexpr int seed = 12345678;

inline float randFloat()
{
	static std::mt19937 e2(seed);
	static std::uniform_real_distribution<float> dist(0.0f, 1.0f);
	return dist(e2);
}

struct SimulationParams
{
	float particleRadius;
	float smoothingRadius;
	float smoothingRadius2;
	float restDensity;
	float gravityMult;

	float particleMass;
	float particleViscosity;
	float particleDrag;

	SimulationParams() = default;

	inline SimulationParams(float pRadius, float sRadius, float restDens, float gravMult, float pMass, float pVisc,
							float pDrag)
	{
		particleRadius = pRadius;
		smoothingRadius = sRadius;
		smoothingRadius2 = sRadius * sRadius;
		restDensity = restDens;
		gravityMult = gravMult;
		particleMass = pMass;
		particleViscosity = pVisc;
		particleDrag = pDrag;
	}
};

class Simulator
{
  public:
	explicit Simulator(int rowSize, vec3 delta = vec3(0.0f));

	~Simulator() = default;

	size_t addParams(const SimulationParams &params);

	void addParticles(size_t N, size_t parameterID);

	void addPlane(Plane plane);

	inline const std::vector<Particle> &getParticles() const { return m_Particles; }

	inline size_t getParticleCount() const { return m_Particles.size(); }

	inline const std::vector<SimulationParams> &getSimParams() const { return m_Params; }

	inline const std::vector<Plane> &getPlanes() const { return m_Collider; }

	inline void update(float timestep)
	{
		computeDensityPressure();
		computeForces();
		integrate(timestep);
		computeCollisions();
	}

	void reset();

  private:
	static bool intersect(const Plane &collider, const vec3 &position, float radius, vec3 &penetrationNormal,
						  vec3 &penetrationPos, float &penetrationLength);
	static vec3 dampVelocity(const Plane &collider, const vec3 &velocity, const vec3 &penetrationNormal, float drag);

	void computeCollisions();
	void integrate(float timestep);

	void computeDensityPressure();
	void computeForces();

  private:
	std::vector<Plane> m_Collider = {};
	std::vector<Particle> m_Particles = {};
	std::vector<SimulationParams> m_Params;
	vec3 m_Delta;
	int m_RowSize;
	ctpl::thread_pool* m_Pool;
	int m_ThreadCount = std::thread::hardware_concurrency();
	std::vector<std::future<void>> m_Jobs;
};
