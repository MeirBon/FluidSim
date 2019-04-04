#pragma once

#include "Buffer.h"
#include "Particle.h"
#include "Plane.h"
#include "VertexArray.h"

#include <cmath>
#include <glm/glm.hpp>
#include <random>
#include <vector>

#include <PolyVoxCore/MarchingCubesSurfaceExtractor.h>
#include <PolyVoxCore/SimpleVolume.h>
#include <PolyVoxCore/SurfaceMesh.h>

#include "ctpl_stl.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

template <>
inline PolyVox::DefaultMarchingCubesController<float>::DensityType
PolyVox::DefaultMarchingCubesController<float>::getThreshold()
{
	return .7f;
}

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
	friend class Simulator;
	// Particle size
	float particleRadius;
	// Radius of influence - from SPH paper
	float smoothingRadius;
	float restDensity;
	float gravityMult;

	float particleMass;
	float particleViscosity;
	float particleDrag;

	vec3 gravity = vec3(0.0f, 9.81, 0.0f);

	SimulationParams() = default;

	inline SimulationParams(float pRadius, float sRadius, float restDens, float gravMult, float pMass, float pVisc,
							float pDrag)
	{
		particleRadius = pRadius;
		smoothingRadius = sRadius;

		smoothingRadiusPow2 = powf(sRadius, 2.0f);
		smoothingRadiusPow6 = powf(sRadius, 6.0f);
		smoothingRadiusPow9 = powf(sRadius, 9.0f);

		restDensity = restDens;
		gravityMult = gravMult;
		particleMass = pMass;
		particleViscosity = pVisc;
		particleDrag = pDrag;
	}

	// Optimization variables
	float particleRadiusPow2;
	float smoothingRadiusPow2;
	float smoothingRadiusPow6;
	float smoothingRadiusPow9;
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

	inline std::vector<SimulationParams> &getSimParams() { return m_Params; }

	inline const std::vector<Plane> &getPlanes() const { return m_Collider; }

	inline void update(float timestep)
	{
		// Update every frame to accomodate for possible changing values.
		for (auto &params : m_Params)
		{
			params.particleRadiusPow2 = params.particleRadius * params.particleRadius;
			params.smoothingRadiusPow2 = params.smoothingRadius * params.smoothingRadius;
			params.smoothingRadiusPow6 =
				params.smoothingRadiusPow2 * params.smoothingRadiusPow2 * params.smoothingRadiusPow2;
			params.smoothingRadiusPow9 = params.smoothingRadiusPow2 * params.smoothingRadiusPow2 *
										 params.smoothingRadiusPow2 * params.smoothingRadiusPow2 *
										 params.smoothingRadius;
		}

		buildGrid();
		computeDensityPressure();
		computeForces();
		integrate(timestep);
		computeCollisions();
	}

	void reset();

	void setParticleGridBounds(glm::vec3 minPoint, glm::vec3 maxPoint);

	float calculateDensity(const vec3 &pos);

	void extractSurface(Shader &shader);

	inline const vec3 &getWorldMin() const { return worldMin; }

	inline float getVoxelScale() { return voxelResScale; }

  private:
	static bool intersect(const Plane &collider, const vec3 &position, float radius, vec3 &penetrationNormal,
						  vec3 &penetrationPos, float &penetrationLength);
	static vec3 dampVelocity(const Plane &collider, const vec3 &velocity, const vec3 &penetrationNormal, float drag);

	void computeCollisions();
	void integrate(float timestep);

	void computeDensityPressure();
	void computeForces();
	void buildGrid();

	i32vec3 getParticleGridPosition(glm::vec3 position);

	vec3 voxelIndexToWorldPos(int voxelX, int voxelY, int voxelZ) const;

	void fillVoxelVolume();

  private:
	static constexpr int gridDimX = 25, gridDimY = 15, gridDimZ = 25;

	std::vector<Plane> m_Collider = {};
	std::vector<Particle> m_Particles = {};
	std::vector<SimulationParams> m_Params;
	std::vector<int> particleGrid[gridDimX][gridDimY][gridDimZ];

	// Effectively AABB of grid
	vec3 worldMin{INFINITY, INFINITY, INFINITY};
	vec3 worldMax{-INFINITY, -INFINITY, -INFINITY};
	vec3 worldLengths{};

	vec3 m_Delta;
	int m_RowSize;
	ctpl::thread_pool *m_Pool;
	int m_ThreadCount = std::thread::hardware_concurrency();
	std::vector<std::future<void>> m_Jobs;

	PolyVox::SimpleVolume<float> *voxelVolume;
	PolyVox::SurfaceMesh<PolyVox::PositionMaterialNormal> surfaceMesh;
	PolyVox::MarchingCubesSurfaceExtractor<PolyVox::SimpleVolume<float>> *surfaceExtractor;

	GLuint fluidVBO, fluidEBO, VAO;
	bool firstLaunch = true;
	const float voxelResScale = 2.0f;

	std::vector<float> poly6LookupTable{};
};
