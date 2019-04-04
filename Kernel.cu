#include "Kernel.cuh"

#define GLM_FORCE_PURE
#include <device_launch_parameters.h>
#include <glm/glm.hpp>

using namespace glm;

#ifndef __CUDACC__
int atomicAdd(int *i, int c){};
#endif

__device__ bool plane_intersect(const Plane &collider, const vec3 &position, float radius, vec3 &penetrationNormal,
								vec3 &penetrationPos, float &penetrationLength)
{
	vec3 projection = collider.position - position;

	penetrationNormal = cross(collider.right, collider.forward);
	penetrationLength = abs(dot(projection, penetrationNormal)) - radius / 2.0f;
	penetrationPos = collider.position - projection;

	return penetrationLength < 0.0f && abs(dot(projection, collider.right)) < collider.size.x &&
		   abs(dot(projection, collider.forward)) < collider.size.y;
}

__device__ vec3 dampVelocity(const Plane &collider, const vec3 &velocity, const vec3 &penetrationNormal, float drag)
{
	const vec3 newV = dot(velocity, penetrationNormal) * penetrationNormal * BOUND_DAMPING +
					  dot(velocity, collider.right) * collider.right * drag +
					  dot(velocity, collider.forward) * collider.forward * drag;
	const vec3 forward = vec3(0.0f, 0.0f, 1.0f);
	const vec3 right = vec3(1.0f, 0.0f, 0.0f);
	const vec3 up = vec3(0.0f, 1.0f, 0.0);

	return dot(newV, forward) * forward + dot(newV, right) * right + dot(newV, up) * up;
}

__device__ inline int getGridPos(int x, int y, int z) { return x + gridDimX * (y + gridDimZ * z); }
__device__ inline int getGridIdx(int gridPos, int i) { return gridPos + i * bucketCapacity; }

__device__ i32vec3 getParticleGridPosition(const vec3 &position, const vec3 &worldMin, const vec3 &worldMax)
{
	const vec3 normalized_position = (position - worldMin) / (worldMax - worldMin);

	// Find buckets to place particles in
	// Could be written vectorized but explicit is easier to understand
	int bucketX = min(max(0, int(normalized_position.x * gridDimX)), gridDimX - 1);
	int bucketY = min(max(0, int(normalized_position.y * gridDimY)), gridDimY - 1);
	int bucketZ = min(max(0, int(normalized_position.z * gridDimZ)), gridDimZ - 1);
	return {bucketX, bucketY, bucketZ};
}

__global__ void clearGrid(int *gridCounter, int *gridIndices)
{
	unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned idy = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned idz = blockIdx.z * blockDim.z + threadIdx.z;
	if (idx >= gridDimX || idy >= gridDimY || idz >= gridDimZ)
		return;

	const int gridPos = getGridPos(idx, idy, idz);
	gridCounter[gridPos] = 0;
}

__global__ void buildGrid(Particle *particles, int particleCount, int *gridCounter, int *gridIndices, vec3 worldMin,
						  vec3 worldMax)
{
	const int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx >= particleCount)
		return;

	const auto &particle = particles[idx];
	const auto pGridPos = getParticleGridPosition(particle.position, worldMin, worldMax);

	const int gridPos = getGridPos(pGridPos.x, pGridPos.y, pGridPos.z);
	const int i = atomicAdd(&gridCounter[gridPos], 1);

	if (i < bucketCapacity)
		gridIndices[getGridIdx(gridPos, i)] = idx;
}

__global__ void computeDensityPressure(Particle *particles, SimulationParams params, Plane *planes, int particleCount,
									   int planeCount, float timestep, int *gridCounter, int *gridIndices,
									   vec3 worldMin, vec3 worldMax)
{
	const int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx >= particleCount)
		return;

	auto &pi = particles[idx];
	const i32vec3 piGridPos = getParticleGridPosition(pi.position, worldMin, worldMax);
	const int gridIdx = getGridPos(piGridPos.x, piGridPos.y, piGridPos.z);
	const int counter = gridCounter[gridIdx];

	for (int i = 0; i < counter && i < bucketCapacity; i++)
	{
		const auto &pj = particles[gridIndices[getGridIdx(gridIdx, i)]];
		const vec3 rij = pj.position - pi.position;
		const float r2 = dot(rij, rij);
		if (r2 < params.smoothingRadiusPow2)
		{
			const float temp = params.smoothingRadiusPow2 - r2;
			const float smoothingRadiusMinR2pow3 = temp * temp * temp;
			pi.density += params.particleMass * (315.0f / (64.0f * glm::pi<float>() * params.smoothingRadiusPow9)) *
						  smoothingRadiusMinR2pow3;
		}
	}

	pi.pressure = GAS_CONST * (pi.density * params.restDensity);
}
__global__ void computeForces(Particle *particles, SimulationParams params, Plane *planes, int particleCount,
							  int planeCount, float timestep, int *gridCounter, int *gridIndices, vec3 worldMin,
							  vec3 worldMax)
{
	const int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx >= particleCount)
		return;

	auto &pi = particles[idx];
	const i32vec3 piGridPos = getParticleGridPosition(pi.position, worldMin, worldMax);

	vec3 forcePressure = vec3(0.0f);
	vec3 forceViscosity = vec3(0.0f);

	const i32vec3 begin = max(i32vec3{0, 0, 0}, piGridPos - 1);
	const i32vec3 end = min(i32vec3{gridDimX - 1, gridDimY - 1, gridDimZ - 1}, piGridPos + 1);

	for (int i = begin.x; i <= end.x; i++)
	{
		for (int j = begin.y; j <= end.y; j++)
		{
			for (int k = begin.z; k <= end.z; k++)
			{
				const int gridj = getGridPos(i, j, k);
				for (int pjIdx = 0; pjIdx < gridCounter[gridj] && pjIdx < bucketCapacity; pjIdx++)
				{
					const auto &pj = particles[gridIndices[getGridIdx(gridj, pjIdx)]];
					if (&pi == &pj)
						continue;

					const vec3 rij = pj.position - pi.position;
					const float r2 = dot(rij, rij);
					if (r2 < params.smoothingRadiusPow2)
					{
						const float r = sqrtf(r2);
						const float sRmin_r = params.smoothingRadius - r;
						const vec3 rijNorm = rij / r;
						const float pressureSum = pi.pressure + pj.pressure;
						const vec3 delta_v = pj.velocity - pi.velocity;
						const float fourtyFiveOverPI_SR6 = 45.0f / (glm::pi<float>() * params.smoothingRadiusPow6);

						forcePressure += -rijNorm * params.particleMass * pressureSum / (2.0f * pj.density) *
										 fourtyFiveOverPI_SR6 * sRmin_r * sRmin_r;
						forceViscosity += params.particleViscosity * params.particleMass * delta_v / pj.density *
										  fourtyFiveOverPI_SR6 * sRmin_r;
					}
				}
			}
		}
	}

	const vec3 forceGravity = -params.gravity * pi.density * params.gravityMult;
	pi.forcePhysic = forcePressure + forceViscosity + forceGravity;
}

__global__ void integrateAndCollisions(Particle *particles, SimulationParams params, Plane *planes, int particleCount,
									   int planeCount, float timestep, int *gridCounter, int *gridIndices,
									   vec3 worldMin, vec3 worldMax)
{
	const int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx >= particleCount)
		return;
	auto &pi = particles[idx];

	pi.velocity += timestep * pi.forcePhysic / pi.density;
	pi.position += timestep * pi.velocity;

	for (int i = 0; i < planeCount; i++)
	{
		vec3 penetrationNormal, penetrationPosition;
		float penetrationLength;
		if (plane_intersect(planes[i], pi.position, params.particleRadius, penetrationNormal, penetrationPosition,
							penetrationLength))
		{
			pi.velocity = dampVelocity(planes[i], pi.velocity, penetrationNormal, 1.0f - params.particleDrag);
			pi.position = penetrationPosition - penetrationNormal * fabs(penetrationLength);
		}
	}
}

// buildGrid();

void LaunchKernels(Particle *particles, int particleCount, SimulationParams params, Plane *planes, int planeCount,
				   float timestep, int *gridCounts, int *gridIndices, vec3 worldMin, vec3 worldMax)
{
	dim3 dimBlock = dim3(4, 4, 4);
	dim3 dimGrid = dim3(gridDimX + dimBlock.x - 1, gridDimY + dimBlock.y - 1, gridDimZ + dimBlock.z - 1);

	clearGrid<<<dimGrid, dimBlock>>>(gridCounts, gridIndices);

	dimBlock = dim3(16 * 16, 1, 1);
	dimGrid = dim3(particleCount + dimBlock.x - 1, 1, 1);
	buildGrid<<<dimGrid, dimBlock>>>(particles, particleCount, gridCounts, gridIndices, worldMin, worldMax);

	computeDensityPressure<<<dimGrid, dimBlock>>>(particles, params, planes, particleCount, planeCount, DT, gridCounts,
												  gridIndices, worldMin, worldMax);
	computeForces<<<dimGrid, dimBlock>>>(particles, params, planes, particleCount, planeCount, DT, gridCounts,
										 gridIndices, worldMin, worldMax);
	cudaDeviceSynchronize();

	integrateAndCollisions<<<dimGrid, dimBlock>>>(particles, params, planes, particleCount, planeCount, DT, gridCounts,
												  gridIndices, worldMin, worldMax);

	cudaDeviceSynchronize();
}