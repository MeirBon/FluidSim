#include "Simulator.h"

constexpr float PI = float(M_PI);

#define THREADED 1
#define USE_GRID 1

// Get the position in a grid of a particle. Returns {grid.x,grid.y,grid.z} vector of indices
// of the grid in which the particle resides.
i32vec3 Simulator::getParticleGridPosition(glm::vec3 position)
{
	glm::vec3 normalized_position = (position - worldMin) / (worldMax - worldMin);

	// Find buckets to place particles in
	// Could be written vectorized but explicit is easier to understand
	int bucketX = std::max(0, int(normalized_position.x * gridDimX));
	int bucketY = std::max(0, int(normalized_position.y * gridDimY));
	int bucketZ = std::max(0, int(normalized_position.z * gridDimZ));

	// Handle edgecase where normalized_position = 1.
	if (bucketX >= gridDimX)
		bucketX = gridDimX - 1;
	if (bucketY >= gridDimY)
		bucketY = gridDimY - 1;
	if (bucketZ >= gridDimZ)
		bucketZ = gridDimZ - 1;

	return {bucketX, bucketY, bucketZ};
}

Simulator::Simulator(int rowSize, vec3 delta) : m_RowSize(rowSize), m_Delta(delta)
{
	m_Pool = new ctpl::thread_pool(m_ThreadCount);
}

size_t Simulator::addParams(const SimulationParams &params)
{
	const auto pid = m_Params.size();
	m_Params.push_back(params);
	return pid;
}

void Simulator::buildGrid()
{
	for (int i = 0; i < gridDimX; ++i)
	{
		for (int j = 0; j < gridDimY; ++j)
		{
			for (int k = 0; k < gridDimZ; ++k)
			{
				particleGrid[i][j][k].clear();
			}
		}
	}

	for (int i = 0; i < m_Particles.size(); ++i)
	{
		// Get data between 0 and 1
		const auto &particle = m_Particles[i];
		i32vec3 buckets = getParticleGridPosition(m_Particles[i].position);

#if 0
		if (buckets.x < 0 || buckets.x >= gridDimX || buckets.y < 0 || buckets.y >= gridDimY || buckets.z < 0 ||
			buckets.z >= gridDimZ)
		{
			perror("particle grid position out of bounds");
		}
#endif
		particleGrid[buckets.x][buckets.y][buckets.z].push_back(i);
	}
}

void Simulator::setParticleGridBounds(glm::vec3 minPoint, glm::vec3 maxPoint)
{
	worldMin = minPoint;
	worldMax = maxPoint;
}

void Simulator::addParticles(size_t N, size_t parameterID)
{
	const auto &params = m_Params[parameterID];
	const float oneOverRowSize = 1.0f / float(m_RowSize);

	for (int i = 0; i < N; i++)
	{
		// TODO(Dan): Jitter not used?
		const float jitter = randFloat() * params.particleRadius * 0.1f;
		const float x = (i % m_RowSize) + (randFloat() * 0.2f - 0.1f);

#if 0 // Limit y so that grid can have well specified boundaries
	  // const float y = 2.0f + (float(i) * oneOverRowSize * oneOverRowSize) * 1.1f;
#else
		const float y = 10.0f;
#endif
		const auto temp = float(mod(double(i) * (double)oneOverRowSize, double(m_RowSize)));
		const float z = temp * (randFloat() * 0.2f - 0.1f);
		m_Particles.emplace_back(vec3(x, y, z) + m_Delta, parameterID);
	}
}

void Simulator::addPlane(Plane plane) { m_Collider.push_back(plane); }

void Simulator::reset()
{
	if (m_Particles.empty())
		return;

	auto oldParticles = m_Particles;
	m_Particles.clear();

	int paramsId = m_Particles[0].parameterID;
	int i = 0;

	while (i < oldParticles.size())
	{
		size_t count = 0;

		while (i < oldParticles.size() && oldParticles[i].parameterID == paramsId)
		{
			i++;
			count++;
		}

		addParticles(count, paramsId);
		paramsId = oldParticles[i].parameterID;
	}
}

bool Simulator::intersect(const Plane &collider, const vec3 &position, float radius, vec3 &penetrationNormal,
						  vec3 &penetrationPos, float &penetrationLength)
{
	vec3 projection = collider.position - position;

	penetrationNormal = cross(collider.right, collider.forward);
	penetrationLength = abs(dot(projection, penetrationNormal)) - radius / 2.0f;
	penetrationPos = collider.position - projection;

	return penetrationLength < 0.0f && abs(dot(projection, collider.right)) < collider.size.x &&
		   abs(dot(projection, collider.forward)) < collider.size.y;
}

vec3 Simulator::dampVelocity(const Plane &collider, const vec3 &velocity, const vec3 &penetrationNormal, float drag)
{
	const vec3 newV = dot(velocity, penetrationNormal) * penetrationNormal * BOUND_DAMPING +
					  dot(velocity, collider.right) * collider.right * drag +
					  dot(velocity, collider.forward) * collider.forward * drag;
	const vec3 forward = vec3(0.0f, 0.0f, 1.0f);
	const vec3 right = vec3(1.0f, 0.0f, 0.0f);
	const vec3 up = vec3(0.0f, 1.0f, 0.0);

	return dot(newV, forward) * forward + dot(newV, right) * right + dot(newV, up) * up;
}

void Simulator::computeCollisions()
{
#if THREADED
	for (int tId = 0; tId < m_ThreadCount; tId++)
	{
		m_Jobs.push_back(m_Pool->push([this, tId](int) {
			for (int i = tId; i < m_Particles.size(); i += m_ThreadCount)
			{
				auto &p = m_Particles[i];
				for (const auto &plane : m_Collider)
				{
					vec3 penetrationNormal, penetrationPosition;
					float penetrationLength;
					const auto &params = m_Params[p.parameterID];

					if (intersect(plane, p.position, params.particleRadius, penetrationNormal, penetrationPosition,
								  penetrationLength))
					{
						p.velocity = dampVelocity(plane, p.velocity, penetrationNormal, 1.0f - params.particleDrag);
						p.position = penetrationPosition - penetrationNormal * abs(penetrationLength);
					}
				}
			}
		}));
	}

	for (auto &job : m_Jobs)
		job.get();
	m_Jobs.clear();
#else
	for (auto &p : m_Particles)
	{
		for (const auto &plane : m_Collider)
		{
			vec3 penetrationNormal, penetrationPosition;
			float penetrationLength;
			const auto &params = m_Params[p.parameterID];

			if (intersect(plane, p.position, params.particleRadius, penetrationNormal, penetrationPosition,
						  penetrationLength))
			{
				p.velocity = dampVelocity(plane, p.velocity, penetrationNormal, 1.0f - params.particleDrag);
				p.position = penetrationPosition - penetrationNormal * abs(penetrationLength);
			}
		}
	}
#endif
}

void Simulator::integrate(float timestep)
{
#if THREADED
	for (int tId = 0; tId < m_ThreadCount; tId++)
	{
		m_Jobs.push_back(m_Pool->push([this, tId](int) {
			for (int i = tId; i < m_Particles.size(); i += m_ThreadCount)
			{
				auto &p = m_Particles[i];
				p.velocity += DT * p.forcePhysic / p.density;
				p.position += DT * p.velocity;
			}
		}));
	}

	for (auto &job : m_Jobs)
		job.get();
	m_Jobs.clear();
#else
	for (int i = 0; i < m_Particles.size(); i++)
	{
		auto &p = m_Particles[i];
		p.velocity += DT * p.forcePhysic / p.density;
		p.position += DT * p.velocity;
	}
#endif
}

void Simulator::computeDensityPressure()
{
#if THREADED
	for (int tId = 0; tId < m_ThreadCount; tId++)
	{
		m_Jobs.push_back(m_Pool->push([this, tId](int) {
			for (int i = tId; i < m_Particles.size(); i += m_ThreadCount)
			{
				auto &pi = m_Particles[i];
				pi.density = 0.0f;
				const auto &paramsi = m_Params[pi.parameterID];
				const float smoothingRadiusPow9 = paramsi.smoothingRadiusPow2 * paramsi.smoothingRadiusPow2 *
												  paramsi.smoothingRadiusPow2 * paramsi.smoothingRadiusPow2 *
												  paramsi.smoothingRadius;
				i32vec3 particleGridSlot = getParticleGridPosition(pi.position);

				i32vec3 begin = max(i32vec3{0, 0, 0}, particleGridSlot - 1);
				i32vec3 end = min(i32vec3{gridDimX - 1, gridDimY - 1, gridDimZ - 1}, particleGridSlot + 1);

				for (int i = begin.x; i <= end.x; ++i)
				{
					for (int j = begin.y; j <= end.y; ++j)
					{
						for (int k = begin.z; k <= end.z; ++k)
						{
							for (const auto &pjIndex : particleGrid[i][j][k])
							{
								// TODO(Dan): Why don't we ignore  pj == pi here?
								const auto &pj = m_Particles[pjIndex];
								const vec3 rij = pj.position - pi.position;
								const float r2 = glm::dot(rij, rij);
								if (r2 < paramsi.smoothingRadiusPow2)
								{
									const float temp = paramsi.smoothingRadiusPow2 - r2;
									const float smoothingRadiusMinR2pow3 = temp * temp * temp;
									pi.density += paramsi.particleMass * (315.0f / (64.0f * PI * smoothingRadiusPow9)) *
												  smoothingRadiusMinR2pow3;
								}
							}
						}
					}
				}
				pi.pressure = GAS_CONST * (pi.density * paramsi.restDensity);
			}
		}));
	}

	for (auto &job : m_Jobs)
		job.get();
	m_Jobs.clear();
#else
	for (auto &pi : m_Particles)
	{
		pi.density = 0.0f;
		const auto &paramsi = m_Params[pi.parameterID];
		const float smoothingRadius9 = paramsi.smoothingRadiusPow2 * paramsi.smoothingRadiusPow2 *
									   paramsi.smoothingRadiusPow2 * paramsi.smoothingRadiusPow2 *
									   paramsi.smoothingRadius;

		i32vec3 particleGridSlot = getParticleGridPosition(pi.position);

		i32vec3 begin = max(i32vec3{0, 0, 0}, particleGridSlot - 1);
		i32vec3 end = min(i32vec3{gridDimX - 1, gridDimY - 1, gridDimZ - 1}, particleGridSlot + 1);

#if !USE_GRID
		for (int i = 0; i < m_Particles.size(); ++i)
		{
			const auto &pj = m_Particles[i];
			const vec3 rij = pj.position - pi.position;
			const float r2 = glm::dot(rij, rij);

			if (r2 < paramsi.smoothingRadiusPow2)
			{
				const float temp = paramsi.smoothingRadiusPow2 - r2;
				const float smoothingRadiusMinR2pow3 = temp * temp * temp;
				pi.density +=
					paramsi.particleMass * (315.0f / (64.0f * PI * smoothingRadius9)) * smoothingRadiusMinR2pow3;
			}
		}
#else
		for (int i = begin.x; i <= end.x; ++i)
		{
			for (int j = begin.y; j <= end.y; ++j)
			{
				for (int k = begin.z; k <= end.z; ++k)
				{
					for (const auto &pjIndex : particleGrid[i][j][k])
					{
						const auto &pj = m_Particles[pjIndex];
						const vec3 rij = pj.position - pi.position;
						const float r2 = glm::dot(rij, rij);
						if (r2 < paramsi.smoothingRadiusPow2)
						{
							const float temp = paramsi.smoothingRadiusPow2 - r2;
							const float smoothingRadiusMinR2pow3 = temp * temp * temp;
							pi.density += paramsi.particleMass * (315.0f / (64.0f * PI * smoothingRadius9)) *
										  smoothingRadiusMinR2pow3;
						}
					}
				}
			}
		}
#endif
		pi.pressure = GAS_CONST * (pi.density * paramsi.restDensity);
	}
#endif
}

void Simulator::computeForces()
{
#if THREADED
	// TODO(Dan): Recheck grid forr this & threading
	for (int tId = 0; tId < m_ThreadCount; tId++)
	{
		m_Jobs.push_back(m_Pool->push([this, tId](int) {
			for (int i = tId; i < m_Particles.size(); i += m_ThreadCount)
			{
				auto &pi = m_Particles[i];
				vec3 forcePressure = vec3(0.0f);
				vec3 forceViscosity = vec3(0.0f);
				const auto &paramsi = m_Params[pi.parameterID];

				const float smoothingRadius6 =
					paramsi.smoothingRadiusPow2 * paramsi.smoothingRadiusPow2 * paramsi.smoothingRadiusPow2;
				i32vec3 particleGridSlot = getParticleGridPosition(pi.position);

				i32vec3 begin = max(i32vec3{0, 0, 0}, particleGridSlot - 1);
				i32vec3 end = min(i32vec3{gridDimX - 1, gridDimY - 1, gridDimZ - 1}, particleGridSlot + 1);

				for (int i = begin.x; i <= end.x; ++i)
				{
					for (int j = begin.y; j <= end.y; ++j)
					{
						for (int k = begin.z; k <= end.z; ++k)
						{
							for (const auto &pjIndex : particleGrid[i][j][k])
							{
								// TODO(Dan): Why don't we ignore  pj == pi here?
								const auto &pj = m_Particles[pjIndex];
								if (&pi == &pj)
									continue;

								const vec3 rij = pj.position - pi.position;
								const float r2 = glm::dot(rij, rij);

								if (r2 < paramsi.smoothingRadiusPow2)
								{
									const float r = sqrtf(r2);
									const float sRmin_r = paramsi.smoothingRadius - r;
									const vec3 rijNorm = rij / r;
									const float pressureSum = pi.pressure + pj.pressure;
									const vec3 delta_v = pj.velocity - pi.velocity;
									const float fourtyFiveOverPI_SR6 = 45.0f / (PI * smoothingRadius6);

									forcePressure += -rijNorm * paramsi.particleMass * pressureSum /
													 (2.0f * pj.density) * fourtyFiveOverPI_SR6 * sRmin_r * sRmin_r;
									forceViscosity += paramsi.particleViscosity * paramsi.particleMass * delta_v /
													  pj.density * fourtyFiveOverPI_SR6 * sRmin_r;
								}
							}
						}
					}
				}

				const vec3 forceGravity = vec3(0.0f, -9.81f, 0.0f) * pi.density * paramsi.gravityMult;
				pi.forcePhysic = forcePressure + forceViscosity + forceGravity;
			}
		}));
	}

	for (auto &job : m_Jobs)
		job.get();

	m_Jobs.clear();
#else
	for (int i = 0; i < m_Particles.size(); i++)
	{
		auto &pi = m_Particles[i];
		vec3 forcePressure = vec3(0.0f);
		vec3 forceViscosity = vec3(0.0f);
		const auto &paramsi = m_Params[pi.parameterID];

		const float smoothingRadius6 =
			paramsi.smoothingRadiusPow2 * paramsi.smoothingRadiusPow2 * paramsi.smoothingRadiusPow2;

		for (const auto &pj : m_Particles)
		{
			if (&pi == &pj)
				continue;

			const vec3 rij = pj.position - pi.position;
			const float r2 = glm::dot(rij, rij);

			if (r2 < paramsi.smoothingRadiusPow2)
			{
				const float r = sqrtf(r2);
				const float sRmin_r = paramsi.smoothingRadius - r;
				const vec3 rijNorm = rij / r;
				const float pressureSum = pi.pressure + pj.pressure;
				const vec3 delta_v = pj.velocity - pi.velocity;
				const float fourtyFiveOverPI_SR6 = 45.0f / (PI * smoothingRadius6);

				forcePressure += -rijNorm * paramsi.particleMass * pressureSum / (2.0f * pj.density) *
								 fourtyFiveOverPI_SR6 * sRmin_r * sRmin_r;
				forceViscosity += paramsi.particleViscosity * paramsi.particleMass * delta_v / pj.density *
								  fourtyFiveOverPI_SR6 * sRmin_r;
			}
		}

		const vec3 forceGravity = vec3(0.0f, -9.81f, 0.0f) * pi.density * paramsi.gravityMult;
		pi.forcePhysic = forcePressure + forceViscosity + forceGravity;
	}
#endif
}
