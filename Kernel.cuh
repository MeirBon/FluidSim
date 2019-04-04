#pragma once

#include "Particle.h"
#include "Plane.h"
#include "Simulator.h"

#include <cuda_runtime.h>

void LaunchKernels(Particle *particles, int particleCount, SimulationParams params, Plane *planes, int planeCount,
				   float timestep, int *gridCounts, int *gridIndices, vec3 worldMin, vec3 worldMax);