#pragma once

#include <glm/glm.hpp>

struct SimulationParams;

using namespace glm;

struct Particle {
    vec3 position; // 12
    vec3 velocity; // 24
    vec3 forcePhysic; // 36
    vec3 forceHeading; // 48

    float density; // 52
    float pressure; // 56
    int parameterID; // 60
    float dummy; // 64 Cache-aligned :D

    Particle(vec3 pos, int pId);
};