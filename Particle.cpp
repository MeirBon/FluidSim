#include "Particle.h"

Particle::Particle(vec3 pos, int pId)
    : position(pos)
    , parameterID(pId)
{
    auto zero = vec3(0.0f);
    velocity = zero;
    forcePhysic = zero;
    forceHeading = zero;
    density = 0.0f;
    pressure = 0.0f;
}
