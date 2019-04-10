#pragma once

#include <glm/ext/matrix_transform.hpp>
#include <glm/glm.hpp>
#include <vector>

#include <GL/glew.h>

#include "Shader.h"

using namespace glm;

struct Plane
{
	vec3 position;
	vec3 right;
	vec3 forward;
	vec2 size;

	GLuint VBOV, VBON, VAO;

	explicit Plane(vec3 pos, vec3 r, vec3 f, vec2 dims);
	void draw(Shader &shader) const;
	void translate(const vec3 &pos);
};
