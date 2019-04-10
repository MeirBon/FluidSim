#include "Plane.h"

using namespace glm;

Plane::Plane(vec3 pos, vec3 r, vec3 f, vec2 dims) : position(pos), right(r), forward(f), size(dims)
{
	std::vector<vec3> vertices;
	std::vector<vec3> normals;
	const vec3 offsetRight = r * dims.x;
	const vec3 offsetForward = f * dims.y;
	const vec3 normal = glm::normalize(glm::cross(r, f));

	vertices.push_back(vec3(0.0f) - offsetRight - offsetForward);
	vertices.push_back(vec3(0.0f) + offsetRight - offsetForward);
	vertices.push_back(vec3(0.0f) - offsetRight + offsetForward);
	vertices.push_back(vec3(0.0f) + offsetRight - offsetForward);
	vertices.push_back(vec3(0.0f) - offsetRight + offsetForward);
	vertices.push_back(vec3(0.0f) + offsetRight + offsetForward);

	for (int i = 0; i < 6; i++)
		normals.push_back(normal);

	// Load data onto GPU
	glGenBuffers(1, &VBOV);
	glBindBuffer(GL_ARRAY_BUFFER, VBOV);
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * vertices.size(), vertices.data(), GL_STATIC_DRAW);
	glGenBuffers(1, &VBON);
	glBindBuffer(GL_ARRAY_BUFFER, VBON);
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * normals.size(), normals.data(), GL_STATIC_DRAW);

	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);

	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, VBOV);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void *)0); // Vertices

	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, VBON);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void *)0); // Normals

	glBindVertexArray(0);
	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
}

void Plane::draw(Shader &shader) const
{
	shader.enable();
	shader.setUniform3f("color", vec3(0.9f));
	shader.setUniform3f("translation", position);
	glBindVertexArray(VAO);
	glDrawArrays(GL_TRIANGLES, 0, 6);
	glBindVertexArray(0);
	shader.disable();
}

void Plane::translate(const vec3 &pos) { position += pos; }
