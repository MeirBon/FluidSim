#include "Plane.h"

using namespace glm;

Plane::Plane(vec3 pos, vec3 r, vec3 f, vec2 dims) : position(pos), right(r), forward(f), size(dims)
{
	std::vector<vec3> positions;
	const vec3 offsetRight = r * (dims.x / 2.0f);
	const vec3 offsetForward = f * (dims.y / 2.0f);

	positions.push_back(pos - offsetRight - offsetForward);
	positions.push_back(pos + offsetRight - offsetForward);
	positions.push_back(pos - offsetRight + offsetForward);
	positions.push_back(pos + offsetRight - offsetForward);
	positions.push_back(pos - offsetRight + offsetForward);
	positions.push_back(pos + offsetRight + offsetForward);

	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);

	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * positions.size(), positions.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void *)0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glEnableVertexAttribArray(0);
}

void Plane::draw(Shader &shader) const
{
	shader.enable();
	shader.setUniform1f("radius", 1.0f);
	shader.setUniform3f("color", vec3(0.0f, 1.0f, 0.0f));
	glBindVertexArray(VAO);
	glDrawArrays(GL_TRIANGLES, 0, 6);
	glBindVertexArray(0);
}