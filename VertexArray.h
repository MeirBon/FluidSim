#pragma once

#include <GL/glew.h>

#include "Buffer.h"

class VertexArray
{
  public:
	VertexArray();
	~VertexArray();

	void bind();
	static void unbind();

	void assignBuffer(GLuint location, Buffer &buffer);

  private:
	GLuint m_ID;
};

VertexArray::VertexArray() { glGenVertexArrays(1, &m_ID); }
VertexArray::~VertexArray() { glDeleteVertexArrays(1, &m_ID); }

void VertexArray::bind() { glBindVertexArray(m_ID); }
void VertexArray::unbind() { glBindVertexArray(0); }

void VertexArray::assignBuffer(GLuint location, Buffer &buffer)
{
	bind();
	glEnableVertexAttribArray(location);
	buffer.bind();
	glVertexAttribPointer(location, buffer.m_Components, buffer.m_DataType, GL_FALSE, 0, (void *)0);
	unbind();
}