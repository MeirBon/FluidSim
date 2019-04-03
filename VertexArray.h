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