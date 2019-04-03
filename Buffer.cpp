#include "Buffer.h"

Buffer::Buffer(unsigned int type, size_t count, size_t bytesPerElement, void *data, GLuint numberOfComponents,
			   GLuint dataType, GLuint usage)
	: m_Type(type), m_DataType(dataType), m_Components(numberOfComponents)
{
	glGenBuffers(1, &m_ID);
	bind();
	glBufferData(type, count * bytesPerElement, data, usage);
	unbind();
}
Buffer::~Buffer() { glDeleteBuffers(1, &m_ID); }

void Buffer::bind() { glBindBuffer(m_Type, m_ID); }
void Buffer::unbind() { glBindBuffer(m_Type, 0); }