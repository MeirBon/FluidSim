#include "Shader.h"

#include <fstream>
#include <iostream>
#include <istream>
#include <sstream>
#include <string>

std::string readFile(const char *filepath)
{
	std::string buffer;
	std::ifstream file;

	file.exceptions(std::ifstream::failbit | std::ifstream::badbit);

	try
	{
		file.open(filepath);
		std::stringstream fileStream;

		fileStream << file.rdbuf();

		file.close();

		buffer.append(fileStream.str());
	}
	catch (const std::ifstream::failure &e)
	{
		std::cout << "Could not read file: " << filepath << std::endl;
	}

	return buffer;
}

Shader::Shader(const char *vertexPath, const char *fragmentPath, const char *geometryPath)
	: m_VertPath(vertexPath), m_FragPath(fragmentPath), m_GeoPath(geometryPath)
{
	m_ShaderId = load();
}

GLuint Shader::load()
{
	const GLuint program = glCreateProgram();
	const GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
	const GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);

	std::string vert_string = readFile(m_VertPath);
	std::string frag_string = readFile(m_FragPath);

	const char *vert_source = vert_string.c_str();
	const char *frag_source = frag_string.c_str();

	glShaderSource(vert_shader, 1, &vert_source, nullptr);
	glCompileShader(vert_shader);
	checkCompileErrors(vert_shader, "VERTEX");

	glShaderSource(frag_shader, 1, &frag_source, nullptr);
	glCompileShader(frag_shader);
	checkCompileErrors(frag_shader, "FRAGMENT");

	glAttachShader(program, vert_shader);
	glAttachShader(program, frag_shader);

	if (m_GeoPath != nullptr)
	{
		const GLuint geo_shader = glCreateShader(GL_GEOMETRY_SHADER);
		std::string geo_string = readFile(m_GeoPath);
		const char *geo_source = geo_string.c_str();

		glShaderSource(geo_shader, 1, &geo_source, nullptr);
		glCompileShader(geo_shader);
		checkCompileErrors(geo_shader, "GEOMETRY");
		glAttachShader(program, geo_shader);
	}

	glLinkProgram(program);
	glValidateProgram(program);
	// checkCompileErrors(program, "PROGRAM");

	glDeleteShader(vert_shader);
	glDeleteShader(frag_shader);

	return program;
}

void Shader::checkCompileErrors(GLuint shader, std::string type)
{
	GLint success;
	GLchar infoLog[1024];
	if (type != "PROGRAM")
	{
		glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
		if (!success)
		{
			glGetShaderInfoLog(shader, 1024, NULL, infoLog);
			std::cout << "ERROR::SHADER_COMPILATION_ERROR of type: " << type << "\n"
					  << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
		}
	}
	else
	{
		glGetProgramiv(shader, GL_LINK_STATUS, &success);
		if (!success)
		{
			glGetProgramInfoLog(shader, 1024, NULL, infoLog);
			std::cout << "ERROR::PROGRAM_LINKING_ERROR of type: " << type << "\n"
					  << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
		}
	}
}

void Shader::enable() const { glUseProgram(m_ShaderId); }

void Shader::disable() const { glUseProgram(0); }

Shader::~Shader() { glDeleteProgram(m_ShaderId); }

const GLuint &Shader::getShaderId() const { return m_ShaderId; }

void Shader::setUniform2f(const char *name, glm::vec2 value)
{
	glUniform2f(getUniformLocation(name), value.x, value.y);
}

void Shader::setUniform2f(const char *name, float x, float y) { glUniform2f(getUniformLocation(name), x, y); }

void Shader::setUniform3f(const char *name, glm::vec3 value)
{
	glUniform3f(getUniformLocation(name), value.x, value.y, value.z);
}

void Shader::setUniform3f(const char *name, float x, float y, float z)
{
	glUniform3f(getUniformLocation(name), x, y, z);
}

void Shader::setUniform4f(const char *name, glm::vec4 value)
{
	glUniform4f(getUniformLocation(name), value.x, value.y, value.z, value.w);
}

void Shader::setUniform4f(const char *name, float x, float y, float z, float w)
{
	glUniform4f(getUniformLocation(name), x, y, z, w);
}

GLint Shader::getUniformLocation(const char *name) { return glGetUniformLocation(m_ShaderId, name); }

void Shader::setUniform1f(const char *name, float value) { glUniform1f(getUniformLocation(name), value); }

void Shader::setUniform1i(const char *name, int value) { glUniform1i(getUniformLocation(name), value); }

void Shader::setUniform1b(const char *name, bool value) { glUniform1i(getUniformLocation(name), (int)value); }

void Shader::setUniformMatrix4fv(const char *name, glm::mat4 value)
{
	glUniformMatrix4fv(getUniformLocation(name), 1, GL_FALSE, glm::value_ptr(value));
}