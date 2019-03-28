#pragma once

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <string>

class Shader {
private:
    GLuint m_ShaderId;
    const char* m_VertPath;
    const char* m_FragPath;
    const char* m_GeoPath;

    GLuint load();
    void checkCompileErrors(GLuint shader, std::string type);

public:
    Shader(const char* vertexPath, const char* fragmentPath,
        const char* geometryPath = nullptr);
    ~Shader();

    const GLuint& getShaderId() const;

    void enable() const;
    void disable() const;

    GLint getUniformLocation(const char* name);
    void setUniform1f(const char* name, float value);
    void setUniform1i(const char* name, int value);
    void setUniform1b(const char* name, bool value);
    void setUniform2f(const char* name, glm::vec2 value);
    void setUniform2f(const char* name, float x, float y);
    void setUniform3f(const char* name, glm::vec3 value);
    void setUniform3f(const char* name, float x, float y, float z);
    void setUniform4f(const char* name, glm::vec4 value);
    void setUniform4f(const char* name, float x, float y, float z, float w);
    void setUniformMatrix4fv(const char* name, glm::mat4 value);

    inline GLint getAttributeLocation(const char* name)
    {
        return glGetAttribLocation(m_ShaderId, name);
    }
};