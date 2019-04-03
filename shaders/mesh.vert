#version 410 core

layout(location = 0) in vec3 pos;
layout(location = 1) in vec3 normal;
layout(location = 2) in float material;

out vec3 Normal;
out vec3 Pos;

uniform mat4 view; // Matrix for camera view direction
uniform mat4 projection; // Matrix to project world onto screen
uniform mat4 model; // Matrix which translates model
uniform vec3 correction;

uniform float radius;

void main()
{
    Normal = normal;
    Pos = pos;

    gl_PointSize = radius;
    gl_Position = projection * view * model * vec4(pos / vec3(2.0f, 3.0f / 2.0f, 2.0f), 1.0f);
}