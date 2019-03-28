#version 330 core

layout(location = 0) in vec3 pos;
layout(location = 1) in vec3 normal;

out vec3 Normal;

uniform mat4 view; // Matrix for camera view direction
uniform mat4 projection; // Matrix to project world onto screen
uniform mat4 model; // Matrix which translates model
uniform vec3 position;
uniform float radius;

void main()
{
    Normal = normal;
    gl_PointSize = radius;
    gl_Position = projection * view * model * vec4(pos * (radius * 0.5f) + position, 1.0f);
}