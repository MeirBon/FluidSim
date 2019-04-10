#version 410 core

layout(location = 0) in vec3 pos;
layout(location = 1) in vec3 normal;

uniform mat4 view; // Matrix for camera view direction
uniform mat4 projection; // Matrix to project world onto screen
uniform mat4 model; // Matrix which translates model

uniform vec3 translation;

out vec3 Normal;

void main()
{
	Normal = normal;
	gl_Position = projection * view * model * vec4(pos + translation, 1.0f);
}