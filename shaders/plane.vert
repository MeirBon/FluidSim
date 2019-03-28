#version 330 core

layout(location = 0) in vec3 pos;

uniform mat4 view; // Matrix for camera view direction
uniform mat4 projection; // Matrix to project world onto screen
uniform mat4 model; // Matrix which translates model
uniform float radius;

void main()
{
    gl_PointSize = radius;
    gl_Position = projection * view * model * vec4(pos, 1.0f);
}