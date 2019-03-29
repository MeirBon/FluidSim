#version 410 core

out vec4 Color;

in vec3 Normal;

uniform vec3 color;
uniform vec3 ambient;

uniform vec3 lightIntensity;
uniform vec3 lightDirection;

void main()
{
    Color = vec4(color * max(0.0f, dot(-lightDirection, Normal)) + ambient * color, 0.3f);
}