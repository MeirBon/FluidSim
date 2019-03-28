#version 330 core

out vec3 Color;

in vec3 Normal;
uniform vec3 color;
uniform vec3 ambient;

uniform vec3 lightIntensity;
uniform vec3 lightDirection;

void main()
{
    Color = color * max(0.0f, dot(-lightDirection, Normal)) + ambient;
}