#version 410 core

out vec4 Color;

in vec3 Normal;
in vec3 Pos;

uniform vec3 color;
uniform vec3 ambient;
uniform float pressure;

uniform vec3 lightIntensity;
uniform vec3 lightDirection;

void main()
{
    // float colorFactor = min(1.0f, max(1.0f - (pressure - 3000.0f) / 15000.0f, 0.01f));
    float colorFactor = 1.0f;
    vec3 col = vec3(1.0f);

    const float Ka = 1.0f;
    const float Kd = 1.0f;
    const float Ks = 0.5f;

    vec3 N = normalize(Normal);
    vec3 L = normalize(lightDirection);

    float lambertian = max(dot(N, L), 0.0f);
    float specular = 0.0f;

    if (lambertian > 0.0f)
    {
        vec3 R = reflect(L, N);
        vec3 V = normalize(-Pos);

        float specularAngle = max(dot(R, V), 0.0f);
        specular = pow(specularAngle, 16.0f);
    }

    col = Ka * ambient + Kd * lambertian * color + Ks * specular * vec3(1.0f);

    Color = vec4(col * colorFactor, 1.0f);
}