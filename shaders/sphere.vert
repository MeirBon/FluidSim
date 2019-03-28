#version 330 core

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

attribute vec3 position;
attribute float velocity;
attribute float radius;

varying float speed;

void main(){
    gl_PointSize = radius;
    speed = velocity;

    gl_Position = MVP * vec4(position, 1.0);
}