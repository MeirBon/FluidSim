#include "Camera.h"

Camera::Camera(glm::vec3 position, glm::vec3 up, float yaw, float pitch)
	: Front(glm::vec3(0.0f, 0.0f, -1.0f))
	, MovementSpeed(SPEED)
	, MouseSensitivity(SENSITIVITY)
	, Zoom(ZOOM)
{
	Position = position;
	WorldUp = up;
	Yaw = yaw;
	Pitch = pitch;
	updateCameraVectors();
}

Camera::Camera(float posX, float posY, float posZ, float upX, float upY,
			   float upZ, float yaw, float pitch)
	: Front(glm::vec3(0.0f, 0.0f, -1.0f))
	, MovementSpeed(SPEED)
	, MouseSensitivity(SENSITIVITY)
	, Zoom(ZOOM)
{
	Position = { posX, posY, posZ };
	WorldUp = { upX, upY, upZ };
	Yaw = yaw;
	Pitch = pitch;
	updateCameraVectors();
}

glm::mat4 Camera::GetViewMatrix()
{
	return glm::lookAt(Position, Position + Front, Up);
}

glm::mat4 Camera::GetProjectionMatrix(float width, float height, float near, float far) {
	const float fov = glm::radians(Zoom);
	const float aspectRatio = width / height;
	return glm::perspective(fov, aspectRatio, near, far);
}


void Camera::ProcessKeyboard(Camera_Movement direction, float deltaTime, float speed)
{
	float velocity = MovementSpeed * deltaTime * speed;
	if (direction == FORWARD)
		Position += Front * velocity;
	if (direction == BACKWARD)
		Position -= Front * velocity;
	if (direction == LEFT)
		Position -= Right * velocity;
	if (direction == RIGHT)
		Position += Right * velocity;
	if (direction == UP)
		Position += Up * velocity;
	if (direction == DOWN)
		Position -= Up * velocity;
}

void Camera::ProcessMouseMovement(float xoffset, float yoffset, GLboolean constrainPitch)
{
	xoffset *= MouseSensitivity;
	yoffset *= MouseSensitivity;

	Yaw += xoffset;
	Pitch += yoffset;

	if (constrainPitch) {
		if (Pitch > 89.0f)
			Pitch = 89.0f;
		if (Pitch < -89.0f)
			Pitch = -89.0f;
	}

	updateCameraVectors();
}

void Camera::ProcessMouseScroll(int xAxis, int yAxis)
{
	if (Zoom >= 1.0f && Zoom <= 45.0f) {
		Zoom += yAxis * MouseSensitivity;
	} else if (Zoom <= 1.0f) {
		Zoom = 1.0f;
	} else if (Zoom >= 45.0f) {
		Zoom = 45.0f;
	}
}

void Camera::updateCameraVectors()
{
	Front = glm::normalize(glm::vec3(
		cos(glm::radians(Yaw)) * cos(glm::radians(Pitch)),
		sin(glm::radians(Pitch)),
		sin(glm::radians(Yaw)) * cos(glm::radians(Pitch))));
	Right = glm::normalize(glm::cross(Front, WorldUp));
	Up = glm::normalize(glm::cross(Right, Front));
}
