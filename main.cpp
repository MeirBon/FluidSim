#include <GL/glew.h>

#include <GLFW/glfw3.h>
#include <iostream>

#include "Camera.h"
#include "Plane.h"
#include "Shader.h"
#include "Simulator.h"
#include "Timer.h"

using namespace glm;

static bool *keys;
static double lastX, lastY;
static Camera camera = Camera(vec3(8.0f, 10.0f, 30.0f));

static bool runSim = false;

#define PARTICLE_COUNT 5000

void handle_keys(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS)
		keys[key] = true;
	else if (action == GLFW_RELEASE)
		keys[key] = false;
}

#define SCRWIDTH 1024
#define SCRHEIGHT 768

inline void CheckGL(int line)
{
	GLenum err;
	while ((err = glGetError()) != GL_NO_ERROR)
		std::cout << "LINE: " << line << ", OpenGL Error: " << glewGetErrorString(err) << ", " << err << std::endl;
}

#define CHECKGL() CheckGL(__LINE__)

int main(int argc, char *argv[])
{
	Timer timer{};
	if (!glfwInit())
		std::cout << "Could not init GLFW." << std::endl, exit(1);

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	auto *window = glfwCreateWindow(SCRWIDTH, SCRHEIGHT, "Fluid", nullptr, nullptr);
	if (!window)
		std::cout << "Could not init GLFW window." << std::endl, exit(1);

	glfwMakeContextCurrent(window);
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);
	glfwSetKeyCallback(window, handle_keys);

	if (glewInit() != GL_NO_ERROR)
		std::cout << "Could not init GLEW." << std::endl, exit(1);

	keys = new bool[512];
	memset(keys, 0, 512 * sizeof(bool));

	auto *shader = new Shader("shaders/plane.vert", "shaders/plane.frag");

	SimulationParams simulationParams;
	simulationParams.particleRadius = 1.0f;
	simulationParams.smoothingRadius = 1.0f;
	simulationParams.smoothingRadius2 = 1.0f;
	simulationParams.restDensity = 15.0f;
	simulationParams.gravityMult = 2000.0f;
	simulationParams.particleMass = 0.1f;
	simulationParams.particleViscosity = 1.0f;
	simulationParams.particleDrag = 0.025f;
	Simulator simulator(16);
	const auto pid = simulator.addParams(simulationParams);
	simulator.addParticles(PARTICLE_COUNT, pid);
	simulator.addPlane(
		Plane(vec3(0.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f), vec3(0.0f, 0.0f, 1.0f), vec2(100.0f, 100.0f)));

	const auto &particles = simulator.getParticles();
	std::vector<vec3> positions;
	positions.reserve(particles.size());
	for (const auto &p : particles)
	{
		positions.push_back(p.position);
	}

	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
	glEnable(GL_DEPTH_TEST);

	GLuint positionBuffer, VAO;

	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);

	glGenBuffers(1, &positionBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, positionBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * positions.size(), positions.data(), GL_DYNAMIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void *)0);
	glEnableVertexAttribArray(0);

	shader->enable();
	shader->setUniformMatrix4fv("view", glm::mat4(1.0f));
	shader->setUniformMatrix4fv("projection", glm::mat4(1.0f));
	shader->setUniformMatrix4fv("model", glm::mat4(1.0f));
	shader->disable();

	timer.reset();
	float elapsed = 0.1f, elapsedSum = 0.0f;
	while (!glfwWindowShouldClose(window))
	{
		elapsedSum += elapsed;
		if (runSim)
			simulator.update(elapsed);
		positions.clear();
		positions.reserve(particles.size());

		for (const auto &p : particles)
		{
			positions.push_back(p.position);
		}

		glBindVertexArray(VAO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * positions.size(), positions.data(), GL_DYNAMIC_DRAW);
		CHECKGL();

		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void *)0);
		CHECKGL();

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

		shader->enable();
		CHECKGL();

		shader->setUniformMatrix4fv("view", camera.GetViewMatrix());
		shader->setUniformMatrix4fv("projection", camera.GetProjectionMatrix(SCRWIDTH, SCRHEIGHT, 0.1f, 1e34f));
		shader->setUniformMatrix4fv("model", glm::mat4(1.0f));
		shader->setUniform1f("radius", 4.f);
		shader->setUniform3f("color", vec3(1.0f, 0.0f, 0.0f));

		glDrawArrays(GL_POINTS, 0, GLsizei(positions.size()));
		CHECKGL();
		glBindVertexArray(0);

		for (const auto &plane : simulator.getPlanes())
		{
			plane.draw(*shader);
			CHECKGL();
		}

		if (keys[GLFW_KEY_LEFT_SHIFT])
			elapsed *= 5.0f;
		if (keys[GLFW_KEY_Q] || keys[GLFW_KEY_ESCAPE])
			glfwSetWindowShouldClose(window, true);
		if (keys[GLFW_KEY_W])
			camera.ProcessKeyboard(FORWARD, elapsed);
		if (keys[GLFW_KEY_S])
			camera.ProcessKeyboard(BACKWARD, elapsed);
		if (keys[GLFW_KEY_A])
			camera.ProcessKeyboard(LEFT, elapsed);
		if (keys[GLFW_KEY_D])
			camera.ProcessKeyboard(RIGHT, elapsed);
		if (keys[GLFW_KEY_SPACE])
			camera.ProcessKeyboard(UP, elapsed);
		if (keys[GLFW_KEY_LEFT_CONTROL])
			camera.ProcessKeyboard(DOWN, elapsed);
		if (keys[GLFW_KEY_UP])
			camera.ProcessMouseMovement(0.0f, 3.0f);
		if (keys[GLFW_KEY_DOWN])
			camera.ProcessMouseMovement(0.0f, -3.0f);
		if (keys[GLFW_KEY_LEFT])
			camera.ProcessMouseMovement(-3.0f, 0.0f);
		if (keys[GLFW_KEY_RIGHT])
			camera.ProcessMouseMovement(3.0f, 0.0f);
		if (keys[GLFW_KEY_R] && elapsedSum > 100.0f)
			runSim = !runSim, elapsedSum = 0.0f;

		glfwPollEvents();
		glfwSwapBuffers(window);

		while ((elapsed = timer.elapsed()) <= (1000.0f / 60.0f)) // Lock to 60 fps
		{
		}
		timer.reset();
	}

	delete[] keys;
	return 0;
}