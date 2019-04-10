#include <GL/glew.h>

#include <GLFW/glfw3.h>
#include <chrono>
#include <iostream>
#include <thread>

#include "Buffer.h"
#include "Camera.h"
#include "Plane.h"
#include "Shader.h"
#include "Simulator.h"
#include "Timer.h"
#include "VertexArray.h"
#include "Window.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

using namespace glm;

#define PARTICLE_COUNT 30000
#define SCRWIDTH 1024
#define SCRHEIGHT 768
#define DRAW_MESH 1

static bool firstMouse = false, drawMesh = false;
static double lastMouseX, lastMouseY;

inline void CheckGL(int line)
{
	GLenum err;
	while ((err = glGetError()) != GL_NO_ERROR)
		std::cout << "LINE: " << line << ", OpenGL Error: " << glewGetErrorString(err) << ", " << err << std::endl;
}

#define CHECKGL() CheckGL(__LINE__)

// TODO(Dan): Temporarily moved this as a global variable because grid array is to big to fit on stack.
Simulator simulator(16, vec3(-6.0f, 0.0f, 0.0f));
static Camera camera = Camera(vec3(0.0f, 25.0f, 30.0f));
static Window window = Window("FluidSim", SCRWIDTH, SCRHEIGHT);

int main(int argc, char *argv[])
{
	Timer timer{};
	window.setMouseCallback([](double x, double y) {
		if (window.keys[GLFW_MOUSE_BUTTON_RIGHT])
			camera.ProcessMouseMovement(float(x), -float(y));
	});

	bool runSim = false;
	auto shader = Shader("shaders/sphere.vert", "shaders/sphere.frag");
	auto meshShader = Shader("shaders/mesh.vert", "shaders/mesh.frag");
	auto planeShader = Shader("shaders/plane.vert", "shaders/plane.frag");

	SimulationParams params{};
	params.particleRadius = .7f;
	params.smoothingRadius = 1.0f;
	params.smoothingRadiusPow2 = 1.0f;
	params.restDensity = 15.0f;
	params.gravityMult = 2000.0f;
	params.particleMass = 0.1f;
	params.particleViscosity = 1.0f;
	params.particleDrag = 0.025f;

	const auto pid = simulator.addParams(params);
	SimulationParams &simulationParams = simulator.getSimParams()[0];
	simulator.addParticles(PARTICLE_COUNT, pid);

	// Bottom plane
	simulator.addPlane(
		Plane(vec3(0.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f), vec3(0.0f, 0.0f, 1.0f), vec2(20.0f, 20.0f)));

	// Left plane
	simulator.addPlane(
		Plane(vec3(-20.0f, 15.0f, 0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f), vec2(20.0f, 15.0f)));

	// Right plane
	simulator.addPlane(
		Plane(vec3(20.0f, 15.0f, 0.0f), vec3(0.0f, 0.0f, -1.0f), vec3(0.0f, 1.0f, 0.0f), vec2(20.0f, 15.0f)));

	// Far Plane
	simulator.addPlane(
		Plane(vec3(0.0f, 15.0f, -20.0f), vec3(-1.0f, 0.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f), vec2(20.0f, 15.0f)));

	// Near plane
	simulator.addPlane(
		Plane(vec3(0.0f, 15.0f, 20.0f), vec3(1.0f, 0.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f), vec2(20.0f, 15.0f)));

	simulator.setParticleGridBounds();

	const auto &particles = simulator.getParticles();

	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	std::string warn, err;
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, "models/sphere.obj");

	if (!err.empty())
		std::cout << "Model error: " << err << std::endl, exit(1);
	if (!ret)
		exit(1);

	std::vector<vec3> vertices, normals;
	for (auto &s : shapes)
	{
		int index_offset = 0;
		for (unsigned int f = 0; f < s.mesh.num_face_vertices.size(); f++)
		{
			int fv = s.mesh.num_face_vertices[f];
			for (int v = 0; v < fv; v++)
			{
				const auto &idx = s.mesh.indices[index_offset + v];
				const auto &vx = attrib.vertices[3 * idx.vertex_index + 0];
				const auto &vy = attrib.vertices[3 * idx.vertex_index + 1];
				const auto &vz = attrib.vertices[3 * idx.vertex_index + 2];
				if (!attrib.normals.empty())
				{
					const auto &nx = attrib.normals[3 * idx.normal_index + 0];
					const auto &ny = attrib.normals[3 * idx.normal_index + 1];
					const auto &nz = attrib.normals[3 * idx.normal_index + 2];

					normals.emplace_back(nx, ny, nz);
				}

				const vec3 vertex = vec3(vx, vy, vz);

				vertices.push_back(glm::normalize(vertex));
			}

			index_offset += fv;
		}
	}

	auto vBuffer = Buffer(GL_ARRAY_BUFFER, vertices.size(), sizeof(vec3), vertices.data(), 3);
	auto nBuffer = Buffer(GL_ARRAY_BUFFER, normals.size(), sizeof(vec3), normals.data(), 3);
	auto sphereVAO = VertexArray();
	sphereVAO.assignBuffer(0, vBuffer);
	sphereVAO.assignBuffer(1, nBuffer);

	shader.enable();
	shader.setUniformMatrix4fv("view", glm::mat4(1.0f));
	shader.setUniformMatrix4fv("projection", glm::mat4(1.0f));
	shader.setUniformMatrix4fv("model", glm::mat4(1.0f));
	shader.setUniform3f("lightIntensity", vec3(1.0f));
	shader.setUniform3f("lightDirection", glm::normalize(vec3(-1.0f, 1.0f, 0.0f)));
	shader.setUniform3f("ambient", vec3(0.1f));
	shader.disable();

	meshShader.enable();
	meshShader.setUniformMatrix4fv("view", glm::mat4(1.0f));
	meshShader.setUniformMatrix4fv("projection", glm::mat4(1.0f));
	meshShader.setUniformMatrix4fv("model", glm::mat4(1.0f));
	meshShader.setUniform3f("lightIntensity", vec3(1.0f));
	meshShader.setUniform3f("lightDirection", glm::normalize(vec3(-1.0f, 1.0f, 0.0f)));
	meshShader.setUniform3f("ambient", vec3(0.1f));
	meshShader.disable();

	planeShader.enable();
	planeShader.setUniformMatrix4fv("view", glm::mat4(1.0f));
	planeShader.setUniformMatrix4fv("projection", glm::mat4(1.0f));
	planeShader.setUniformMatrix4fv("model", glm::mat4(1.0f));
	planeShader.setUniform3f("lightIntensity", vec3(1.0f));
	planeShader.setUniform3f("lightDirection", glm::normalize(vec3(-1.0f, -1.0f, 0.0f)));
	planeShader.setUniform3f("ambient", vec3(0.1f));
	planeShader.disable();

	timer.reset();
	// TODO(Dan): Not obvious what elapsedSum does (Meir): It's used when a button is pressed, makes sure you don't keep
	// switching toggles
	float elapsed = 0.1f, elapsedSum = 0.0f;
	simulator.update(0.0f);
	while (!window.shouldClose())
	{
		const auto size = window.getSize();
		const int width = std::get<0>(size);
		const int height = std::get<1>(size);

		elapsedSum += elapsed;
		if (runSim)
			simulator.update(elapsed);

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

		if (drawMesh)
		{
			meshShader.enable();
			meshShader.setUniformMatrix4fv("view", camera.GetViewMatrix());
			meshShader.setUniformMatrix4fv("projection", camera.GetProjectionMatrix(width, height, 0.1f, 1e34f));
			meshShader.setUniform1f("radius", simulationParams.particleRadius);
			meshShader.setUniform3f("color", vec3(0.40f, 0.75f, 1.0f));

			glm::mat4 model(1.0f);
			const float scale = 1.0f / simulator.getVoxelScale() * 2.0f;
			const auto worldMin = simulator.getWorldMin();
			model = glm::translate(model, {worldMin.x, 0.0f, worldMin.z});
			model = glm::scale(model, vec3(scale));
			meshShader.setUniformMatrix4fv("model", model);

			simulator.extractSurface(meshShader);
			meshShader.disable();
			CHECKGL();
		}
		else
		{
			shader.enable();
			shader.setUniformMatrix4fv("view", camera.GetViewMatrix());
			shader.setUniformMatrix4fv("projection", camera.GetProjectionMatrix(width, height, 0.1f, 1e34f));
			shader.setUniform1f("radius", simulationParams.particleRadius);
			shader.setUniform3f("color", vec3(0.40f, 0.75f, 1.0f));

			sphereVAO.bind();
			for (const auto &p : particles)
			{
				shader.setUniform1f("pressure", p.pressure);
				shader.setUniform3f("position", p.position);
				glDrawArrays(GL_TRIANGLES, 0, vertices.size());
			}
			CHECKGL();
			sphereVAO.unbind();
		}

		planeShader.enable();
		planeShader.setUniformMatrix4fv("view", camera.GetViewMatrix());
		planeShader.setUniformMatrix4fv("projection", camera.GetProjectionMatrix(width, height, 0.1f, 1e34f));

		for (const auto &plane : simulator.getPlanes())
		{
			plane.draw(planeShader);
			CHECKGL();
		}

		const auto *keys = window.keys;

		elapsed *= 3.0f;
		if (keys[GLFW_KEY_LEFT_SHIFT])
			elapsed *= 5.0f;
		if (keys[GLFW_KEY_Q] || keys[GLFW_KEY_ESCAPE])
			window.close();
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
			simulator.moveBounds(glm::vec3(0, 0.10, 0));
		if (keys[GLFW_KEY_DOWN])
			simulator.moveBounds(glm::vec3(0, -0.10, 0));
		if (keys[GLFW_KEY_LEFT])
			simulator.moveBounds(glm::vec3(-0.1, 0, 0));
		if (keys[GLFW_KEY_RIGHT])
			simulator.moveBounds(glm::vec3(0.1, 0, 0));
		if (keys[GLFW_KEY_R] && elapsedSum > 200.0f)
			runSim = !runSim, elapsedSum = 0.0f;
		if (keys[GLFW_KEY_BACKSPACE] && elapsedSum > 200.0f)
			simulator.reset(), elapsedSum = 0.0f;

		ImGui::Begin("Params");
		ImGui::Text("Running: %i", runSim);
		ImGui::Checkbox("Render Surface", &drawMesh);
		ImGui::DragFloat("Mass", &simulationParams.particleMass,0.05,0.05,1);
		ImGui::DragFloat("Radius", &simulationParams.particleRadius, 0.005f, 0.7f, 10.0f);
		ImGui::DragFloat("Collision Drag", &simulationParams.particleDrag, 0.005f, 0.0f, 1.0f);
		ImGui::DragFloat("Viscosity", &simulationParams.particleViscosity, 0.005f, 0.0f, 10.0f);
		ImGui::DragFloat("Gravity", &simulationParams.gravity.y, 0.01f, 0.0f, 100.0f);
		if (ImGui::Button("reset"))
			simulationParams = params, simulator.reset();

		ImGui::End();

		window.pollEvents();
		window.present();

		constexpr float desiredFrametime = 1000.0f / 60.0f;
		elapsed = timer.elapsed();
		if (elapsed != desiredFrametime)
		{
			std::this_thread::sleep_for(std::chrono::milliseconds(int(desiredFrametime - elapsed)));
			elapsed = timer.elapsed();
		}
		timer.reset();
	}

	return 0;
}