#pragma once

#include <iostream>
#include <tuple>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <ImGuiGLFW.h>
#include <ImGuiOpenGL3.h>

#include <functional>

class Window
{
  public:
	Window(const char *title, int width, int height);
	~Window();

	void pollEvents();
	bool shouldClose();
	void present();
	void setTitle(const char *title);
	void close();

	void setMouseCallback(std::function<void(double, double)> callback);
	void setResizeCallback(std::function<void(int, int)> callback);
	std::pair<int, int> getSize();

  public:
	bool *keys;
	GLFWwindow *m_Window;
	ImGuiContext *m_Context;

	std::function<void(double, double)> onMouseMove = [](double, double) {};
	std::function<void(int, int)> onResize = [](int, int) {};
};

static Window *instance = nullptr;

void handle_keys(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS)
		instance->keys[key] = true;
	else if (action == GLFW_RELEASE)
		instance->keys[key] = false;
}

void mouse_button_callback(GLFWwindow *window, int button, int action, int mods)
{
	if (action == GLFW_PRESS)
		instance->keys[button] = true;
	else if (action == GLFW_RELEASE)
		instance->keys[button] = false;
}

void window_size_callback(GLFWwindow *window, int width, int height)
{
	instance->onResize(width, height);
	glViewport(0, 0, width, height);
}

Window::Window(const char *title, int width, int height)
{
	instance = this;
	if (!glfwInit())
		std::cout << "Could not init GLFW." << std::endl, exit(1);

	glfwSetErrorCallback([](int error, const char *description) {
		std::cout << "GLFW error (" << error << "): " << description << std::endl;
	});

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

	m_Window = glfwCreateWindow(width, height, "Fluid", nullptr, nullptr);
	if (!m_Window)
		std::cout << "Could not init GLFW window." << std::endl, exit(1);

	glfwMakeContextCurrent(m_Window);
	glfwSetInputMode(m_Window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);
	glfwSetKeyCallback(m_Window, handle_keys);
	glfwSetMouseButtonCallback(m_Window, mouse_button_callback);
	glfwSetWindowSizeCallback(m_Window, window_size_callback);
	glfwSetCursorPosCallback(m_Window, [](auto *window, auto x, auto y) {
		static bool first = true;
		static double lastX, lastY;
		if (first)
			first = false, lastX = x, lastY = y;

		instance->onMouseMove(x - lastX, y - lastY);
		lastX = x, lastY = y;
	});

	if (glewInit() != GL_NO_ERROR)
		std::cout << "Could not init GLEW." << std::endl, exit(1);

	keys = new bool[512];
	memset(keys, 0, 512 * sizeof(bool));

	IMGUI_CHECKVERSION();
	m_Context = ImGui::CreateContext();
	ImGui::StyleColorsDark();

	ImGui_ImplGlfw_InitForOpenGL(m_Window, true);
	ImGui_ImplOpenGL3_Init("#version 150");

	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();
}

Window::~Window()
{
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext(m_Context);
	delete[] keys;
	glfwDestroyWindow(m_Window);
	glfwTerminate();
}

void Window::pollEvents() { glfwPollEvents(); }
bool Window::shouldClose() { return glfwWindowShouldClose(m_Window) == 1; }
void Window::present()
{
	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

	glfwSwapBuffers(m_Window);

	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();
}
void Window::setTitle(const char *title) { glfwSetWindowTitle(m_Window, title); }
void Window::close() { glfwSetWindowShouldClose(m_Window, 1); }

inline void Window::setMouseCallback(std::function<void(double, double)> callback) { onMouseMove = callback; }

inline void Window::setResizeCallback(std::function<void(int, int)> callback) { onResize = callback; }

std::pair<int, int> Window::getSize()
{
	int width, height;
	glfwGetWindowSize(m_Window, &width, &height);
	return std::make_pair(width, height);
}