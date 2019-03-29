#pragma once

#include <iostream>
#include <tuple>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

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
	std::pair<int, int> getSize();

	bool *keys;

  private:
	GLFWwindow *m_Window;
};

static Window *instance = nullptr;

void handle_keys(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS)
		instance->keys[key] = true;
	else if (action == GLFW_RELEASE)
		instance->keys[key] = false;
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

	if (glewInit() != GL_NO_ERROR)
		std::cout << "Could not init GLEW." << std::endl, exit(1);

	keys = new bool[512];
	memset(keys, 0, 512 * sizeof(bool));
}

Window::~Window()
{
	delete[] keys;
	glfwDestroyWindow(m_Window);
	glfwTerminate();
}

void Window::pollEvents() { glfwPollEvents(); }
bool Window::shouldClose() { return glfwWindowShouldClose(m_Window) == 1; }
void Window::present() { glfwSwapBuffers(m_Window); }
void Window::setTitle(const char *title) { glfwSetWindowTitle(m_Window, title); }
void Window::close() { glfwSetWindowShouldClose(m_Window, 1); }

std::pair<int, int> Window::getSize()
{
	int width, height;
	glfwGetWindowSize(m_Window, &width, &height);
	return std::make_pair(width, height);
}