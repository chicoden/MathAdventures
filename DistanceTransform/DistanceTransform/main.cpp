#include <glad/gl.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

unsigned int compileShaderFromFile(GLenum type, const char* path) {
	std::ifstream file(path);
	if (!file.is_open()) return 0;

	std::stringstream buffer;
	buffer << file.rdbuf();
	std::string code = buffer.str();
	const char* source = code.c_str();

	file.close();

	unsigned int shader = glCreateShader(type);
	glShaderSource(shader, 1, &source, NULL);
	glCompileShader(shader);

	int compileStatus;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &compileStatus);
	if (compileStatus == GL_FALSE) {
		int messageLength;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &messageLength);

		char* errorMessage = new char[messageLength];
		glGetShaderInfoLog(shader, messageLength, NULL, errorMessage);

		std::cout << "Failed to compile " << path << '\n' << errorMessage << std::endl;

		delete errorMessage;
		glDeleteShader(shader);
		return 0;
	}

	return shader;
}

unsigned int compileProgram(std::vector<unsigned int>& shaders) {
	std::vector<unsigned int> shaderList = std::move(shaders);
	unsigned int program = glCreateProgram();

	for (unsigned int shader : shaderList)
		glAttachShader(program, shader);

	glLinkProgram(program);

	for (unsigned int shader : shaderList) {
		glDetachShader(program, shader);
		glDeleteShader(shader);
	}

	int linkStatus;
	glGetProgramiv(program, GL_LINK_STATUS, &linkStatus);
	if (linkStatus == GL_FALSE) {
		glDeleteProgram(program);
		return 0;
	}

	return program;
}

int main() {
	if (!glfwInit()) return -1;

	//glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
	GLFWwindow* window = glfwCreateWindow(800, 600, "Distance Transform", NULL, NULL);
	if (!window) {
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);

	int version = gladLoadGL(glfwGetProcAddress);
	if (version == 0) {
		std::cout << "Failed to initialize OpenGL context" << std::endl;
		glfwMakeContextCurrent(NULL);
		glfwTerminate();
		return -1;
	}

	std::cout << "Loaded OpenGL " << GLAD_VERSION_MAJOR(version) << '.' << GLAD_VERSION_MINOR(version) << std::endl;
	std::cout << "Vendor: " << glGetString(GL_VENDOR) << '\n';
	std::cout << "Renderer: " << glGetString(GL_RENDERER) << '\n';
	std::cout << "Version: " << glGetString(GL_VERSION) << '\n';
	std::cout << "GLSL Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;

	float quad[] = {
		-1.0, -1.0,
		+1.0, -1.0,
		+1.0, +1.0,
		-1.0, +1.0
	};

	unsigned int vbo;
	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(quad), quad, GL_STATIC_DRAW);

	std::vector<unsigned int> shaders = {
		compileShaderFromFile(GL_VERTEX_SHADER, "shaders/passthrough.vert"),
		compileShaderFromFile(GL_FRAGMENT_SHADER, "shaders/splash.frag")
	};

	unsigned int program = compileProgram(shaders);

	while (!glfwWindowShouldClose(window)) {
		glClear(GL_COLOR_BUFFER_BIT);

		glUseProgram(program);
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
		glDrawArrays(GL_TRIANGLE_FAN, 0, 4);

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glfwMakeContextCurrent(NULL);
	glfwTerminate();
	return 0;
}