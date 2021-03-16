// Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>
#include <sstream>

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#include <GLFW/glfw3.h>

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
using namespace glm;

// Include AntTweakBar
#include <AntTweakBar.h>

#include <common/shader.hpp>
#include <common/controls.hpp>
#include <common/objloader.hpp>
#include <common/vboindexer.hpp>

// ATTN 1A is the general place in the program where you have to change the code base to satisfy a Task of Project 1A.
// ATTN 1B for Project 1B. ATTN 1C for Project 1C. Focus on the ones relevant for the assignment you're working on.

typedef struct Vertex {
	float Position[4];
	float Color[4];
	void SetCoords(float *coords) {
		Position[0] = coords[0];
		Position[1] = coords[1];
		Position[2] = coords[2];
		Position[3] = coords[3];
	}
	void SetColor(float *color) {
		Color[0] = color[0];
		Color[1] = color[1];
		Color[2] = color[2];
		Color[3] = color[3];
	}
};

// ATTN: use POINT structs for cleaner code (POINT is a part of a vertex)
// allows for (1-t)*P_1+t*P_2  avoiding repeat for each coordinate (x,y,z)
typedef struct point {
	float x, y, z;
	point(const float x = 0, const float y = 0, const float z = 0) : x(x), y(y), z(z){};
	point(float *coords) : x(coords[0]), y(coords[1]), z(coords[2]){};
	point operator -(const point& a) const {
		return point(x - a.x, y - a.y, z - a.z);
	}
	point operator +(const point& a) const {
		return point(x + a.x, y + a.y, z + a.z);
	}
	point operator *(const float& a) const {
		return point(x * a, y * a, z * a);
	}
	point operator /(const float& a) const {
		return point(x / a, y / a, z / a);
	}
	float* toArray() {
		float array[] = { x, y, z, 1.0f };
		return array;
	}
};

// Function prototypes
int initWindow(void);
void initOpenGL(void);
void createVAOs(Vertex[], GLushort[], int);
void createObjects(void);
void pickVertex(bool method); //changed to store color values
void storeColorFunct(GLuint gPickedIndex);
void moveVertex(void);
void key_1(void);
void key_2(void);
void key_4(void);
void key_5(void);
void renderScene(void);
void cleanup(void);
static void mouseCallback(GLFWwindow*, int, int, int);
static void keyCallback(GLFWwindow*, int, int, int, int);

// GLOBAL VARIABLES
GLFWwindow* window;
const GLuint window_width = 1024, window_height = 768;

glm::mat4 gProjectionMatrix;
glm::mat4 gViewMatrix;

// Program IDs
GLuint programID;
GLuint pickingProgramID;

// Uniform IDs
GLuint MatrixID;
GLuint ViewMatrixID;
GLuint ModelMatrixID;
GLuint PickingMatrixID;
GLuint pickingColorArrayID;
GLuint pickingColorID;

GLuint gPickedIndex;
std::string gMessage;

// ATTN: INCREASE THIS NUMBER AS YOU CREATE NEW OBJECTS
const GLuint NumObjects = 9; // Number of objects types in the scene

// Keeps track of IDs associated with each object
GLuint VertexArrayId[NumObjects];
GLuint VertexBufferId[NumObjects];
GLuint IndexBufferId[NumObjects];

size_t VertexBufferSize[NumObjects];
size_t IndexBufferSize[NumObjects];
size_t NumVerts[NumObjects];	// Useful for glDrawArrays command
size_t NumIdcs[NumObjects];	// Useful for glDrawElements command

// Initialize ---  global objects -- not elegant but ok for this project
const size_t IndexCount = 8;
const int maxsub = 5;
const int MAXSIZE = 320;
const int BIG_N = 256;
Vertex Vertices[IndexCount];
Vertex WhiteVertices[IndexCount];
Vertex SubVertices[maxsub+1][MAXSIZE];
Vertex MidVertices[IndexCount];
Vertex BBVertices[IndexCount * 2];
Vertex travelVert[1];
Vertex tangent[2];
Vertex normal[2];
Vertex binormal[2];
Vertex BBVert_Line[IndexCount * (BIG_N+1)];
GLushort Indices[IndexCount];
GLushort SubIndices[MAXSIZE];
GLushort BBIndices[IndexCount * BIG_N];

float storeColor[4] = { 0.0f, 0.0f, 0.0f, 1.0f };
int storeIndex = 0;
float store_xpos = 0.0f;
bool shift_press = false;
bool startup = true;
bool animation = false;
bool double_view = false;
float timer = 0.0f;
const float PI = 3.14159265358979323846;
int subdiv_k = 0;
bool bezier_render = true;
int sz = (int)IndexCount * pow(2, subdiv_k);

// ATTN: DON'T FORGET TO INCREASE THE ARRAY SIZE IN THE PICKING VERTEX SHADER WHEN YOU ADD MORE PICKING COLORS
// NOTE: currently set to 8 in vertex shader for the 8 points
float pickingColor[IndexCount];

int initWindow(void) {
	// Initialise GLFW
	if (!glfwInit()) {
		fprintf(stderr, "Failed to initialize GLFW\n");
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // FOR MAC

	// ATTN: Project 1A, Task 0 == Change the name of the window
	// Open a window and create its OpenGL context
	window = glfwCreateWindow(window_width, window_height, "Fohrman,Kyle(0514-1178)", NULL, NULL);
	if (window == NULL) {
		fprintf(stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n");
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		return -1;
	}

	// Initialize the GUI display
	TwInit(TW_OPENGL_CORE, NULL);
	TwWindowSize(window_width, window_height);
	TwBar * GUI = TwNewBar("Picking");
	TwSetParam(GUI, NULL, "refresh", TW_PARAM_CSTRING, 1, "0.1");
	TwAddVarRW(GUI, "Last picked object", TW_TYPE_STDSTRING, &gMessage, NULL);

	// Set up inputs
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_FALSE);
	glfwSetCursorPos(window, window_width / 2, window_height / 2);
	glfwSetMouseButtonCallback(window, mouseCallback);
	glfwSetKeyCallback(window, keyCallback);

	return 0;
}

void initOpenGL(void) {
	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);
	// Cull triangles which normal is not towards the camera
	glEnable(GL_CULL_FACE);

	// Projection matrix : 45� Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
	//glm::mat4 ProjectionMatrix = glm::perspective(45.0f, 4.0f / 3.0f, 0.1f, 100.0f);
	// Or, for Project 1, use an ortho camera :
	gProjectionMatrix = glm::ortho(-4.0f, 4.0f, -3.0f, 3.0f, 0.0f, 100.0f); // In world coordinates

	// Camera matrix
	gViewMatrix = glm::lookAt(
		glm::vec3(0, 0, -5), // Camera is at (0,0,-5) below the origin, in World Space
		glm::vec3(0, 0, 0), // and looks at the origin
		glm::vec3(0, 1, 0)  // Head is looking up at the origin (set to 0,-1,0 to look upside-down)
	);

	// Create and compile our GLSL program from the shaders
	programID = LoadShaders("p1_StandardShading.vertexshader", "p1_StandardShading.fragmentshader");
	pickingProgramID = LoadShaders("p1_Picking.vertexshader", "p1_Picking.fragmentshader");

	// Get a handle for our "MVP" uniform
	MatrixID = glGetUniformLocation(programID, "MVP");
	ViewMatrixID = glGetUniformLocation(programID, "V");
	ModelMatrixID = glGetUniformLocation(programID, "M");
	PickingMatrixID = glGetUniformLocation(pickingProgramID, "MVP");
	
	// Get a handle for our "pickingColorID" uniform
	pickingColorArrayID = glGetUniformLocation(pickingProgramID, "PickingColorArray");
	pickingColorID = glGetUniformLocation(pickingProgramID, "PickingColor");

	// Define pickingColor array for picking program
	// use a for-loop here
	for (int i = 0; i < IndexCount; i++)
	{
		pickingColor[i] = i / 255.0f;
	}

	// Define objects
	createObjects();

	// ATTN: create VAOs for each of the newly created objects here:
	// for several objects of the same type use a for-loop
	int obj = 0;  // initially there is only one type of object 
	VertexBufferSize[obj] = sizeof(Vertices);
	VertexBufferSize[1] = sizeof(WhiteVertices);
	VertexBufferSize[2] = sizeof(SubVertices);
	VertexBufferSize[3] = sizeof(MidVertices);
	VertexBufferSize[4] = sizeof(BBVert_Line);
	VertexBufferSize[5] = sizeof(travelVert);
	VertexBufferSize[6] = sizeof(tangent);
	VertexBufferSize[7] = sizeof(normal);
	VertexBufferSize[8] = sizeof(binormal);
	IndexBufferSize[obj] = sizeof(Indices);
	NumIdcs[obj] = IndexCount;
	
	createVAOs(Vertices, Indices, obj);
	createVAOs(WhiteVertices, Indices, 1);
	createVAOs(SubVertices[0], SubIndices, 2);
	createVAOs(MidVertices, Indices, 3);
	createVAOs(BBVert_Line, BBIndices, 4);
	createVAOs(travelVert, Indices, 5);
	createVAOs(tangent, Indices, 6);
	createVAOs(normal, Indices, 7);
	createVAOs(binormal, Indices, 8);
}

// this actually creates the VAO (structure) and the VBO (vertex data buffer)
void createVAOs(Vertex Vertices[], GLushort Indices[], int ObjectId) {
	GLenum ErrorCheckValue = glGetError();
	const size_t VertexSize = sizeof(Vertices[0]);
	const size_t RgbOffset = sizeof(Vertices[0].Position);

	// Create Vertex Array Object
	glGenVertexArrays(1, &VertexArrayId[ObjectId]);
	glBindVertexArray(VertexArrayId[ObjectId]);

	// Create buffer for vertex data
	glGenBuffers(1, &VertexBufferId[ObjectId]);
	glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[ObjectId]);
	glBufferData(GL_ARRAY_BUFFER, VertexBufferSize[ObjectId], Vertices, GL_STATIC_DRAW);

	// Create buffer for indices
	if (Indices != NULL) {
		glGenBuffers(1, &IndexBufferId[ObjectId]);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IndexBufferId[ObjectId]);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, IndexBufferSize[ObjectId], Indices, GL_STATIC_DRAW);
	}

	// Assign vertex attributes
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, VertexSize, 0);
	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, VertexSize, (GLvoid*)RgbOffset);

	glEnableVertexAttribArray(0);	// position
	glEnableVertexAttribArray(1);	// color

	// Disable our Vertex Buffer Object 
	glBindVertexArray(0);

	ErrorCheckValue = glGetError();
	if (ErrorCheckValue != GL_NO_ERROR)
	{
		fprintf(
			stderr,
			"ERROR: Could not create a VBO: %s \n",
			gluErrorString(ErrorCheckValue)
		);
	}
}

void createObjects(void) {
	// ATTN: DERIVE YOUR NEW OBJECTS HERE:  each object has
	// an array of vertices {pos;color} and
	// an array of indices (no picking needed here) (no need for indices)
	// ATTN: Project 1A, Task 1 == Add the points in your scene
	
	if (startup)
	{

		for (int i = 0; i < IndexCount; i++)
		{
			int scale = 1;
			Vertices[i] = { {scale*cos(2*PI*i/IndexCount),scale*sin(2*PI*i/IndexCount),0.0f, 1.0f},{0.0f, 1.0f, 0.0f, 1.0f} };
			WhiteVertices[i] = { {scale*cos(2*PI*i/IndexCount),scale*sin(2*PI*i/IndexCount),0.0f, 1.0f},{1.0f, 1.0f, 1.0f, 1.0f} };
		}

		for (int i = 0; i < IndexCount; i++)
		{
			Indices[i] = i;
		}

		for (int i = 0; i < MAXSIZE; i++)
		{
			SubIndices[i] = i;
		}

		float yellowcolor[4] = { 1.0f, 1.0f, 0.0f, 1.0f };
		float redcolor[4] = { 1.0f, 0.0f, 0.0f, 1.0f };
		float greencolor[4] = { 0.0f, 1.0f, 0.0f, 1.0f };
		float bluecolor[4] = { 0.0f, 0.0f, 1.0f, 1.0f };

		travelVert[0].SetCoords(BBVertices[0].Position);
		travelVert[0].SetColor(yellowcolor);

		for (int i = 0; i < 2; i++)
		{
			tangent[i].SetColor(redcolor);
			normal[i].SetColor(greencolor);
			binormal[i].SetColor(bluecolor);
		}

		startup = false;
	}

	// ATTN: Project 1B, Task 1 == create line segments to connect the control points

	point p1, p2, result;
	float subdivColor[4] = {0.0f, 0.0f, 1.0f, 1.0f};
	int temp_sz;
	sz = (int)IndexCount * pow(2, subdiv_k);

	for (int i = 0; i < IndexCount; i++) //First level of subdivision, k = 0
	{
		SubVertices[0][i].SetCoords(Vertices[i].Position);
		SubVertices[0][i].SetColor(subdivColor);
	}

	for (int j = 1; j <= subdiv_k; j++) //further subdivision, k = 1 to k = 5
	{
		for (int i = 0; i < IndexCount * pow(2, j - 1); i++)
		{
			p1 = point(SubVertices[j - 1][i].Position);
			if (i == IndexCount * pow(2, j - 1) - 1) p2 = point(SubVertices[j - 1][0].Position);
			else p2 = point(SubVertices[j - 1][i + 1].Position);


			result = ((p1 * 3) + p2) / 4;
			SubVertices[j][2 * i].SetCoords(result.toArray());
			SubVertices[j][2 * i].SetColor(subdivColor);

			result = (p1 + (p2 * 3)) / 4;
			temp_sz = (int)IndexCount * pow(2, j);
			temp_sz = (2*i + 1) % temp_sz;
			SubVertices[j][temp_sz].SetCoords(result.toArray());
			SubVertices[j][temp_sz].SetColor(subdivColor);
		}

	}

	// ATTN: Project 1B, Task 2 == create the vertices associated to the smoother curve generated by subdivision

	// ATTN: Project 1B, Task 4 == create the BB control points and apply De Casteljau's for their corresponding for each piece

	float bbcolor[4] = { 1.0f, 1.0f, 0.0f, 1.0f };
	float midcolor[4] = {1.0f, 0.0f, 0.0f, 1.0f};

	for (int i = 0; i < IndexCount; i++) //Initialize BB control points and midpoints
	{
		BBVertices[2*i].SetCoords(Vertices[i].Position);
		p1 = point(Vertices[i].Position);
		if (i == IndexCount - 1) p2 = point(Vertices[0].Position);
		else p2 = point(Vertices[i + 1].Position);

		result = (p1 + p2) / 2;
		BBVertices[2 * i + 1].SetCoords(result.toArray());
		MidVertices[i].SetCoords(result.toArray());

		BBVertices[2 * i].SetColor(bbcolor);
		BBVertices[2 * i + 1].SetColor(bbcolor);
		MidVertices[i].SetColor(midcolor);
	}

	point BBcoefficients[IndexCount][3];

	for (int i = 0; i < IndexCount; i++) //Create BB curve coefficients
	{
		/*
		{c(i,0), c(i,1), c(i,2)} represents a BB coefficient triplet.
		c(i,1) = Vertices[i];   //Vertices[i + 1] is Vertices[0] at end
				                //further implementation requires conversion to point struct first
		c(i,0) = (Vertices[i - 1] + Vertices[i]) / 2;
		c(i,2) = (Vertices[i + 1] + Vertices[i + 2]) / 2;
		// we then re-initialize c(i,0) and c(i,2)
		c(i,0) = (c(i,0) + c(i,1)) / 2;
		c(i,2) = (c(i,1) + c(i,2)) / 2;
		// Because c(i,0) and c(i,2) are based on midpoints adjacent to c(i,1), we can thus see that c(i,2) = c(i + 1,0) for i < n - 1
		*/

		if (i == 0) BBcoefficients[i][0] = (point(Vertices[IndexCount - 1].Position) + point(Vertices[i].Position)) / 2;
		else BBcoefficients[i][0] = (point(Vertices[i - 1].Position) + point(Vertices[i].Position)) / 2;
		BBcoefficients[i][1] = point(Vertices[i].Position);
		BBcoefficients[i][2] = (point(Vertices[i].Position) + point(Vertices[(i + 1) % IndexCount].Position)) / 2;
	}

	point result2, result3;
	for (int i = 0; i < IndexCount; i++)
	{
		for (int j = 0; j <= BIG_N; j++)
		{
			float t = float(j) / (float)BIG_N;
			result = BBcoefficients[i][0] * pow((1-t),2);
			result2 = BBcoefficients[i][1] * 2 * t * (1 - t);
			result3 = BBcoefficients[i][2] * pow(t,2);

			result = result + result2 + result3;
			BBVert_Line[i * BIG_N + j].SetCoords(result.toArray());
			BBVert_Line[i * BIG_N + j].SetColor(bbcolor);
		}
	}

	// ATTN: Project 1C, Task 3 == set coordinates of yellow point based on BB curve and perform calculations to find
	// the tangent, normal, and binormal
}

void pickVertex(bool method) {
	// Clear the screen in white
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(pickingProgramID);
	{
		glm::mat4 ModelMatrix = glm::mat4(1.0); // initialization
		if (double_view)
		{
			ModelMatrix = glm::translate(glm::mat4(), glm::vec3(0.0f, 1.5f, 0.0f)) * ModelMatrix;
		}
		// ModelMatrix == TranslationMatrix * RotationMatrix;
		glm::mat4 MVP = gProjectionMatrix * gViewMatrix * ModelMatrix;
		// MVP should really be PVM...
		// Send the MVP to the shader (that is currently bound)
		// as data type uniform (shared by all shader instances)
		glUniformMatrix4fv(PickingMatrixID, 1, GL_FALSE, &MVP[0][0]);

		// pass in the picking color array to the shader
		glUniform1fv(pickingColorArrayID, IndexCount, pickingColor);

		glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[0]); //added code for new buffer- Kyle

		// --- enter vertices into VBO and draw
		glEnable(GL_PROGRAM_POINT_SIZE);
		glBindVertexArray(VertexArrayId[0]);
		glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[0], Vertices);	// update buffer data
		glDrawElements(GL_POINTS, NumIdcs[0], GL_UNSIGNED_SHORT, (void*)0);

		glBindVertexArray(0);
	}
	glUseProgram(0);
	glFlush();
	// --- Wait until all the pending drawing commands are really done.
	// Ultra-mega-over slow ! 
	// There are usually a long time between glDrawElements() and
	// all the fragments completely rasterized.
	glFinish();

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	// --- Read the pixel at the center of the screen.
	// You can also use glfwGetMousePos().
	// Ultra-mega-over slow too, even for 1 pixel, 
	// because the framebuffer is on the GPU.
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);
	unsigned char data[4];  // 2x2 pixel region
	glReadPixels(xpos, window_height - ypos, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, data);
       	// window_height - ypos;  
	// OpenGL renders with (0,0) on bottom, mouse reports with (0,0) on top

	// Convert the color back to an integer ID
	gPickedIndex = int(data[0]);

	

	if (gPickedIndex >= 0 && gPickedIndex < IndexCount)
	{
		if (method == 0)
		{
			storeColorFunct(gPickedIndex);
			storeIndex = gPickedIndex;

			float newColor[4];
			for (int i = 0; i < 3; i++) newColor[i] = 1.0f;
			newColor[3] = 1.0f;
			Vertices[gPickedIndex].SetColor(newColor);
		}
		else if (method == 1)
		{
			float greencolor[4] = { 0.0f, 1.0f, 0.0f, 1.0f };
			Vertices[storeIndex].SetColor(greencolor);
		}
	}

	createObjects();	// re-evaluate curves in case vertices have been moved
	renderScene();

	// ATTN: Project 1A, Task 2
	// Find a way to change color of selected vertex and
	// store original color

	// Uncomment these lines if you wan to see the picking shader in effect
	//glfwSwapBuffers(window);
	//continue; // skips the visible rendering
}

void storeColorFunct(GLuint gPickedIndex){
	for (int i = 0; i < 3; i++) storeColor[i] = Vertices[gPickedIndex].Color[i];
}

// ATTN: Project 1A, Task 3 == Retrieve your cursor position, get corresponding world coordinate, and move the point accordingly

// ATTN: Project 1C, Task 1 == Keep track of z coordinate for selected point and adjust its value accordingly based on if certain
// buttons are being pressed
void moveVertex(void) {
	glm::mat4 ModelMatrix = glm::mat4(1.0);
	if (double_view)
	{
		ModelMatrix = glm::translate(glm::mat4(), glm::vec3(0.0f, 1.5f, 0.0f)) * ModelMatrix;
	}
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	glm::vec4 vp = glm::vec4(viewport[0], viewport[1], viewport[2], viewport[3]);

	if (gPickedIndex >= IndexCount) { 
		// Any number > vertices-indices is background!
		gMessage = "background";
	}
	else {
		std::ostringstream oss;
		oss << "point " << gPickedIndex;
		gMessage = oss.str();

		double xpos = 0;
		double ypos = 0;
		glfwGetCursorPos(window, &xpos, &ypos);
		glm::vec3 mouseCoords = glm::unProject(glm::vec3(xpos, ypos, 0.0f), ModelMatrix, gProjectionMatrix, vp);
		mouseCoords[0] *= -1;
		mouseCoords[1] *= -1;

		float newCoords[4];

		if (shift_press == true)
		{
			newCoords[0] = Vertices[gPickedIndex].Position[0];
			newCoords[1] = Vertices[gPickedIndex].Position[1];
			newCoords[2] = mouseCoords[0] - Vertices[gPickedIndex].Position[0];
			newCoords[3] = 1.0f;
		}
		else
		{
			newCoords[0] = mouseCoords[0];
			newCoords[1] = mouseCoords[1];
			if (double_view)
			{
				newCoords[1] -= 3.0f;
			}
			newCoords[2] = Vertices[gPickedIndex].Position[2];
			newCoords[3] = 1.0f;
		}

		Vertices[gPickedIndex].SetCoords(newCoords);
		WhiteVertices[gPickedIndex].SetCoords(newCoords);
	}
}

void renderScene(void) {    
	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.4f, 0.0f);
	// Re-clear the screen for visible rendering
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(programID);
	{
		// see comments in pick
		glm::mat4 ModelMatrix = glm::mat4(1.0); 
		if (double_view)
		{
			ModelMatrix = glm::translate(glm::mat4(), glm::vec3(0.0f, 1.5f, 0.0f)) * ModelMatrix;
		}
		glm::mat4 MVP = gProjectionMatrix * gViewMatrix * ModelMatrix;
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
		glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		glUniformMatrix4fv(ViewMatrixID, 1, GL_FALSE, &gViewMatrix[0][0]);
		
		glEnable(GL_PROGRAM_POINT_SIZE);

		glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[0]); //added code for new buffer- Kyle

		glBindVertexArray(VertexArrayId[0]);	// Draw Vertices
		glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[0], Vertices);		// Update green buffer data
		glDrawElements(GL_POINTS, NumIdcs[0], GL_UNSIGNED_SHORT, (void*)0);
		// // If don't use indices
		// glDrawArrays(GL_POINTS, 0, NumVerts[0]);
		
		glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[1]);
		glBindVertexArray(VertexArrayId[1]);
		glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[1], WhiteVertices);		// Update white buffer data
		glDrawArrays(GL_LINE_LOOP, 0, 8);
		
		glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[2]);
		glBindVertexArray(VertexArrayId[2]);
		glBufferSubData(GL_ARRAY_BUFFER, 0, MAXSIZE * sizeof(SubVertices[0][0]), SubVertices[subdiv_k]);		// Update blue buffer data
		glDrawArrays(GL_LINE_LOOP, 0, sz);

		if (bezier_render == true)
		{
			glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[3]);
			glBindVertexArray(VertexArrayId[3]);	// Draw Vertices
			glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[3], MidVertices);		// Update red buffer data
			glDrawArrays(GL_POINTS, 0, 8);

			glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[4]);
			glBindVertexArray(VertexArrayId[4]);
			glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[4], BBVert_Line);		// Update yellow buffer data
			glDrawArrays(GL_LINE_LOOP, 0, (IndexCount*BIG_N));
		}

		if (animation == true)
		{
			glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[5]);
			glBindVertexArray(VertexArrayId[5]);	// Draw Vertices
			glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[5], travelVert);		// Update animation buffer data
			glDrawArrays(GL_POINTS, 0, 1);

			glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[6]);
			glBindVertexArray(VertexArrayId[6]);	// Draw Vertices
			glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[6], tangent);		// Update tangent buffer data
			glDrawArrays(GL_LINE_STRIP, 0, 2);

			glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[8]);
			glBindVertexArray(VertexArrayId[8]);	// Draw Vertices
			glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[8], binormal);		// Update binormal buffer data
			glDrawArrays(GL_LINE_STRIP, 0, 2);

			glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[7]);
			glBindVertexArray(VertexArrayId[7]);	// Draw Vertices
			glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[7], normal);		// Update normal buffer data
			glDrawArrays(GL_LINE_STRIP, 0, 2);
		}
		
		glBindVertexArray(0);

		// ATTN: OTHER BINDING AND DRAWING COMMANDS GO HERE
		// one set per object:
		// glBindVertexArray(VertexArrayId[<x>]); etc etc

		// ATTN: Project 1C, Task 2 == Refer to https://learnopengl.com/Getting-started/Transformations and
		// https://learnopengl.com/Getting-started/Coordinate-Systems - draw all the objects associated with the
		// curve twice in the displayed fashion using the appropriate transformations
		if (double_view == true)
		{
			ModelMatrix = glm::mat4(1.0);
			glm::mat4 TranslationMatrix = glm::translate(glm::mat4(), glm::vec3(0.0f, -1.5f, 0.0f));
			glm::vec3 myRotationAxis(0.0f, 1.0f, 0.0f);
			ModelMatrix = glm::rotate(ModelMatrix, PI / 2, myRotationAxis);;
			ModelMatrix = TranslationMatrix * ModelMatrix;
			
			MVP = gProjectionMatrix * gViewMatrix * ModelMatrix;
			glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
			glUniformMatrix4fv(ViewMatrixID, 1, GL_FALSE, &gViewMatrix[0][0]);

			glEnable(GL_PROGRAM_POINT_SIZE);

			glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[0]); //added code for new buffer- Kyle

			glBindVertexArray(VertexArrayId[0]);	// Draw Vertices
			glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[0], Vertices);		// Update green buffer data
			glDrawElements(GL_POINTS, NumIdcs[0], GL_UNSIGNED_SHORT, (void*)0);
			// // If don't use indices
			// glDrawArrays(GL_POINTS, 0, NumVerts[0]);

			glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[1]);
			glBindVertexArray(VertexArrayId[1]);
			glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[1], WhiteVertices);		// Update white buffer data
			glDrawArrays(GL_LINE_LOOP, 0, 8);

			glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[2]);
			glBindVertexArray(VertexArrayId[2]);
			glBufferSubData(GL_ARRAY_BUFFER, 0, MAXSIZE * sizeof(SubVertices[0][0]), SubVertices[subdiv_k]);		// Update blue buffer data
			glDrawArrays(GL_LINE_LOOP, 0, sz);

			if (bezier_render == true)
			{
				glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[3]);
				glBindVertexArray(VertexArrayId[3]);	// Draw Vertices
				glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[3], MidVertices);		// Update red buffer data
				glDrawArrays(GL_POINTS, 0, 8);

				glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[4]);
				glBindVertexArray(VertexArrayId[4]);
				glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[4], BBVert_Line);		// Update yellow buffer data
				glDrawArrays(GL_LINE_LOOP, 0, (IndexCount * BIG_N));
			}

			if (animation == true)
			{
				glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[5]);
				glBindVertexArray(VertexArrayId[5]);	// Draw Vertices
				glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[5], travelVert);		// Update animation buffer data
				glDrawArrays(GL_POINTS, 0, 1);

				glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[6]);
				glBindVertexArray(VertexArrayId[6]);	// Draw Vertices
				glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[6], tangent);		// Update tangent buffer data
				glDrawArrays(GL_LINE_STRIP, 0, 2);

				glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[8]);
				glBindVertexArray(VertexArrayId[8]);	// Draw Vertices
				glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[8], binormal);		// Update binormal buffer data
				glDrawArrays(GL_LINE_STRIP, 0, 2);

				glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[7]);
				glBindVertexArray(VertexArrayId[7]);	// Draw Vertices
				glBufferSubData(GL_ARRAY_BUFFER, 0, VertexBufferSize[7], normal);		// Update normal buffer data
				glDrawArrays(GL_LINE_STRIP, 0, 2);
			}

			glBindVertexArray(0);
		}
	}
	glUseProgram(0);
	// Draw GUI
	TwDraw();

	// Swap buffers
	glfwSwapBuffers(window);
	glfwPollEvents();
}

void cleanup(void) {
	// Cleanup VBO and shader
	for (int i = 0; i < NumObjects; i++) {
		glDeleteBuffers(1, &VertexBufferId[i]);
		glDeleteBuffers(1, &IndexBufferId[i]);
		glDeleteVertexArrays(1, &VertexArrayId[i]);
	}
	glDeleteProgram(programID);
	glDeleteProgram(pickingProgramID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();
}

// Alternative way of triggering functions on mouse click and keyboard events
static void mouseCallback(GLFWwindow* window, int button, int action, int mods) {	
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		pickVertex(0);
	}
	else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
		pickVertex(1);
	}
}

static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (key == GLFW_KEY_1 && action == GLFW_PRESS) {
		key_1();
	}
	if (key == GLFW_KEY_2 && action == GLFW_PRESS) {
		key_2();
	}
	if (key == GLFW_KEY_4 && action == GLFW_PRESS) {
		key_4();
	}
	if (key == GLFW_KEY_5 && action == GLFW_PRESS) {
		key_5();
	}
	if (key == GLFW_KEY_R && action == GLFW_PRESS) {
		startup = true;
		animation = false;
		double_view = false;
		bezier_render = true;
		createObjects();
		renderScene();
	}
	if ((key == GLFW_KEY_RIGHT_SHIFT || key == GLFW_KEY_LEFT_SHIFT) && action == GLFW_PRESS)
	{
		shift_press = true;
	}
	if ((key == GLFW_KEY_RIGHT_SHIFT || key == GLFW_KEY_LEFT_SHIFT) && action == GLFW_RELEASE)
	{
		shift_press = false;
		float greencolor[4] = { 0.0f, 1.0f, 0.0f, 1.0f };
		for (int i = 0; i < IndexCount; i++)
		{
			Vertices[i].SetColor(greencolor);
		}
	}
}

void key_1() {
	if (subdiv_k >= maxsub) subdiv_k = 0;
	else subdiv_k += 1;

	createObjects();
	renderScene();
}

void key_2() {
	if (bezier_render == true)
	{
		bezier_render = false;
	}
	else
	{
		bezier_render = true;
	}

	createObjects();
	renderScene();
}

void key_4() {
	if (double_view == true)
	{
		double_view = false;
	}
	else
	{
		double_view = true;
	}

	createObjects();
	renderScene();
}

void key_5() {
	if (animation == false)
	{
		point start, end;
		point cross1, cross2;
		float magnitude;
		animation = true;
		timer = 18.0f;
		float timer_index;

		while (animation == true)
		{
			if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
			{
				break;
			}

			travelVert[0].SetCoords(BBVertices[(int)timer].Position);
			tangent[0].SetCoords(BBVertices[(int)timer].Position);
			normal[0].SetCoords(BBVertices[(int)timer].Position);
			binormal[0].SetCoords(BBVertices[(int)timer].Position);

			start = point(BBVertices[(int)timer].Position); //set tangent position
			if ((int)timer >= (IndexCount * (BIG_N + 1) - 1))
			{
				end = point(BBVertices[0].Position) - start;
			}
			else
			{
				end = point(BBVertices[((int)timer + 1)].Position) - start;
			}

			magnitude = pow(((end.x * end.x) + (end.y * end.y) + (end.z * end.z)), 0.5);
			end = (end / magnitude);

			end = end + start;
			tangent[1].SetCoords(end.toArray());

			//setting binormal as cross product of control lines
			
			timer_index = (timer + 1) / BIG_N;
			cross1 = point(MidVertices[(int)timer_index].Position);
			if ((int)timer_index >= IndexCount) cross2 = point(MidVertices[0].Position);
			else cross2 = point(MidVertices[((int)timer_index + 1)].Position);

			cross1 = cross1 - start;
			cross2 = cross2 - start;
			end.x = cross1.y * cross2.z - cross1.z * cross2.y;
			end.y = cross1.z * cross2.x - cross1.x * cross2.z;
			end.z = cross1.x * cross2.y - cross1.y * cross2.x;

			magnitude = pow(((end.x * end.x) + (end.y * end.y) + (end.z * end.z)), 0.5);
			end = (end / magnitude);

			end = end * -1;
			end = end + start;
			binormal[1].SetCoords(end.toArray());
			
			//setting normal as cross product of tangent and binormal
			cross1 = point(tangent[1].Position);
			cross2 = point(binormal[1].Position);

			cross1 = cross1 - start;
			cross2 = cross2 - start;
			end.x = cross1.y * cross2.z - cross1.z * cross2.y;
			end.y = cross1.z * cross2.x - cross1.x * cross2.z;
			end.z = cross1.x * cross2.y - cross1.y * cross2.x;

			magnitude = pow(((end.x * end.x) + (end.y * end.y) + (end.z * end.z)), 0.5);
			end = (end / magnitude);

			end = end + start;
			normal[1].SetCoords(end.toArray());

			renderScene();

			timer += 0.5f;
			if (timer >= (IndexCount * (BIG_N+1) - 1)) timer = 18.0f;
		}
	}
	else
	{
		animation = false;
		timer = 0.0f;
		createObjects();
		renderScene();
	}
}

int main(void) {
	// ATTN: REFER TO https://learnopengl.com/Getting-started/Creating-a-window
	// AND https://learnopengl.com/Getting-started/Hello-Window to familiarize yourself with the initialization of a window in OpenGL

	// Initialize window
	int errorCode = initWindow();
	if (errorCode != 0)
		return errorCode;

	// ATTN: REFER TO https://learnopengl.com/Getting-started/Hello-Triangle to familiarize yourself with the graphics pipeline
	// from setting up your vertex data in vertex shaders to rendering the data on screen (everything that follows)

	// Initialize OpenGL pipeline
	initOpenGL();

	double lastTime = glfwGetTime();
	int nbFrames = 0;
	do {
		// Timing 
		double currentTime = glfwGetTime();
		nbFrames++;
		if (currentTime - lastTime >= 1.0){ // If last prinf() was more than 1sec ago
			printf("%f ms/frame\n", 1000.0 / double(nbFrames));
			nbFrames = 0;
			lastTime += 1.0;
		}

		// DRAGGING: move current (picked) vertex with cursor
		if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT)) {
			moveVertex();
		}

		// ATTN: Project 1B, Task 2 and 4 == account for key presses to activate subdivision and hiding/showing functionality
		// for respective tasks

		// DRAWING the SCENE
		createObjects();	// re-evaluate curves in case vertices have been moved
		renderScene();

	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
	glfwWindowShouldClose(window) == 0);

	cleanup();

	return 0;
}
