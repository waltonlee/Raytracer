#define NUM_OPENGL_LIGHTS 8

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <GL/glui.h>
#include "Shape.h"
#include "Cube.h"
#include "Cylinder.h"
#include "Cone.h"
#include "Sphere.h"
#include "SceneParser.h"
#include "Camera.h"
#include "FlatPrim.h"
#include "raytracer.h"

using namespace std;

/** These are the live variables passed into GLUI ***/
int  isectOnly = 0;

int	 camRotU = 0;
int	 camRotV = 0;
int	 camRotW = 0;
int  viewAngle = 45;
float eyeX = 2;
float eyeY = 2;
float eyeZ = 2;
float lookX = -2;
float lookY = -2;
float lookZ = -2;
int recursion = 0;

std::vector<string> filenames;
std::vector<Image> images;

/** These are GLUI control panel objects ***/
int  main_window;
string filenamePath = "data/tests/all.xml";
GLUI_EditText* filenameTextField = NULL;
GLubyte* pixels = NULL;
int pixelWidth = 0, pixelHeight = 0;
int screenWidth = 0, screenHeight = 0;

/** these are the global variables used for rendering **/
Cube* cube = new Cube();
Cylinder* cylinder = new Cylinder();
Cone* cone = new Cone();
Sphere* sphere = new Sphere();
SceneParser* parser = NULL;
Camera* camera = new Camera();
PrimList* primlist;

void setupCamera();
void updateCamera();

void storeImage(string name) {
	string word;
	float max;
	int iter = 0;
	filenames.push_back(name);
	Image i;
	ifstream file;
	file.open(name);
	if (file.is_open()) {
		//ignore header
		getline(file, word);
		getline(file, word);
		file >> i.w;
		file >> i.h;
		file >> word;
		max = stof(word);
		i.pixels = new float[i.w * i.h * 3];
		while (file >> word) {
			i.pixels[iter] = stof(word) / max;
			iter++;
	    }
		images.push_back(i);
	    file.close();
	}
}

void setupImages() {
	bool found;
	string filename;
    for(std::vector<FlatPrim*>::iterator it = primlist->prims.begin(); it != primlist->prims.end(); ++it) {
        for(std::vector<ScenePrimitive*>::iterator scenePrim = (*it)->prims.begin(); scenePrim != (*it)->prims.end(); ++scenePrim) {
			if ((*scenePrim)->material.textureMap->isUsed) {
				found = false;
				filename = (*scenePrim)->material.textureMap->filename;
		        for(std::vector<string>::iterator name = filenames.begin(); name != filenames.end(); ++name) {
					if ((*name) == filename) {
						found = true;
						break;
					}
				}
				if (!found) {
					storeImage(filename);
				}
			}
		}
	}
	cout << "done storing images" << endl;
}

void setPixel(GLubyte* buf, int x, int y, int r, int g, int b) {
	buf[(y*pixelWidth + x) * 3 + 0] = (GLubyte)r;
	buf[(y*pixelWidth + x) * 3 + 1] = (GLubyte)g;
	buf[(y*pixelWidth + x) * 3 + 2] = (GLubyte)b;
}

float magnitude(Vector v) {
	return sqrt(SQR(v[0]) + SQR(v[1]) + SQR(v[2]));
}

Shape *getShape(int shapeType) {
	switch (shapeType) {
	case SHAPE_CUBE:
		return cube;
	case SHAPE_CYLINDER:
		return cylinder;
	case SHAPE_CONE:
		return cone;
	case SHAPE_SPHERE:
		return sphere;
	default:
		return cube;
	}
}

/* assumes min_t is -1.0 to begin with */
void findNearest(double &min_t, FlatPrim *&nearest, ScenePrimitive *&near_prim, Point eyeP, Vector rayV) {
    for(std::vector<FlatPrim*>::iterator it = primlist->prims.begin(); it != primlist->prims.end(); ++it) {
		Point objEye = (*it)->inverted * eyeP;
		Vector objRay = (*it)->inverted * rayV;
        for(std::vector<ScenePrimitive*>::iterator scenePrim = (*it)->prims.begin(); scenePrim != (*it)->prims.end(); ++scenePrim) {
			Shape *shape = getShape((*scenePrim)->type);
			double t = shape->Intersect(objEye, objRay, Matrix());
			if ((t < min_t || IN_RANGE(min_t, -1.0)) && t >= 0) {
				min_t = t;
				nearest = (*it);
				near_prim = (*scenePrim);
			}
		}
	}
}

// currently stores texture mapped values in Ir, Ig, and Ib when passed in
// otherwise they are 0
void getLightIlluminosity(int lightnum, Vector N, Point intersection, ScenePrimitive *prim, SceneGlobalData globalData, float &Ir, float &Ig, float &Ib, bool is_texture, float blend) {
	float rDiffuse = Ir; float rSpec = Ir;
	float gDiffuse = Ig; float gSpec = Ig;
	float bDiffuse = Ib; float bSpec = Ib;
	SceneLightData light;
	Vector Lm;
	Vector R;
	Vector V = camera->GetEyePoint() - intersection;
	V.normalize();
	double t = -1;
	double light_t;
	FlatPrim *nearest = NULL;
	ScenePrimitive *near_prim = NULL;
	parser->getLightData(lightnum, light);
	if (light.type == LIGHT_POINT) {
		Lm = light.pos - intersection;
		Lm.normalize();
		findNearest(t, nearest, near_prim, intersection + EPSILON * Lm, Lm);
		light_t = magnitude(light.pos-intersection);
		if (light_t > t && !IN_RANGE(t, -1)) {
			Ir = 0; Ig = 0; Ib = 0;
			return;
		}
	}
	else if (light.type == LIGHT_DIRECTIONAL) {
		Lm = -1 * light.dir;
		Lm.normalize();
		findNearest(t, nearest, near_prim, intersection + EPSILON * Lm, Lm);
		if (nearest != NULL) {
			Ir = 0; Ig = 0; Ib = 0;
			return;
		}
	}
	else {
		std::cerr << "not a point or directional light\n";
		//only handling point lights for this assignment
		return;
	}
	float N_dot_Lm = dot(N, Lm);
	Vector temp = -1 * Lm;
	R = temp - 2 * dot(temp, N) * N;
	float R_dot_V = dot(R, V);
	if (prim->material.shininess > 1) {
		R_dot_V = pow(R_dot_V, prim->material.shininess);
	}
	if (!is_texture) {
		rDiffuse = prim->material.cDiffuse.r;
		gDiffuse = prim->material.cDiffuse.g;
		bDiffuse = prim->material.cDiffuse.b;
		rSpec = prim->material.cSpecular.r;
		gSpec = prim->material.cSpecular.g;
		bSpec = prim->material.cSpecular.b;
	}
	else {
		rDiffuse = prim->material.cDiffuse.r * (1.0-blend) + Ir * blend;
		gDiffuse = prim->material.cDiffuse.g * (1.0-blend) + Ig * blend;
		bDiffuse = prim->material.cDiffuse.b * (1.0-blend) + Ib * blend;
		rSpec = prim->material.cSpecular.r * (1.0-blend) + Ir * blend;
		bSpec = prim->material.cSpecular.b * (1.0-blend) + Ib * blend;
		gSpec = prim->material.cSpecular.g * (1.0-blend) + Ig * blend;
	}
	Ir = light.color.r * (globalData.kd * rDiffuse * N_dot_Lm
			+ globalData.ks * rSpec * R_dot_V);
	Ig = light.color.g * (globalData.kd * gDiffuse * N_dot_Lm
			+ globalData.ks * gSpec * R_dot_V);
	Ib = light.color.b * (globalData.kd * bDiffuse * N_dot_Lm
			+ globalData.ks * bSpec * R_dot_V);
}

void calculateIlluminosity(float &Ir, float &Ig, float &Ib, Point eyeP, Vector rayV, double t, SceneGlobalData &globalData, FlatPrim *nearest, ScenePrimitive *near_prim, int layer, float coeff) {
	//change to have the ui control recursion
	if (layer > recursion) {
		Ir = 0; Ig = 0; Ib = 0;
		return;
	}
	float r = -1; float g = -1; float b = -1;
	float text_r= -1; float text_g= -1; float text_b = -1;
	float red = -1; float green = -1; float blue = -1;
	float blend = 0;
	Shape *shape = getShape(near_prim->type);
	Point intersection = eyeP + rayV * t;
	Point objEye = nearest->inverted * eyeP;
	Vector objRay = nearest->inverted * rayV;
	bool is_texture = near_prim->material.textureMap->isUsed;
	if (is_texture) {
		int pointer = 0;
		string filename = near_prim->material.textureMap->filename;
        for(std::vector<string>::iterator name = filenames.begin(); name != filenames.end(); ++name) {
			if (*name == filename) {
				shape->getColor(objEye, objRay, t, text_r, text_g, text_b, images[pointer], near_prim->material.textureMap->repeatU, near_prim->material.textureMap->repeatV);
				break;
			}
			pointer++;
		}
		blend = near_prim->material.blend;
	}
	if (!is_texture) {
		r = near_prim->material.cAmbient.r;
		g = near_prim->material.cAmbient.g;
		b = near_prim->material.cAmbient.b;
	}
	else {
		r = near_prim->material.cAmbient.r * (1.0-blend) + text_r * blend;
		g = near_prim->material.cAmbient.g * (1.0-blend) + text_g * blend;
		b = near_prim->material.cAmbient.b * (1.0-blend) + text_b * blend;
	}
	Ir = globalData.ka * r;
	Ig = globalData.ka * g;
	Ib = globalData.ka * b;
	Matrix m = transpose(nearest->inverted);
	Vector N = shape->findIsectNormal(objEye, objRay, t);
	N = m * N;
	N.normalize();
	SceneLightData light;
	Vector Lm;
	for (int iter = 0; iter < parser->getNumLights(); iter++) {
		red = text_r; green = text_g; blue = text_b;
		getLightIlluminosity(iter, N, intersection, near_prim, globalData, red, green, blue, is_texture, blend);
		Ir += red;
		Ig += green;
		Ib += blue;
	}
	// calculate reflected illuminosity
	if (!is_texture && coeff > pow(10, -2)) {
		red = 0; green = 0; blue = 0;
		Vector reflected_ray = rayV - 2 * dot(rayV, N) * N;
		double reflected_t = -1;
		FlatPrim *reflected_flatprim = NULL;
		ScenePrimitive *reflected_prim = NULL;
		Point isectepsilon = intersection + EPSILON * reflected_ray;
		findNearest(reflected_t, reflected_flatprim, reflected_prim, isectepsilon, reflected_ray);
		if (reflected_t >= 0) {
			float max = fmax(near_prim->material.cReflective.r, near_prim->material.cReflective.g);
			max = fmax(max, near_prim->material.cReflective.b);
			coeff *= max * globalData.ks;
			calculateIlluminosity(red, green, blue, isectepsilon, reflected_ray, reflected_t, globalData, reflected_flatprim, reflected_prim, layer + 1, coeff);
			Ir += globalData.ks * near_prim->material.cReflective.r * red;
			Ig += globalData.ks * near_prim->material.cReflective.g * green;
			Ib += globalData.ks * near_prim->material.cReflective.b * blue;
		}
	}
	if (Ir > 1) Ir = 1;
	if (Ib > 1) Ib = 1;
	if (Ig > 1) Ig = 1;
	if (Ir < 0) Ir = 0;
	if (Ib < 0) Ib = 0;
	if (Ig < 0) Ig = 0;
}

void callback_start(int id) {
	cout << "start button clicked!" << endl;

	if (parser == NULL) {
		cout << "no scene loaded yet" << endl;
		return;
	}
	if (primlist == NULL) {
		cout << "primlist not loaded yet" << endl;
		return;
	}

	pixelWidth = screenWidth;
	pixelHeight = screenHeight;

	updateCamera();

	if (pixels != NULL) {
		delete pixels;
	}
	pixels = new GLubyte[pixelWidth  * pixelHeight * 3];
	memset(pixels, 0, pixelWidth * pixelHeight * 3);

	cout << "(w, h): " << pixelWidth << ", " << pixelHeight << endl;

	Point eyeP = camera->GetEyePoint();
	Vector rayV;
	double min_t;
	FlatPrim *nearest = NULL;
	ScenePrimitive *near_prim = NULL;
	SceneGlobalData globalData;
	parser->getGlobalData(globalData);

	for (int i = 0; i < pixelWidth; i++) {
		for (int j = 0; j < pixelHeight; j++) {
			nearest = NULL; near_prim = NULL;
			min_t = -1;
			rayV = generateRay(i, j, pixelWidth, pixelHeight, eyeP, camera);
			findNearest(min_t, nearest, near_prim, eyeP, rayV);
			if (IN_RANGE(min_t, -1.0)) {
				setPixel(pixels, i, j, 0, 0, 0);
			}
			else if (isectOnly) {
				setPixel(pixels, i, j, 255, 255, 255);
			}
			else {
				float Ir = 0; float Ig = 0; float Ib = 0;
				calculateIlluminosity(Ir, Ig, Ib, eyeP, rayV, min_t, globalData, nearest, near_prim, 0, 1.0);
				Ir *= 255;
				Ig *= 255;
				Ib *= 255;
				setPixel(pixels, i, j, Ir, Ig, Ib);
			}
		}
	}
	glutPostRedisplay();
}

void callback_load(int id) {
	char curDirName [2048];
	if (filenameTextField == NULL) {
		return;
	}
	printf ("%s\n", filenameTextField->get_text());

	if (parser != NULL) {
		delete parser;
	}
	if (primlist != NULL) {
		delete primlist;
	}
	parser = new SceneParser (filenamePath);
	cout << "success? " << parser->parse() << endl;
	primlist = new PrimList(parser->getRootNode());
	setupImages();

	setupCamera();
}

/***************************************** myGlutIdle() ***********/

void myGlutIdle(void)
{
	/* According to the GLUT specification, the current window is
	undefined during an idle callback.  So we need to explicitly change
	it if necessary */
	if (glutGetWindow() != main_window)
		glutSetWindow(main_window);

	glutPostRedisplay();
}


/**************************************** myGlutReshape() *************/

void myGlutReshape(int x, int y)
{
	float xy_aspect;

	xy_aspect = (float)x / (float)y;
	glViewport(0, 0, x, y);
	camera->SetScreenSize(x, y);

	screenWidth = x;
	screenHeight = y;

	glutPostRedisplay();
}


/***************************************** setupCamera() *****************/
void setupCamera()
{
	SceneCameraData cameraData;
	parser->getCameraData(cameraData);

	camera->Reset();
	camera->SetViewAngle(cameraData.heightAngle);
	if (cameraData.isDir == true) {
		camera->Orient(cameraData.pos, cameraData.look, cameraData.up);
	}
	else {
		camera->Orient(cameraData.pos, cameraData.lookAt, cameraData.up);
	}

	viewAngle = camera->GetViewAngle();
	Point eyeP = camera->GetEyePoint();
	Vector lookV = camera->GetLookVector();
	eyeX = eyeP[0];
	eyeY = eyeP[1];
	eyeZ = eyeP[2];
	lookX = lookV[0];
	lookY = lookV[1];
	lookZ = lookV[2];
	camRotU = 0;
	camRotV = 0;
	camRotW = 0;
	GLUI_Master.sync_live_all();
}

void updateCamera()
{
	camera->Reset();

	Point guiEye (eyeX, eyeY, eyeZ);
	Point guiLook(lookX, lookY, lookZ);
	Vector guiUp = camera->GetUpVector();
	camera->SetViewAngle(viewAngle);
	camera->Orient(guiEye, guiLook, guiUp);
	camera->RotateU(camRotU);
	camera->RotateV(camRotV);
	camera->RotateW(camRotW);
}

/***************************************** myGlutDisplay() *****************/

void myGlutDisplay(void)
{
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	if (parser == NULL) {
		return;
	}

	if (pixels == NULL) {
		return;
	}

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glDrawPixels(pixelWidth, pixelHeight, GL_RGB, GL_UNSIGNED_BYTE, pixels);
	glutSwapBuffers();
}

void onExit()
{
	delete cube;
	delete cylinder;
	delete cone;
	delete sphere;
	delete camera;
	if (parser != NULL) {
		delete parser;
	}
	if (pixels != NULL) {
		delete pixels;
	}
	if (primlist != NULL) {
		delete primlist;
	}
}

/**************************************** main() ********************/

int main(int argc, char* argv[])
{
	atexit(onExit);

	/****************************************/
	/*   Initialize GLUT and create window  */
	/****************************************/

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowPosition(50, 50);
	glutInitWindowSize(500, 500);

	main_window = glutCreateWindow("COMP 175 Assignment 4");
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutReshape);

	/****************************************/
	/*         Here's the GLUI code         */
	/****************************************/

	GLUI* glui = GLUI_Master.create_glui("GLUI");

	filenameTextField = new GLUI_EditText( glui, "Filename:", filenamePath);
	filenameTextField->set_w(300);
	glui->add_button("Load", 0, callback_load);
	glui->add_button("Start!", 0, callback_start);
	glui->add_checkbox("Isect Only", &isectOnly);

	GLUI_Panel *recursion_panel = glui->add_panel("Recursion");
	(new GLUI_Spinner(recursion_panel, "Recursion:", &recursion))
		->set_int_limits(0, 4);
	glui->add_column_to_panel(recursion_panel, true);



	GLUI_Panel *camera_panel = glui->add_panel("Camera");
	(new GLUI_Spinner(camera_panel, "RotateV:", &camRotV))
		->set_int_limits(-179, 179);
	(new GLUI_Spinner(camera_panel, "RotateU:", &camRotU))
		->set_int_limits(-179, 179);
	(new GLUI_Spinner(camera_panel, "RotateW:", &camRotW))
		->set_int_limits(-179, 179);
	(new GLUI_Spinner(camera_panel, "Angle:", &viewAngle))
		->set_int_limits(1, 179);

	glui->add_column_to_panel(camera_panel, true);

	GLUI_Spinner* eyex_widget = glui->add_spinner_to_panel(camera_panel, "EyeX:", GLUI_SPINNER_FLOAT, &eyeX);
	eyex_widget->set_float_limits(-10, 10);
	GLUI_Spinner* eyey_widget = glui->add_spinner_to_panel(camera_panel, "EyeY:", GLUI_SPINNER_FLOAT, &eyeY);
	eyey_widget->set_float_limits(-10, 10);
	GLUI_Spinner* eyez_widget = glui->add_spinner_to_panel(camera_panel, "EyeZ:", GLUI_SPINNER_FLOAT, &eyeZ);
	eyez_widget->set_float_limits(-10, 10);

	GLUI_Spinner* lookx_widget = glui->add_spinner_to_panel(camera_panel, "LookX:", GLUI_SPINNER_FLOAT, &lookX);
	lookx_widget->set_float_limits(-10, 10);
	GLUI_Spinner* looky_widget = glui->add_spinner_to_panel(camera_panel, "LookY:", GLUI_SPINNER_FLOAT, &lookY);
	looky_widget->set_float_limits(-10, 10);
	GLUI_Spinner* lookz_widget = glui->add_spinner_to_panel(camera_panel, "LookZ:", GLUI_SPINNER_FLOAT, &lookZ);
	lookz_widget->set_float_limits(-10, 10);

	glui->add_button("Quit", 0, (GLUI_Update_CB)exit);

	glui->set_main_gfx_window(main_window);

	/* We register the idle callback with GLUI, *not* with GLUT */
	GLUI_Master.set_glutIdleFunc(myGlutIdle);

	glutMainLoop();

	return EXIT_SUCCESS;
}
