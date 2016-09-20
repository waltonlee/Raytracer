#ifndef SHAPE_H
#define SHAPE_H

#include <GL/glui.h>
#include "Algebra.h"
#include <iostream>

// TODO move later to another file
struct Image {
	int w;
	int h;
	float *pixels;
};

class Shape {
public:
	Shape() {
		crossection = NULL;
		crossNormals = NULL;
		bodyPoints = NULL;
		bodyNormals = NULL;
		isRotate = true;
		calculated = false;
		fans = NULL;
		fanLine = NULL;
		numFans = 0;
		fanLineNormals = NULL;
		fanNormals = NULL;
	};
	~Shape() {
		clearBodyPoints();
		clearCrossection();
		clearFans();
	};

	virtual double Intersect(Point eyePointP, Vector rayV, Matrix transformMatrix) = 0;
	virtual Vector findIsectNormal(Point eyePoint, Vector ray, double dist) = 0;
	virtual void getColor(Point eye, Vector ray, double time, float &red, float &green, float &blue, Image &image, int i, int j) = 0;

	void setSegments(int x, int y) {
		if (isRotate) {
			if (m_segmentsY != y) {
				m_segmentsX = x;
				m_segmentsY = y;
				calculateCrossection();
				rotateCrossection();
				calculateFanLine();
				rotateFans();
				calculated = true;
			}
			else if (m_segmentsX != x) {
				m_segmentsX = x;
				rotateCrossection();
				rotateFans();
				calculated = true;
			}
		} else {
			if (m_segmentsX != x || m_segmentsY != y) {
				m_segmentsX = x;
				m_segmentsY = y;
				calculateSides();
				calculateTop();
				calculated = true;
			}
		}
	};

	virtual void draw() {
		if (isRotate) {
			if (!calculated) {
				calculateCrossection();
				rotateCrossection();
				calculateFanLine();
				rotateFans();
				calculated = true;
			}
			drawBody();
			drawFans();
			//exit(1);
		} else {
			if (!calculated) {
				calculateSides();
				calculateTop();
				calculated = true;
			}
			drawSides();
			drawTop();
		}
	};

	virtual void drawNormal() {
		if (isRotate) {
			if (!calculated) {
				calculateCrossection();
				rotateCrossection();
				calculateFanLine();
				rotateFans();
				calculated = true;
			}
			drawBodyNormals();
			drawFanNormals();
		}
		else {
			if (!calculated) {
				calculateSides();
				calculateTop();
				calculated = true;
			}
			drawSideNormals();
			drawTopNormals();
		}
	};

protected:

	double quadraticFormula(double &a, double &b, double &c) {
		double discriminate = SQR(b) - 4 * a * c;
		if (discriminate < 0) {
			return -1;
		}
		double d1 = (-b + sqrt(discriminate)) / (2 * a);
		double d2 = (-b - sqrt(discriminate)) / (2 * a);
		double result = -1;
		if (d1 >= 0 && d2 >= 0) {
			result = fmin(d1, d2);
		} else if (d1 >= 0) {
			result = d1;
		} else if (d2 >= 0) {
			result = d2;
		}
		return result;
	};

	double planeIntersection(Point eyePointP, Vector rayV, int plane, double planeVal) {
		double eyePointVal = eyePointP[plane];
		double rayVal = rayV[plane];
		if (IN_RANGE(rayVal, 0) || (planeVal - eyePointVal) * rayVal < 0) {
			return -1;
		}
		return (planeVal - eyePointVal) / rayVal;
	};

	void clearCrossection() {
		if (crossection != NULL) {
			delete [] crossection;
		}
		if (crossNormals != NULL) {
			delete [] crossNormals;
		}
	};

	void clearBodyPoints() {
		if (bodyPoints != NULL) {
			delete [] bodyPoints;
		}
		if (bodyNormals != NULL) {
			delete [] bodyNormals;
		}
	};

	void clearFanLine() {
		if (fanLine != NULL) {
			delete [] fanLine;
		}
		if (fanLineNormals != NULL) {
			delete [] fanLineNormals;
		}
	};

	void clearFans() {
		if (fans != NULL) {
			delete [] fans;
		}
		if (fanNormals != NULL) {
			delete [] fanNormals;
		}
	};

	void rotateCrossection() {
		clearBodyPoints();
		bodyPoints = new Point[m_segmentsX * crossLength()];
		bodyNormals = new Vector[m_segmentsX * crossLength()];
		float theta = 2 * PI / m_segmentsX;
		for (int i = 0; i < m_segmentsX; i++) {
			for (int j = 0; j < crossLength(); j++) {
				bodyPoints[i * crossLength() + j] = rotY_mat(theta * i) * crossection[j];
				bodyNormals[i * crossLength() + j] = rotY_mat(theta * i) * crossNormals[j];
			}
		}
	};

	void rotateFans() {
		clearFans();
		fans = new Point[numFans * (m_segmentsX + 1)];
		fanNormals = new Vector[numFans * (m_segmentsX + 1)];
		float theta = 2 * PI / m_segmentsX;
		for (int i = 0; i < numFans; i++) {
			fans[i * (m_segmentsX + 1)] = fanLine[i * 2];
			fanNormals[i * (m_segmentsX + 1)] = fanLineNormals[i * 2];
			Point p = fanLine[i * 2 + 1];
			Vector v = fanLineNormals[i * 2 + 1];
			for (int j = 0; j < m_segmentsX; j++) {
				fans[i * (m_segmentsX + 1) + j + 1] = rotY_mat(theta * j) * p;
				fanNormals[i * (m_segmentsX + 1) + j + 1] = rotY_mat(theta * j) * v;
			}
		}
	};

	/* size of crossection array */
	virtual int crossLength() {
		return 0;
	};

	virtual void calculateCrossection() {};
	virtual void calculateFanLine() {};
	virtual void calculateSides() {};
	virtual void calculateTop() {};

	/* crossection to rotate, points from top to bottom */
	Point *crossection;
	Point *bodyPoints;
	Vector *crossNormals;
	Vector *bodyNormals;
	Point *fanLine;
	Point *fans;
	Vector *fanLineNormals;
	Vector *fanNormals;
	int numFans;
	int numSides;
	bool isRotate;

	/* WARNING, we never reset this to false except for the constructor */
	bool calculated;

	int m_segmentsX, m_segmentsY;

	int width = 1;

private:
	void drawBody() {
		int numPoints = crossLength() * m_segmentsX;
		Vector v;
		for (int i = 0; i < crossLength() * m_segmentsX; i = i + crossLength()) {
			glBegin(GL_TRIANGLE_STRIP);
			v = bodyNormals[i];
			glNormal3dv(v.unpack());
			glVertex3d(bodyPoints[i][0], bodyPoints[i][1], bodyPoints[i][2]);
			v = bodyNormals[(i+crossLength())%numPoints];
			glNormal3dv(v.unpack());
			glVertex3d(bodyPoints[(i+crossLength())%numPoints][0], bodyPoints[(i+crossLength())%numPoints][1], bodyPoints[(i+crossLength())%numPoints][2]);
			for (int j = i+1; j < i + crossLength(); j++) {
				v = bodyNormals[j];
				glNormal3dv(v.unpack());
				Point p = bodyPoints[j];
				glVertex3dv(p.unpack());
				v = bodyNormals[(j+crossLength())%numPoints];
				glNormal3dv(v.unpack());
				p = bodyPoints[(j+crossLength())%numPoints];
				glVertex3dv(p.unpack());
			}
			glEnd();
		}
	};

	void drawBodyNormals() {
		for (int i = 0; i < m_segmentsX * crossLength(); i++) {
			glBegin(GL_LINES);
				Point p = bodyPoints[i];
				glVertex3dv(p.unpack());
				Point p_1 = p + .1 * bodyNormals[i];
				glVertex3dv(p_1.unpack());
			glEnd();
		}
	};

	void drawFans() {
		Vector v;
		for (int i = 0; i < numFans; i++) {
			glBegin(GL_TRIANGLE_FAN);
				v = fanNormals[i * (m_segmentsX + 1)];
				glNormal3dv(v.unpack());
				Point p = fans[i * (m_segmentsX + 1)];
				glVertex3dv(p.unpack());
				for (int j = 0; j < m_segmentsX; j++) {
					v = fanNormals[i * (m_segmentsX + 1) + j + 1];
					glNormal3dv(v.unpack());
					p = fans[i * (m_segmentsX + 1) + j + 1];
					glVertex3dv(p.unpack());
					//Vector v = getGLNormal(fanNormals[i * (m_segmentsX + 1) + j + 1], fanNormals[i * (m_segmentsX + 1) + 1 + (j+1) % m_segmentsX]);
					//glNormal3dv(v.unpack());
				}
				v = fanNormals[i * (m_segmentsX + 1)];
				glNormal3dv(v.unpack());
				p = fans[i * (m_segmentsX + 1) + 1];
				glVertex3dv(p.unpack());
			glEnd();
		}
	};

	void drawFanNormals() {
		glBegin(GL_LINES);
		for (int i = 0; i < numFans; i++) {
			float theta = 2 * PI / m_segmentsX;
			for (int j = 0; j < m_segmentsX; j++) {
				Point p = fans[i * (m_segmentsX + 1)];
				glVertex3dv(p.unpack());
				Point p_1 = p + .1 * fanNormals[i * (m_segmentsX + 1)];
				p_1 = rotY_mat(theta * j) * p_1;
				glVertex3dv(p_1.unpack());
			}
			for (int j = 1; j < m_segmentsX + 1; j++) {
				Point p = fans[i * (m_segmentsX + 1) + j];
				glVertex3dv(p.unpack());
				Point p_1 = p + .1 * fanNormals[i * (m_segmentsX + 1) + j];
				glVertex3dv(p_1.unpack());
			}
		}
		glEnd();
	};

	void drawSides() {
		int sidenum = (m_segmentsX+1) * (m_segmentsY+1);
		for (int s = 0; s < numSides; s++) {
			Vector v = bodyNormals[s];
			for (int col = 0; col < m_segmentsX; col++) {
				glBegin(GL_TRIANGLE_STRIP);
				glNormal3dv(v.unpack());
				for (int row = 0; row < m_segmentsY+1; row++) {
					Point p = bodyPoints[s * sidenum + col * (m_segmentsY + 1) + row];
					glVertex3dv(p.unpack());
					p = bodyPoints[s * sidenum + (col+1) * (m_segmentsY + 1) + row];
					glVertex3dv(p.unpack());
				}
				glEnd();
			}
		}
	};

	void drawSideNormals() {
		glBegin(GL_LINES);
		int sidecount = (m_segmentsX+1)*(m_segmentsY+1);
		for (int i = 0; i < numSides; i++) {
			for (int j = 0; j < sidecount; j++) {
				Point p = bodyPoints[i * sidecount + j];
				glVertex3dv(p.unpack());
				Point p_1 = p + .1 * bodyNormals[i];
				glVertex3dv(p_1.unpack());
			}
		}
		glEnd();
	};

	void drawTopNormals() {
		glBegin(GL_LINES);
		int sidecount = SQR((m_segmentsX+1));
		for (int i = 0; i < numFans; i++) {
			Vector v = fanNormals[i];
			for (int j = 0; j < sidecount; j++) {
				Point p = fans[i * sidecount + j];
				glVertex3dv(p.unpack());
				Point p_1 = p + .1 * v;
				glVertex3dv(p_1.unpack());
			}
		}
		glEnd();
	};

	void drawTop() {
		for (int s = 0; s < numFans; s++) {
			Vector v = fanNormals[s];
			for (int col = 0; col < m_segmentsX; col++) {
				glBegin(GL_TRIANGLE_STRIP);
				glNormal3dv(v.unpack());
				for (int row = 0; row < m_segmentsX+1; row++) {
					Point p = fans[s * SQR((m_segmentsX+1)) + col * (m_segmentsX + 1) + row];
					glVertex3dv(p.unpack());
					p = fans[s * SQR((m_segmentsX+1)) + (col+1) * (m_segmentsX + 1) + row];
					glVertex3dv(p.unpack());
				}
				glEnd();
			}
		}
	};

	Vector getGLNormal(Vector v1, Vector v2) {
		Vector v = v1 + v2;
		Vector x = normalize(v);
		return normalize(v);
	};
};

#endif
