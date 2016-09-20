#ifndef CUBE_H
#define CUBE_H

#include "Shape.h"
#include <iostream>

class Cube : public Shape {
public:
	Cube() {
		isRotate = false;
	};
	~Cube() {};

	double Intersect(Point eyePointP, Vector rayV, Matrix transformMatrix) {
		// test for intersection at each plane that makes up a cube
		double min_t = -1;
		for (int i = 0; i < 3; i++) {
			for (double j = -.5; j < 1; j += 1) {
				double t = planeIntersection(eyePointP, rayV, i, j);
				if (t > 0) {
					Point p = eyePointP + rayV * t;
					if ((p[(i+1)%3] <= width/2.0 && p[(i+1)%3] >= width/(-2.0) && p[(i+2)%3] <= width/2.0 && p[(i+2)%3] >= width/(-2.0)) && (min_t > t || min_t < 0) && t >= 0) {
						min_t = t;
					}
				}
			}
		}
		return min_t;
	};

	Vector findIsectNormal(Point eyePoint, Vector ray, double dist) {
		Point p = eyePoint + ray * dist;
		Vector v;
		for (int i = 0; i < 3; i++) {
			for (double j = -.5; j < 1; j += 1) {
				if (IN_RANGE(p[i], j)) {
					v[i] = j;
					v.normalize();
					return v;
				}
			}
		}
	};


	void getColor(Point eye, Vector ray, double time, float &red, float &green, float &blue, Image &image, int i, int j) {
		Point intersect = eye + ray * time;
		int axis = -1;
		int value = -1;
		float a, b, u, v;
		int pixel, s, t;
		for (int i = 0; i < 3; i++) {
			for (double j = -.5; j < 1; j += 1) {
				if (IN_RANGE(intersect[i], j)) {
					axis = i;
					value = j;
					break;
				}
			}
			if (axis != -1) {
				break;
			}
		}
		if (axis != 1) {
			b = intersect[1];
			a = (axis == 0) ? intersect[2] : intersect[0];
		}
		else {
			b = intersect[2];
			a = intersect[0];
		}
		b *= -1;
		u = a + 0.5;
		v = b + 0.5;
		s = (int)fmod((u * image.w * i), image.w);
		t = (int)fmod((v * image.h * j), image.h);
		pixel = t * image.w * 3 + s * 3;
		red = image.pixels[pixel];
		green = image.pixels[pixel + 1];
		blue = image.pixels[pixel + 2];
	};

protected:

	void calculateTop() {
		clearFans();
		numFans = 2;
		int oneside = SQR((m_segmentsX+1));
		fans = new Point[2 * oneside];
		for (int i = 0; i < m_segmentsX+1; i++) {
			for (int j = 0; j < m_segmentsX+1; j++) {
				fans[i*(m_segmentsX+1) + j] = Point(-.5 + (float)i/m_segmentsX, .5, -.5 + (float)j/m_segmentsX);
				fans[oneside + i*(m_segmentsX+1) + j] = Point(-.5 + (float)i/m_segmentsX, -.5, .5 - (float)j/m_segmentsX);
			}
		}
		fanNormals = new Vector[numFans];
		fanNormals[0] = Vector(0, 1, 0);
		fanNormals[1] = Vector(0, -1, 0);
	}

	void calculateSides() {
		clearBodyPoints();
		numSides = 4;
		int sidecount = (m_segmentsY + 1) * (m_segmentsX + 1);
		bodyPoints = new Point[sidecount * numSides];
		bodyNormals = new Vector[numSides];
		bodyNormals[0] = Vector(-.5, 0, 0);
		for (int i = 0; i < m_segmentsX + 1; i++) {
			for (int j = 0; j < m_segmentsY + 1; j++) {
				bodyPoints[i*(m_segmentsY+1) + j] = Point(-.5, .5 - (float)j/m_segmentsY, -.5 + (float)i/m_segmentsX);
			}
		}
		float theta = 2 * PI / numSides;
		for (int i = 1; i < numSides; i++) {
			Matrix r = rotY_mat(theta * i);
			bodyNormals[i] = r * bodyNormals[0];
			for (int j = 0; j < sidecount; j++) {
				Point p = bodyPoints[j];
				bodyPoints[i*sidecount + j] = r * p;
			}
		}
	}
};

#endif
