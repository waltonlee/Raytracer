#ifndef CYLINDER_H
#define CYLINDER_H

#include "Shape.h"

class Cylinder : public Shape {
public:
	Cylinder() {};
	~Cylinder() {};

	double Intersect(Point eyePointP, Vector rayV, Matrix transformMatrix) {
		double min_t = -1;
		double t;
		for (double i = -.5; i < 1; i += 1) {
			t = planeIntersection(eyePointP, rayV, 1, i);
			if (t > 0) {
				Point p = eyePointP + rayV * t;
				if ((SQR(p[0]) + SQR(p[2]) < SQR(width/2.0)) && t >= 0 && (min_t < 0 || min_t > t)) {
					min_t = t;
				}
			}
		}
		double a = SQR(rayV[0]) + SQR(rayV[2]);
		double b = 2 * (eyePointP[0] * rayV[0] + eyePointP[2] * rayV[2]);
		double c = SQR(eyePointP[0]) + SQR(eyePointP[2]) - SQR(.5);
		double result = quadraticFormula(a, b, c);
		if (result == -1) {
			return min_t;
		}
		Point p = eyePointP + rayV * result;
		if (p[1] > width/2.0 || p[1] < -1 * width/2.0) {
			result = -1;
		}
		return ((result < min_t || min_t < 0) && result >= 0) ? result : min_t;
	};

	virtual Vector findIsectNormal(Point eyePoint, Vector ray, double dist) {
		Point p = eyePoint + ray * dist;
		if (IN_RANGE(p[1], width/2.0) && eyePoint[1] > width/2.0) {
			return Vector(0, 1, 0);
		}
		else if (IN_RANGE(p[1], width/-2.0) && eyePoint[1] < width/-2.0) {
			return Vector(0, -1, 0);
		}
		else {
			Vector v(p[0], 0, p[2]);
			v.normalize();
			return v;
		}
	};

	void getColor(Point eye, Vector ray, double time, float &red, float &green, float &blue, Image &image, int i, int j) {
		Point p = eye + ray * time;
		int axis = -1;
		int value = -1;
		float a, b, u, v;
		int pixel, s, t;
		if (IN_RANGE(p[1], width/2.0) && eye[1] > width/2.0 || IN_RANGE(p[1], width/-2.0) && eye[1] < width/-2.0) {
			a = p[0];
			b = p[2];
			u = a + .5;
			v = b * -1 + .5;
		}
		else {
			a = p[2]/p[0];
			b = p[1] * -1;
			if (p[0] < 0) {
				u = (atan(a) + PI/2 + PI) / (2 * PI);
			}
			else {
				u = (atan(a) + PI/2) / (2 * PI);
			}
			v = b + .5;
		}
		s = (int)fmod((u * image.w * i), image.w);
		t = (int)fmod((v * image.h * j), image.h);
		pixel = t * image.w * 3 + s * 3;
		red = image.pixels[pixel];
		green = image.pixels[pixel + 1];
		blue = image.pixels[pixel + 2];
	};

protected:
	void calculateCrossection() {
		clearCrossection();
		crossection = new Point[crossLength()];
		crossNormals = new Vector[crossLength()];
		Vector v(0, 0, 1);
		for (int i = 0; i < crossLength(); i++) {
			Point p(0,
				((float)m_segmentsY-i)/m_segmentsY-.5,
				 .5);
			crossection[i] = p;
			crossNormals[i] = v;
		}
	}

	int crossLength() {
		return m_segmentsY+1;
	}

	void calculateFanLine() {
		clearFanLine();
		numFans = 2;
		fanLine = new Point[4];
		fanLineNormals = new Vector[4];
		Vector v(0, 1, 0);

		fanLine[0] = Point(0, .5, 0);
		fanLineNormals[0] = v;

		fanLine[1] = Point(0, .5, .5);
		fanLineNormals[1] = v;

		fanLine[2] = Point(0, -.5, 0);
		v = Vector(0, -1, 0);
		fanLineNormals[2] = v;

		fanLine[3] = Point(0, -.5, .5);
		fanLineNormals[3] = v;
	}
};
#endif
