#ifndef SPHERE_H
#define SPHERE_H

#include "Shape.h"
#include <iostream>

class Sphere : public Shape {
public:
	Sphere() {
	};
	~Sphere() {};

	double Intersect(Point eyePointP, Vector rayV, Matrix transformMatrix) {
		double a = dot(rayV, rayV);
		Vector dist = eyePointP - Point(0, 0, 0);
		double b = 2 * dot(rayV, dist);
		double c = dot(dist, dist) - SQR(.5);
		return quadraticFormula(a, b, c);
	};

	Vector findIsectNormal(Point eyePoint, Vector ray, double dist) {
		Point p = eyePoint + ray * dist;
		Vector v(p[0], p[1], p[2]);
		v.normalize();
		return v;
	};

	void getColor(Point eye, Vector ray, double time, float &red, float &green, float &blue, Image &image, int i, int j) {
		Point p = eye + ray * time;
		int axis = -1;
		int value = -1;
		float a, b, u, v;
		int pixel, s, t;
			a = p[2]/p[0];
			b = p[1] * -1;
			if (p[0] < 0) {
				u = (atan(a) + PI/2 + PI) / (2 * PI);
			}
			else {
				u = (atan(a) + PI/2) / (2 * PI);
			}
			v = b + .5;
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
		for (int i = 0; i < crossLength(); i++) {
			float y = .5 * cos(PI * (i+1) / m_segmentsY);
			float z = .5 * sin(PI * (i+1) / m_segmentsY);
			Point p(0, y, z);
			crossection[i] = p;
			Vector v(0, y, z);
			crossNormals[i] = normalize(v);
		}
	}

	void calculateFanLine() {
		clearFanLine();
		numFans = 2;
		fanLine = new Point[4];
		fanLineNormals = new Vector[4];

		fanLine[0] = Point(0, .5, 0);
		Vector v(0, .5, 0);
		fanLineNormals[0] = normalize(v);

		fanLine[1] = Point(0, .5 * cos(PI / m_segmentsY), .5 * sin(PI / m_segmentsY));
		v = Vector(0, .5 * cos(PI / m_segmentsY), .5 * sin(PI / m_segmentsY));
		fanLineNormals[1] = normalize(v);

		fanLine[2] = Point(0, -.5, 0);
		v = Vector(0, -.5, 0);
		fanLineNormals[2] = normalize(v);

		fanLine[3] = Point(0, -.5 * cos(PI / m_segmentsY), .5 * sin(PI / m_segmentsY));
		v = Vector(0, -.5 * cos(PI / m_segmentsY), .5 * sin(PI / m_segmentsY));
		fanLineNormals[3] = normalize(v);
	}

	int crossLength() {
		return m_segmentsY-1;
	}
};

#endif
