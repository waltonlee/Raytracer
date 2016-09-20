#ifndef CONE_H
#define CONE_H

#include "Shape.h"

class Cone : public Shape {
public:
	Cone() {};
	~Cone() {};

	double Intersect(Point eyePointP, Vector rayV, Matrix transformMatrix) {
		double min_t = -1;
		double t;
		t = planeIntersection(eyePointP, rayV, 1, -.5);
		if (t > 0) {
			Point p = eyePointP + rayV * t;
			if ((SQR(p[0]) + SQR(p[2]) < SQR(width/2.0)) && t >= 0 && (min_t < 0 || min_t > t)) {
				min_t = t;
			}
		}

		double yk = eyePointP[1] - 0.5;
		double k = width/2.0;
		double a = SQR(rayV[0]) + SQR(rayV[2]) - k * k * SQR(rayV[1]);
		double b = 2 * (eyePointP[0] * rayV[0] + eyePointP[2] * rayV[2] - k * k * yk * rayV[1]);
		double c = SQR(eyePointP[0]) + SQR(eyePointP[2]) - k * k * yk * yk;
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

	Vector findIsectNormal(Point eyePoint, Vector ray, double dist) {
		Point p = eyePoint + ray * dist;
		if (IN_RANGE(p[1], width/-2.0) && eyePoint[1] < width/-2.0) {
			return Vector(0, -1, 0);
		}
		Vector ho(p[0], 0, p[2]);
		ho.normalize();
		ho = ho * 2.0;
		Vector vert(0, 1, 0);
		Vector v = ho + vert;
		v.normalize();
		return v;
	};

	void getColor(Point eye, Vector ray, double time, float &red, float &green, float &blue, Image &image, int i, int j) {
		Point p = eye + ray * time;
		int axis = -1;
		int value = -1;
		float a, b, u, v;
		int pixel, s, t;
		if (IN_RANGE(p[1], width/-2.0) && eye[1] < width/-2.0) {
			u = .5 * sqrt(2 + SQR(a) - SQR(b) + 2 * a * sqrt(2)) - .5 * sqrt(2 + SQR(a) - SQR(b) - 2 * a * sqrt(2));
			v = .5 * sqrt(2 - SQR(a) + SQR(b) + 2 * a * sqrt(2)) - .5 * sqrt(2 - SQR(a) + SQR(b) - 2 * a * sqrt(2));
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
		/* slope is -2 so -1/slope */
		float inverse_slope = .5;
		Vector v(0, 1, 2);
		v = normalize(v);
		for (int i = 0; i < m_segmentsY; i++) {
			/* cone has fan for top component */
			Point p(0,
				((float)m_segmentsY-(i+1))/m_segmentsY-.5,
				.5*(i+1)/((float)m_segmentsY));
			crossection[i] = p;
			crossNormals[i] = v;
		}
	}

	int crossLength() {
		return m_segmentsY;
	}

	void calculateFanLine() {
		clearFanLine();
		numFans = 2;
		fanLine = new Point[4];
		fanLineNormals = new Vector[4];
		Vector v(0, 1, 2);
		v = normalize(v);

		fanLine[0] = Point(0, .5, 0);
		fanLineNormals[0] = v;

		fanLine[1] = Point(0, .5 - 1.0/m_segmentsY, .5 * 1.0/m_segmentsY);
		fanLineNormals[1] = v;

		fanLine[2] = Point(0, -.5, 0);
		v = Vector(0, -1, 0);
		fanLineNormals[2] = normalize(v);

		fanLine[3] = Point(0, -.5, .5);
		v = Vector(0, -1, 0);
		fanLineNormals[3] = normalize(v);
	}
};

#endif
