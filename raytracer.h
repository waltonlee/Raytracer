#ifndef RAYTRACER_H
#define RAYTRACER_H

#include "Camera.h"
#include <iostream>

Vector generateRay(int x, int y, int w, int h, Point eye, Camera *camera) {
	Point P(eye);
	Point Q = P + camera->GetNearPlane() * camera->GetLookVector();

	float theta = DEG_TO_RAD(camera->GetViewAngle());
	float H = camera->GetNearPlane() * tan(theta/2.0);
	float W = H * camera->GetScreenWidthRatio();
	float a = -W + 2 * W * x/w;
	float b = -H + 2 * H * y/h;

	Point S = Q + a * camera->GetUVector() + b * camera->GetVVector();
	Vector ray = S - P;
	ray.normalize();
	return ray;
}

#endif
