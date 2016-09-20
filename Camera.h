

#ifndef CAMERA_H
#define CAMERA_H

#include "Algebra.h"

class Camera {
public:
Camera();
~Camera();
void Orient(Point eye, Point focus, Vector up);
void Orient(Point eye, Vector look, Vector up);
void SetViewAngle (double viewAngle);
void SetNearPlane (double nearPlane);
void SetFarPlane (double farPlane);
void SetScreenSize (int screenWidth, int screenHeight);

Matrix GetProjectionMatrix();
Matrix GetModelViewMatrix();

Matrix GetProjectionMatrixOrthographic();


void RotateV(double angle);
void RotateU(double angle);
void RotateW(double angle);
void Rotate(Point p, Vector axis, double degree);
void Translate(const Vector &v);

Point GetEyePoint();
Vector GetLookVector();
Vector GetUpVector();
Vector GetUVector();
Vector GetVVector();
double GetViewAngle();
double GetNearPlane();
double GetFarPlane();
int GetScreenWidth();
int GetScreenHeight();

double GetFilmPlanDepth();
double GetScreenWidthRatio();

void Reset();

private:
    double viewAngle;
    double nearPlane;
    double farPlane;
    float screenWidth, screenHeight;
    float pos_x, pos_y, pos_z;
    Point eye;
    Point focus;
    Vector look_vec;
    Vector   up_vec;
    Vector    u_vec;
    Vector    v_vec;
    Vector    w_vec;

    //Matrices
    Matrix projviewMatrix;
    Matrix modelviewMatrix;
    Matrix orthoProjviewMatrix;
};
#endif

