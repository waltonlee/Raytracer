#include "Camera.h"

Camera::Camera() {

}

Camera::~Camera() {
}

void Camera::Orient(Point eye, Point focus, Vector up) {
	this->eye = eye;
	this->pos_x = eye[0];
	this->pos_y = eye[1];
	this->pos_z = eye[2];
	this->focus = focus;
	this->up_vec = up;
	Vector look(focus[0] - eye[0],focus[1] - eye[1],focus[2] - eye[2]);
	this->w_vec  = -1 * look;
	this->w_vec.normalize();
	this->u_vec  = cross(up, w_vec);
	this->u_vec.normalize();
	this->v_vec  = cross(w_vec, u_vec);
	this->v_vec.normalize();
}


void Camera::Orient(Point eye, Vector look, Vector up) {
	this->eye = eye;
	this->pos_x = eye[0];
	this->pos_y = eye[1];
	this->pos_z = eye[2];
	this->look_vec = look;
	this->up_vec =  up;
	this->w_vec  = -1 * look;
	this->w_vec.normalize();
	this->u_vec  = cross(up, w_vec);
	this->u_vec.normalize();
	this->v_vec  = cross(w_vec, u_vec);
	this->v_vec.normalize();
}

/* */
Matrix Camera::GetProjectionMatrix() {

	//unhinge

	float c = -1 * nearPlane / farPlane;
	Matrix unhinge(   1, 0, 0, 0,
					  0, 1, 0, 0,
					  0, 0, -1/(c + 1), c/( c + 1),
					  0,0,-1,0);

	//scale
	Matrix normalize( 1/ (tan(viewAngle/2.0f) * farPlane * GetScreenWidthRatio()), 0, 0, 0,
					  0, 1 / (tan(viewAngle/2.0f) * farPlane), 0, 0,
					  0, 0, 1/farPlane, 0,
					  0,0,0,1);

	this->projviewMatrix = unhinge * normalize;
	return this->projviewMatrix;
}

Matrix Camera::GetProjectionMatrixOrthographic() {

	//unhinge
	//scale

	Matrix normalize( 2/screenWidth, 0, 0, 0,
					  0, 2/screenHeight, 0, 0,
					  0, 0, 1/farPlane, 0,
					  0,0,0,1);
	this->orthoProjviewMatrix = normalize;
	return normalize;
}


void Camera::SetViewAngle (double viewAngle) {
	this->viewAngle = DEG_TO_RAD(viewAngle)	;
}

void Camera::SetNearPlane (double nearPlane) {
	this->nearPlane = nearPlane;
}

void Camera::SetFarPlane (double farPlane) {
	this->farPlane = farPlane;
}

void Camera::SetScreenSize (int screenWidth, int screenHeight) {
	this->screenHeight = screenHeight;
	this->screenWidth = screenWidth;
}


Matrix Camera::GetModelViewMatrix() {

	//rotate
	Matrix m_rotate(    u_vec[0],    u_vec[1],    u_vec[2], 0,
					    v_vec[0],    v_vec[1],    v_vec[2], 0,
					    w_vec[0],    w_vec[1],    w_vec[2], 0,
					           0,           0,           0, 1);

	//translate
	Vector modelView_translate( -pos_x, -pos_y, -pos_z);
	Matrix m_translate = trans_mat( modelView_translate);

	this->modelviewMatrix = m_rotate * m_translate;

	return this->modelviewMatrix;
}

void Camera::RotateV(double angle) {
	// yaw
	angle = DEG_TO_RAD(angle);
	Matrix rot_v = rot_mat(v_vec, angle);
	w_vec = rot_v * w_vec;
	u_vec = rot_v * u_vec;

}

void Camera::RotateU(double angle) {
	// pitch
	angle = DEG_TO_RAD(angle);
	Matrix rot_u = rot_mat(u_vec, angle);
	w_vec = rot_u * w_vec;
	v_vec = rot_u * v_vec;
}

void Camera::RotateW(double angle) {
	// roll
	angle = DEG_TO_RAD(angle);
	Matrix rot_w = rot_mat(w_vec, -angle);
	v_vec = rot_w * v_vec;
	u_vec = rot_w * u_vec;
}

void Camera::Translate(const Vector &v) {
	//change x y z
	eye = eye + v;
	pos_x = v[0];
	pos_y = v[1];
	pos_z = v[2];
}


void Camera::Rotate(Point p, Vector axis, double degrees) {

	Matrix rot_vec = rot_mat(p, axis, DEG_TO_RAD(degrees));
	v_vec = rot_vec * v_vec;
	u_vec = rot_vec * u_vec;
	w_vec = rot_vec * w_vec;
}


Point Camera::GetEyePoint() {
	return eye;
}

Vector Camera::GetLookVector() {
	return -1 * w_vec;
}

Vector Camera::GetUpVector() {
	return v_vec;
}

Vector Camera::GetUVector(){
	return u_vec;
}

Vector Camera::GetVVector(){
	return v_vec;
}

double Camera::GetViewAngle() {
	return RAD_TO_DEG(viewAngle);
}

double Camera::GetNearPlane() {
	return nearPlane;
}

double Camera::GetFarPlane() {
	return farPlane;
}

int Camera::GetScreenWidth() {
	return screenWidth;
}

int Camera::GetScreenHeight() {
	return screenHeight;
}

double Camera::GetFilmPlanDepth() {
	return nearPlane;
}

double Camera::GetScreenWidthRatio() {
	return (double)screenWidth/(double)screenHeight;
}

void Camera::Reset()
{
	// this->eye = Point(2,2,2);
	// this->pos_x = eye[0];
	// this->pos_y = eye[1];
	// this->pos_z = eye[2];
	// this->focus = Point(-2,-2,-2);
	// this->up_vec = Vector(0,1,0);
	// Vector look(focus[0] - eye[0],focus[1] - eye[1],focus[2] - eye[2]);
	// this->w_vec  = -1 * look;
	// this->w_vec.normalize();
	// this->u_vec  = cross(this->up_vec, w_vec);
	// this->u_vec.normalize();
	// this->v_vec  = cross(w_vec, u_vec);
	// this->v_vec.normalize();
	this->nearPlane = 0.3f;
	this->farPlane = 100.0f;
	this->viewAngle = 45.0f;
	// this->screenWidth = 640.0f;
	// this->screenHeight = 480.0f;
}
