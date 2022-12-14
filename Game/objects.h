#pragma once
#include <glm/glm.hpp>
using glm::vec3;
using glm::vec4;

enum MaterialType {
	Normal, Reflective, Transparent
};

class Sphere
{
public:
	Sphere(vec3 center, float radius, vec3 color, float shininess, MaterialType matType) :
			center(center), radius(radius), color(color), shininess(shininess), matType(matType) {}

	vec3 center;
	float radius;
	vec3 color;
	float shininess;
	MaterialType matType;
};

class Plane
{
public:
	Plane(vec4 params, vec3 color, float shininess, MaterialType matType) :
		params(params), color(color), shininess(shininess), matType(matType) {}

	vec4 params;
	vec3 color;
	float shininess;
	MaterialType matType;
};

class DirLight
{
public:
	DirLight(vec3 direction, vec3 intensity) :
		direction(direction), intensity(intensity) {}

	vec3 direction;
	vec3 intensity;
};

class SpotLight
{
public:
	SpotLight(vec3 position, vec3 direction, vec3 intensity, float angle) :
		position(position), direction(direction), intensity(intensity), angle(angle) {}

	vec3 position;
	vec3 direction;
	vec3 intensity;
    float angle;
};