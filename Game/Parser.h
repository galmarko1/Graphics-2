#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "objects.h"

using std::ifstream;
using std::vector;
using std::string;
using std::stringstream;

void parse_scene(vec4& camera, vec3& ambient, vector<Sphere>& spheres, vector<Plane>& planes,
				 vector<DirLight>& dirLights, vector<SpotLight>& spotLights)
{
	
	ifstream file("C:\\Users\\itayr\\Documents\\scene2.txt");
	string lineStr = "";

	string objectType = "";
	int colorIndex = 0;
	int positionIndex = 0;
	int intensityIndex = 0;
	int sphereIndex = 0;
	int planeIndex = 0;
	int dirLightIndex = 0;
	int spotLightIndex = 0;
	vector<string> shapesTypes;
	vector<string> lightsTypes;
	
	while (getline(file, lineStr))
	{
		stringstream line(lineStr);
		string object;
		line >> object;

		if (object == "e") // camera
		{
			line >> camera.x;	line >> camera.y;	line >> camera.z;	line >> camera.w;
		}
		else if (object == "a") // camera
		{
			line >> ambient.r;	line >> ambient.g;	line >> ambient.b;
		}
		else if (object == "o" || object == "r" || object == "t") // shape
		{
			MaterialType matType = Normal;

			if (object == "o") matType = Normal;
			else if (object == "r") matType = Reflective;
			else if (object == "t") matType = Transparent;

			float x, y, z, w;
			line >> x;	line >> y;	line >> z;	line >> w;

			if (w > 0)
			{
				vec3 center(x, y, z);
				float radius = w;
				spheres.push_back(Sphere(center, radius, vec3(0, 0, 0), 0, matType));
				shapesTypes.push_back("sphere");
			}
			else
			{
				planes.push_back(Plane(vec4(x, y, z, w), vec3(0, 0, 0), 0, matType));
				shapesTypes.push_back("plane");
			}
		}
		else if (object == "c")
		{
			float x, y, z, w;
			line >> x;	line >> y;	line >> z;	line >> w;

			if (shapesTypes[colorIndex] == "sphere")
			{
				spheres[sphereIndex].color = vec3(x, y, z);
				spheres[sphereIndex].shininess = w;
				sphereIndex++;
			}
			else
			{
				planes[planeIndex].color = vec3(x, y, z);
				planes[planeIndex].shininess = w;
				planeIndex++;
			}

			colorIndex++;
		}
		else if (object == "d")
		{
			float x, y, z, w;
			line >> x;	line >> y;	line >> z;	line >> w;

			if (w == 0)
			{
				dirLights.push_back(DirLight(vec3(x, y, z), vec3(0)));
				lightsTypes.push_back("dir");
			}
			else
			{
				spotLights.push_back(SpotLight(vec3(0), vec3(x,y,z), vec3(0), 0));
				lightsTypes.push_back("spot");
			}
		}
		else if (object == "p")
		{
			float x, y, z, w;
			line >> x;	line >> y;	line >> z;	line >> w;
			
			spotLights[positionIndex].position = vec3(x, y, z);
			spotLights[positionIndex].angle = w;
			
			positionIndex++;
		}
		else if (object == "i")
		{
			float x, y, z, w;
			line >> x;	line >> y;	line >> z;	line >> w;

			if (lightsTypes[intensityIndex] == "dir")
			{
				dirLights[dirLightIndex].intensity = vec3(x, y, z);
				dirLightIndex++;
			}
			else
			{
				spotLights[spotLightIndex].intensity = vec3(x, y, z);
				spotLightIndex++;
			}

			intensityIndex++;
		}
	}

	file.close();
}