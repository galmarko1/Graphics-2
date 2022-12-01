#include "game.h"
#include "Ray.h"
#include <iostream>
#include <algorithm>
#include <glm/gtc/matrix_transform.hpp>

using namespace glm;

static void printMat(const glm::mat4 mat)
{
	std::cout<<" matrix:"<<std::endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			std::cout<< mat[j][i]<<" ";
		std::cout<<std::endl;
	}
}

Game::Game() : Scene()
{
}

Game::Game(float angle ,float relationWH, float near1, float far1) : Scene(angle,relationWH,near1,far1)
{ 	
}

float hit_sphere(vec3 center, double radius, Ray r) {
    vec3 oc = r.origin() - center;
    float a = dot(r.direction(), r.direction());
    float b = 2.0 * dot(oc, r.direction());
    float c = dot(oc, oc) - radius*radius;
    float discriminant = b*b - 4*a*c;
    if (discriminant < 0) //no hit
        return -1;
    return (-b -sqrt(discriminant))/(2.0 * a);
}


vec3 ray_color(Ray& r){
    vec3 center = vec3(0, 0, -1); // sphere center
    vec3 sphere_color = vec3(1,0,0);
    vec3 light_dir = normalize(vec3(0, 15, -3)); // light dir
    vec3 light_color = vec3(1, 1, 1);
    float diffuseK = 1;
    float specularK = 0.7;

    float t = hit_sphere(center, 0.5, r);
    if(t > 0)
    {
        vec3 normal = normalize(r.at(t) - center) ;
        vec3 R = light_dir - 2 * dot(light_dir, normal) * normal;
        float p = pow(max(dot(normalize(-1.0f * r.direction()), normalize(R)), 0.0f) , 256);
        vec3 specular = vec3 (0,0,0);
        specular = specularK * p * light_color ;

        vec3 diffuse = (float) std::max((double) dot(normalize(light_dir), normal), 0.0) * sphere_color * diffuseK;
        vec3 color = diffuse + specular;
        return color;
    }
    //background
    vec3 unit_dir = normalize(r.direction());
    t = 0.5 * unit_dir[1] + 1.0;
    float t2 = 1.0 - t;
    return t2 * vec3 (1.0,1.0,1.0) + t * vec3 (0.5, 0.7, 1.0);
}

unsigned char* part1(){
    float mat[256][256][3];
//    float ratio = 16.0/9.0;
//    int width = 400;
//    int height = width/ ratio;

    vec3 origin = vec3(0, 0, 4.0);

    for(int i = 0; i <256 ; i ++){
        for(int j = 0; j < 256; j ++){
            double u = double(j)/ 255;
            double v = double(i)/ 255;
            vec3 xDir = vec3(2 * u,0,0);
            vec3 yDir = vec3(0, 2 * v, 0);
            vec3 left_corner = vec3(-1, -1, 0);

            vec3 dir = left_corner + xDir + yDir - origin;
            Ray r = Ray(origin, dir);
            vec3 color = ray_color(r);
            mat[i-1][j-1][0] = color[0];
            mat[i-1][j-1][1] = color[1];
            mat[i-1][j-1][2] = color[2];

        }
    }

    int index = 0;
    unsigned char *data = new unsigned char [256 * 256 * 4];

    for (int i =0; i< 256; i++){
        for(int j =0; j < 256; j++){
            data[index] = mat[i][j][0] * 255.0;
            data[index + 1] = mat[i][j][1] * 255.0;
            data[index + 2] = mat[i][j][2] * 255.0;
            data[index + 3] = 255; //A of RGBA
            index += 4;
        }
    }
    return data;
}

void Game::Init()
{		

	AddShader("../res/shaders/pickingShader");	
	AddShader("../res/shaders/basicShader");

//    part1();

//	AddTexture("../res/textures/box0.bmp",false);

    unsigned char *zibi = part1();
    AddTexture(256,256,zibi);

	AddShape(Plane,-1,TRIANGLES);
	
	pickedShape = 0;
	
	SetShapeTex(0,0);
	MoveCamera(0,zTranslate,10);
	pickedShape = -1;
	
	//ReadPixel(); //uncomment when you are reading from the z-buffer
}

void Game::Update(const glm::mat4 &MVP,const glm::mat4 &Model,const int  shaderIndx)
{
	Shader *s = shaders[shaderIndx];
	int r = ((pickedShape+1) & 0x000000FF) >>  0;
	int g = ((pickedShape+1) & 0x0000FF00) >>  8;
	int b = ((pickedShape+1) & 0x00FF0000) >> 16;
	s->Bind();
	s->SetUniformMat4f("MVP", MVP);
	s->SetUniformMat4f("Normal",Model);
	s->SetUniform4f("lightDirection", 0.0f , 0.0f, -1.0f, 0.0f);
	if(shaderIndx == 0)
		s->SetUniform4f("lightColor",r/255.0f, g/255.0f, b/255.0f,1.0f);
	else 
		s->SetUniform4f("lightColor",0.7f,0.8f,0.1f,1.0f);
	s->Unbind();
}

void Game::WhenRotate()
{
}

void Game::WhenTranslate()
{
}

void Game::Motion()
{
	if(isActive)
	{
	}
}

Game::~Game(void)
{
}
