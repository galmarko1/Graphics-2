#include "game.h"
#include "Ray.h"
#include <iostream>
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

bool hit_sphere(vec3 center, double radius, Ray r) {
    vec3 oc = r.origin() - center;
    auto a = dot(r.direction(), r.direction());
    auto b = 2.0 * dot(oc, r.direction());
    auto c = dot(oc, oc) - radius*radius;
    auto discriminant = b*b - 4*a*c;
    return (discriminant > 0);
}


vec3 ray_color(Ray& r){
    if(hit_sphere(vec3(128, 128, -10), 15, r)){
        return vec3 (1.0, 0.0, 0.0);
    }
    vec3 unit_dir = normalize(r.direction());
    float t = 0.5 * unit_dir[1] + 1.0;
    float t2 = 1.0 - t;
    return t2 * vec3 (1.0,1.0,1.0) + t * vec3 (0.5, 0.7, 1.0);
//    return vec3(255.0, 0.0, 0.0);
}

unsigned char* part1(){
    float mat[256][256][3];
//    float ratio = 16.0/9.0;
//    int width = 400;
//    int height = width/ ratio;

    vec3 origin = vec3(128, 128, 10);

    for(int i = 1; i <= 256; i ++){
        for(int j = 1; j <= 256; j ++){
            vec3 dir = vec3( i - origin[0] , j - origin[1], -10);
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
