#include <iostream>
#include <algorithm>
#include <random>
#include "game.h"
#include "Ray.h"
#include "objects.h"
#include "Parser.h"

using namespace glm;

struct Hit
{
    float t = -1;
    vec3 pos;
    vec3 normal;
    int id = -1;
    int shape_type = -1; // 1 = sphere, 2 = plane
};

std::vector<Sphere> spheres;
std::vector<Plane> planes;
std::vector<DirLight> dirLights;
std::vector<SpotLight> spotLights;

glm::vec3 ambient;

int samples_per_pixel;

#define EPSILON 0.01f

Game::Game() : Scene() {

}

Game::Game(float angle, float relationWH, float near1, float far1) : Scene(angle, relationWH, near1, far1) {
}

inline double random_double() {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

vec3 ray_color(Ray& r, int depth);// forward deceleration

vec3 refract(vec3 ray, vec3 normal, float refract_ratio) {
    //    float cos_theta = min(dot(-ray, normal), 1.0f);
    //    float sin_theta = sqrt(1.0 - cos_theta * cos_theta);
    //
    //    bool cannot_refract = (refract_ratio * sin_theta) > 1.0f;
    //
    //    if (!cannot_refract) {
    //        vec3 r_out_perp = refract_ratio * (ray + cos_theta * normal);
    //        float magnitude = dot(r_out_perp, r_out_perp);
    //        vec3 r_out_parallel = -normal * (float) (sqrtf(
    //                abs(1.0 - pow(magnitude, 2.0f))));
    vec3 uv = normalize(ray);
    float dt = dot(uv, normalize(normal));
    if (refract_ratio == 1.5) {
        normal = -normal;
    }
    float discriminant = 1.0f - refract_ratio * refract_ratio * pow((1.0f - dt * dt), 1);
    if (discriminant > 0) {
        return refract_ratio * (uv - dt * normal) - normal * sqrt(discriminant);
    }
    else {
        vec3 reflect_ray = 2 * dot(normalize(ray), -normal) * -normal;
        return
            reflect_ray;
    }
}

float hit_sphere(vec3 center, double radius, Ray r) {
    vec3 oc = r.origin() - center;
    float a = dot(r.direction(), r.direction());
    float half_b = dot(oc, r.direction());
    float c = dot(oc, oc) - radius * radius;
    float discriminant = half_b * half_b - a * c;
    if (discriminant < 0) //no hit
        return -1;
    float sqrtd = sqrt(discriminant);

    float root = (-half_b - sqrtd) / a;
    if (root < 0.01) {
        root = (-half_b + sqrtd) / a;
        if (root < 0.01)
            return -1;
    }

    return root;
}


float hit_plane(vec4 params, Ray r) {
    float a = params.x, b = params.y, c = params.z, d = params.w;
    vec3 normal = normalize(vec3(a, b, c));
    float denominator = dot(normal, r.direction());
    if (abs(denominator) > 0.0001) {
        float t = -(dot(r.origin(), normal) + d) / denominator;
        if (t > 0.0001) {
            return t;
        }
    }
    return -1;
}

Hit minHit(Ray& r, float minT, float maxT)
{
    Hit hit;
    hit.t = maxT;

    for (int i = 0; i < spheres.size(); i++)
    {
        float t = hit_sphere(spheres[i].center, spheres[i].radius, r);

        if (t > minT && t < hit.t)
        {
            hit.t = t;
            hit.pos = r.at(hit.t);
            hit.id = i;
            hit.normal = normalize(hit.pos - spheres[hit.id].center);
            hit.shape_type = 1;
        }
    }

    for (int i = 0; i < planes.size(); i++)
    {
        float t = hit_plane(planes[i].params, r);

        if (t > minT && t < hit.t)
        {
            hit.t = t;
            hit.pos = r.at(hit.t);
            hit.id = i;
            hit.normal = -1.0f * normalize(vec3(planes[i].params.x, planes[i].params.y, planes[i].params.z));
            hit.shape_type = 2;
        }
    }
    
    return hit;
}

vec3 ray_color(Ray& r, int depth)
{
    // maximum recursion
    if (depth == 0) {
        return vec3(0, 0, 0);
    }
    
    Hit hit = minHit(r, 0, 10000);

    if (hit.id == -1) { //background
        return vec3(0, 0, 0);
    }

    float diffuseK = 1;
    float specularK = 0.7;
    
    vec3 color = hit.shape_type == 1 ? spheres[hit.id].color : planes[hit.id].color;
    float shininess = hit.shape_type == 1 ? spheres[hit.id].shininess : planes[hit.id].shininess;
    MaterialType matType = hit.shape_type == 1 ? spheres[hit.id].matType : planes[hit.id].matType;

    if (hit.shape_type == 2) // plane checkboard
    {
        float chess = floor(hit.pos.x * 1.5) + floor(hit.pos.z * 1.5);
        chess = 2 * ((chess * 0.5) - int(chess * 0.5));
        diffuseK = chess == 0 ? 1 : 0.5;
    }

    vec3 sum_lights = ambient * color;

    for (int i = 0; i < dirLights.size(); i++) {
        DirLight light = dirLights[i];
        vec3 lightDir = -normalize(light.direction);

        Ray rayToLight = Ray(hit.pos + hit.normal*EPSILON, lightDir);
        Hit rayToLightHit = minHit(rayToLight, EPSILON, 10000);
            
        if (rayToLightHit.id != -1 && rayToLightHit.shape_type == 1) // shadow from sphere
        {
            continue;
        }

        vec3 diffuse = (float)std::max((double)dot(lightDir, hit.normal), 0.0) * color * light.intensity * diffuseK;
        vec3 R = normalize( - lightDir + 2 * dot(lightDir, hit.normal) * hit.normal );
        float p = pow(max(dot(-r.direction(), R), 0.0f), shininess);
        vec3 specular = specularK * p * light.intensity;
        
        sum_lights += (diffuse + specular);
    }

    for (int i = 0; i < spotLights.size(); i++) {
        SpotLight light = spotLights[i];
        vec3 lightDir = normalize(light.position - hit.pos);
        vec3 spotlightDir = normalize(light.direction);

        Ray rayToLight = Ray(hit.pos + hit.normal * EPSILON, lightDir);
        Hit rayToLightHit = minHit(rayToLight, EPSILON, length(light.position - hit.pos));

        if (rayToLightHit.id != -1 && rayToLightHit.shape_type == 1) // shadow from sphere
        {
            continue;
        }

        float spotCosine = dot(spotlightDir, -lightDir);

        if (spotCosine >= light.angle)
        {
            vec3 diffuse = (float)std::max((double)dot(lightDir, hit.normal), 0.0) * color * light.intensity * diffuseK;
            vec3 R = normalize(-lightDir + 2 * dot(lightDir, hit.normal) * hit.normal);
            float p = pow(max(dot(-r.direction(), R), 0.0f), shininess);
            vec3 specular = specularK * p * light.intensity;

            sum_lights += (diffuse + specular);
        }
    }

    if (matType == Reflective) {
        Ray reflect_ray = Ray(hit.pos + hit.normal * EPSILON,
            r.direction() + 2 * dot(-r.direction(), hit.normal) * hit.normal);
        sum_lights = specularK * ray_color(reflect_ray, depth - 1); // multiply with specularK 0.7 like in class? not the same as images
    }
    else if (matType == Transparent)
    {
        if (hit.shape_type == 1) // sphere
        {
            vec3 refract_ray_in = refract(r.direction(), hit.normal, (1.0f / 1.5f));
            Ray ray = Ray(hit.pos, refract_ray_in);
            float outPoint = hit_sphere(spheres[hit.id].center, spheres[hit.id].radius, ray);
            vec3 newNormal = normalize(spheres[hit.id].center - ray.at(outPoint));
            vec3 refract_ray_out = refract(normalize(ray.direction()), newNormal, (1.5f / 1.0f));
            Ray refract_ray;
            refract_ray = Ray(ray.at(outPoint), refract_ray_out);
            sum_lights += 1.0f * ray_color(refract_ray, depth - 1);
        }
        else // plane
        {
            vec3 refract_ray_in = refract(r.direction(), hit.normal, (1.0f / 1.5f));
            Ray refRay(hit.pos, refract_ray_in);
            sum_lights = ray_color(refRay, depth - 1);
        }
    }

    sum_lights = clamp(sum_lights, vec3(0, 0, 0), vec3(1, 1, 1));
    return sum_lights;
}

unsigned char* part1() {

    float mat[256][256][3];

    vec4 camera(0);
    parse_scene(camera, ambient, spheres, planes, dirLights, spotLights);
    vec3 origin(camera.x, camera.y, camera.z);
    samples_per_pixel = 32;

    /*vec3 origin = vec3(0.0, 0.0, 4.0);
    ambient = vec3(0.1, 0.2, 0.3);
    samples_per_pixel = 8;*/

    // Scene from PDF
    /*spheres.push_back(Sphere(vec3(-0.7, -0.7, -2.0), 0.5, vec3(1, 0, 0), 10, Normal));
    spheres.push_back(Sphere(vec3(0.6, -0.5, -1.0), 0.5, vec3(0.6, 0.0, 0.8), 10, Normal));
    planes.push_back(Plane(vec4(0, -0.5, -1.0, -3.5), vec3(0.0, 1.0, 1.0), 10, Normal));
    dirLights.push_back(DirLight(vec3(0.0, 0.5, -1.0), vec3(0.7, 0.5, 0.0)));
    spotLights.push_back(SpotLight(vec3(2.0, 1.0, 3.0), vec3(0.5, 0.0, -1), vec3(0.2, 0.5, 0.7), 0.6));*/


    // Scene 1
    //planes.push_back(Plane(vec4(1, 0, -0.1, -3.0), vec3(0.6, 0.0, 0.8), 20, Reflective));
    //planes.push_back(Plane(vec4(0, 0, -1.0, -3.5), vec3(0.7, 0.7, 0), 10, Reflective)); // floor
    //planes.push_back(Plane(vec4(-1, 0, -0.1, -3.0), vec3(0.0, 0.9, 0.5), 15, Reflective));
    //planes.push_back(Plane(vec4(0, 1, -0.1, -3.0), vec3(0.0, 0.8, 0.8), 10, Reflective));
    //planes.push_back(Plane(vec4(0, -1, -0.1, -3.0), vec3(0.9, 0.0, 0.1), 15, Reflective));
    //spheres.push_back(Sphere(vec3(-0.7, 0.7, -2.0), 1, vec3(1, 0, 0), 15, Normal));
    //spheres.push_back(Sphere(vec3(0.8, -0.5, -1.0), 0.7, vec3(0.0, 1.0, 0.8), 10, Normal));
    //spotLights.push_back(SpotLight(vec3(0.0, 0.0, 0.0), vec3(0.0, 0.5, -1), vec3(0.3, 0.9, 0.2), 0.8));
    //spotLights.push_back(SpotLight(vec3(0.0, 0.0, 0.0), vec3(0.5, 0.0, -1), vec3(0.8, 0.5, 0.7), 0.9));
    //spotLights.push_back(SpotLight(vec3(-0.2, 0.0, 0.0), vec3(-0.4, -0.3, -1), vec3(0.8, 0.5, 0.7), 0.7));
    //dirLights.push_back(DirLight(vec3(0.3, 0.5, -1.0), vec3(0.7, 0.8, 0.3)));

    // Scene 3

    //planes.push_back(Plane(vec4(1, 0, -0.1, -3.0), vec3(0.6, 0.0, 0.8), 20, Normal));
    //planes.push_back(Plane(vec4(0, 0, -1.0, -3.5), vec3(0.7, 0.7, 0), 10, Normal));
    //planes.push_back(Plane(vec4(-1, 0, -0.1, -3.0), vec3(0.0, 0.9, 0.5), 15, Normal));
    //planes.push_back(Plane(vec4(0, 1, -0.1, -3.0), vec3(0.0, 0.8, 0.8), 10, Normal));
    //planes.push_back(Plane(vec4(0, -1, -0.1, -3.0), vec3(0.9, 0.0, 0.1), 15, Normal));
    //
    ////spotLights.push_back(SpotLight(vec3(0, 0, 4), vec3(0, 0, -1), vec3(1, 1, 1), 0.98));
    //
    //spotLights.push_back(SpotLight(vec3(0.0, 0.0, 0.0), vec3(0.0, 0.5, -1), vec3(0.3, 0.9, 0.2), 0.8));
    //spotLights.push_back(SpotLight(vec3(0.0, 0.0, 0.0), vec3(0.5, 0.0, -1), vec3(0.9, 0.5, 0.5), 0.9));
    //spotLights.push_back(SpotLight(vec3(-0.2, 0.0, 0.0), vec3(-0.4, -0.3, -1), vec3(0.3, 0.5, 0.9), 0.7));
    //

    // Scene 4
    /*spheres.push_back(Sphere(vec3(-0.0, -0.5, -1.0), 0.3, vec3(0, 1, 1), 10, Normal));
    spheres.push_back(Sphere(vec3(-1.7, -0.7, -3.0), 0.7, vec3(1, 0.0, 0), 20, Normal));
    spheres.push_back(Sphere(vec3(2.6, -1.5, -10.0), 2.5, vec3(0.6, 0, 0.8), 15, Normal));
    spheres.push_back(Sphere(vec3(1.3, 1.5, -7.0), 1.5, vec3(0.9, 0.0, 0.0), 10, Normal));
    spheres.push_back(Sphere(vec3(-0.6, -0.5, -5.0), 1.0, vec3(0.0, 0.0, 0.8), 10, Normal));
    planes.push_back(Plane(vec4(0, -1, -1, -8.5), vec3(0.7, 0.7, 0.0), 10, Normal));
    dirLights.push_back(DirLight(vec3(0.0, -0.7, -1.0), vec3(0.9, 0.5, 0.0)));
    spotLights.push_back(SpotLight(vec3(-2.0, -1.0, 3.0), vec3(1.5, 0.9, -1), vec3(0.2, 0.5, 0.7), 0.6));*/

    // Scene 5
    /*spheres.push_back(Sphere(vec3(0.6, 0.5, -1.0), 0.5, vec3(0.6, 0, 0.8), 10, Normal));
    spheres.push_back(Sphere(vec3(-0.7, 0.7, -2.0), 0.5, vec3(1, 0.0, 0), 10, Normal));;
    planes.push_back(Plane(vec4(0, -0.5, -1, -3.5), vec3(0.0, 1.0, 1.0), 10, Reflective));
    dirLights.push_back(DirLight(vec3(0.0, 0.5, -1.0), vec3(0.7, 0.5, 0.0)));
    spotLights.push_back(SpotLight(vec3(2.0, 1.0, 3.0), vec3(0.5, 0.0, -1), vec3(0.2, 0.5, 0.7), 0.6));*/

    //Scene Google reflect

    /*spheres.push_back(Sphere(vec3(0, -0.3, -4), 0.5, vec3(0, 0, 1.0), 10, Transparent));
    planes.push_back(Plane(vec4(0, -1, 0, -1), vec3(1.0, 1.0, 1.0), 200, Normal));
    dirLights.push_back(DirLight(vec3(0, -1, 0), vec3(1.0, 1.0, 1.0)));*/
    //spotLights.push_back(SpotLight(vec3(0, 0, -4), vec3(0, 0, -1), vec3(1, 1, 1), 0.96));

    for (int y = 0; y < 256; y++) {
        for (int x = 0; x < 256; x++) {
            vec3 pixelColor = vec3(0, 0, 0);
            for (int i = 0; i < samples_per_pixel; i++) {
                double u = (double(x) + random_double()) / 255;
                double v = (double(y) + random_double()) / 255;
                vec3 xDir = vec3(2 * u, 0, 0);
                vec3 yDir = vec3(0, 2 * v, 0);
                vec3 left_corner = vec3(-1, -1, 0);

                vec3 dir = left_corner + xDir + yDir - origin;
                Ray r = Ray(origin, dir);
                vec3 color = ray_color(r, 4);
                pixelColor += color;
            }
            
            vec3 avgPixelColor = clamp((1.0f/samples_per_pixel) * pixelColor, 0.0f, 1.0f);
 
            mat[y][x][0] = avgPixelColor[0];
            mat[y][x][1] = avgPixelColor[1];
            mat[y][x][2] = avgPixelColor[2];
        }
    }

    int index = 0;
    unsigned char* data = new unsigned char[256 * 256 * 4];

    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            data[index] = mat[i][j][0] * 255.0;
            data[index + 1] = mat[i][j][1] * 255.0;
            data[index + 2] = mat[i][j][2] * 255.0;
            data[index + 3] = 255; //A of RGBA
            index += 4;
        }
    }
    return data;
}

void Game::Init() {

    AddShader("../res/shaders/pickingShader");
    AddShader("../res/shaders/basicShader");

    unsigned char* rayTraceTexture = part1();
    AddTexture(256, 256, rayTraceTexture);

    AddShape(Plane, -1, TRIANGLES);

    pickedShape = 0;

    SetShapeTex(0, 0);
    MoveCamera(0, zTranslate, 10);
    pickedShape = -1;

    //ReadPixel(); //uncomment when you are reading from the z-buffer
}

void Game::Update(const glm::mat4& MVP, const glm::mat4& Model, const int shaderIndx) {
    Shader* s = shaders[shaderIndx];
    int r = ((pickedShape + 1) & 0x000000FF) >> 0;
    int g = ((pickedShape + 1) & 0x0000FF00) >> 8;
    int b = ((pickedShape + 1) & 0x00FF0000) >> 16;
    s->Bind();
    s->SetUniformMat4f("MVP", MVP);
    s->SetUniformMat4f("Normal", Model);
    s->SetUniform4f("lightDirection", 0.0f, 0.0f, -1.0f, 0.0f);
    if (shaderIndx == 0)
        s->SetUniform4f("lightColor", r / 255.0f, g / 255.0f, b / 255.0f, 1.0f);
    else
        s->SetUniform4f("lightColor", 0.7f, 0.8f, 0.1f, 1.0f);
    s->Unbind();
}

void Game::WhenRotate() {
}

void Game::WhenTranslate() {
}

void Game::Motion() {
    if (isActive) {
    }
}

Game::~Game(void) {
}
