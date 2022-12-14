#include <iostream>
#include <algorithm>
#include <random>
#include "game.h"
#include "Ray.h"
#include "objects.h"

using namespace glm;

std::vector<Sphere> spheres;
std::vector<Plane> planes;
std::vector<DirLight> dirLights;
std::vector<SpotLight> spotLights;

glm::vec3 ambient;

int samples_per_pixel;

static void printMat(const glm::mat4 mat) {
    std::cout << " matrix:" << std::endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++)
            std::cout << mat[j][i] << " ";
        std::cout << std::endl;
    }
}

Game::Game() : Scene() {

}

Game::Game(float angle, float relationWH, float near1, float far1) : Scene(angle, relationWH, near1, far1) {
}

inline double random_double() {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

vec3 ray_color(Ray& r, int depth, int planeOrSphere, int objectIndex);// forward deceleration

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
    vec3 center = vec3(-a * d, -b * d, -c * d);
    float denominator = dot(normal, r.direction());
    if (abs(denominator) > 0.0001) {
        vec3 diff = center - r.origin();
        float t = dot(diff, normal) / denominator;
        if (t > 0.0001) {
            return t;
        }
    }
    return -1;
}

std::pair<int, int> minHit(Ray& r, int planeOrSphere, int objectIndex) {
    float min_t_sphere = 100000;
    float min_t_plane = 100000;
    int min_sphere_index = -1;
    int min_plane_index = -1;
    int i = 0;

    for (Sphere sphere : spheres) {
        if (planeOrSphere == 1) {
            if (objectIndex == i) {
                i++;
                continue;
            }
        }
        float curr_t = hit_sphere(sphere.center, sphere.radius, r);
        if (curr_t > -1) {
            if (min_t_sphere > curr_t) {
                min_sphere_index = i;
                min_t_sphere = curr_t;
            }
        }
        i++;
    }

    i = 0;
    for (Plane plane : planes) {
        if (planeOrSphere == 2) {
            if (objectIndex == i) {
                i++;
                continue;
            }
        }
        float curr_t = hit_plane(plane.params, r);
        if (curr_t > -1) {
            if (min_t_plane > curr_t) {
                min_plane_index = i;
                min_t_plane = curr_t;
            }
        }
        i++;
    }
    std::pair<int, int> p;
    if (min_plane_index == -1 && min_sphere_index == -1) {
        p.first = -1;
        p.second = 0;
    }
    else if (min_t_plane < min_t_sphere) {
        p.first = min_plane_index;
        p.second = 2; // plane vector
    }
    else {
        p.first = min_sphere_index;
        p.second = 1; // sphere vector
    }
    return p;
}

vec3 ray_color_plane(Ray& r, int index, int depth) {
    Plane plane = planes[index];
    vec3 plane_color = plane.color;

    float diffuseK = 1;
    float specularK = 0.7;
    float spot_factor = 1.0;
    float t = hit_plane(plane.params, r);

    if (t > 0) {
        vec3 normal = -1.0f * normalize(vec3(plane.params.x, plane.params.y, plane.params.z));
        vec3 sum_lights = ambient * plane.color;
        vec3 diffuse_and_specular = vec3(0, 0, 0);
        vec3 coordinate = r.at(t);

        float chess = floor(coordinate.x) + floor(coordinate.y);
        chess = 2 * ((chess * 0.5) - int(chess * 0.5));
        diffuseK = chess == 0 ? 1 : 0.5;

        for (int i = 0; i < dirLights.size(); i++) {
            DirLight light = dirLights[i];
            Ray ray = Ray(coordinate, -1.0f * light.direction);
            std::pair<int, int> hit = minHit(ray, 2, index);
            if (hit.second == 1) {
                sum_lights += vec3(0, 0, 0);
                continue;
            }
            
            vec3 R = light.direction - 2 * dot(-light.direction, normal) * normal;
            float p = pow(max(dot(normalize(-1.0f * r.direction()), normalize(R)), 0.0f), plane.shininess);
            vec3 specular = vec3(0, 0, 0);
            specular = specularK * p * light.intensity;
            vec3 diffuse =
                (float)std::max((double)dot(normalize(-light.direction), normal), 0.0) * plane.color * light.intensity *
                diffuseK;
            diffuse_and_specular += diffuse + specular;
            sum_lights += diffuse_and_specular;
        }
        for (int i = 0; i < spotLights.size(); i++) {
            SpotLight light = spotLights[i];
            vec3 spot_light_dir = coordinate - light.position; // spotlight dir
            Ray ray = Ray(coordinate, -1.0f * spot_light_dir);
            std::pair<int, int> hit = minHit(ray, 2, index);
            if (hit.second == 1) {
                sum_lights += vec3(0, 0, 0);
                continue;
            }
            
            vec3 L = normalize(light.direction);

            vec3 D = normalize(r.at(t) - light.position);
            float spotCosine = dot(D, L);

            if (spotCosine >= light.angle) {
                vec3 R = 1.0f * spot_light_dir - 2 * dot(-1.0f * spot_light_dir, normal) * normal;
                float p = pow(max(dot(normalize(-1.0f * r.direction()), normalize(R)), 0.0f), plane.shininess);
                vec3 specular = vec3(0, 0, 0);
                specular = specularK * p * light.intensity;
                vec3 diffuse = light.intensity * plane.color *
                    (float)std::max((double)dot(normalize(1.0f * spot_light_dir), normal), 0.0) *
                    diffuseK;
                sum_lights += specular + diffuse;
            }
        }
        if (plane.matType == Reflective) {
            Ray reflect_ray = Ray(r.at(t),
                normalize(r.direction()) - 2 * dot(normalize(r.direction()), -normal) * -normal);
            sum_lights += specularK * ray_color(reflect_ray, depth - 1, 1, index);
        }

        if (plane.matType == Transparent) {
            vec3 refract_ray_in = refract(normalize(r.direction()), normal, (1.0f / 1.5f));
            Ray ray = Ray(r.at(t), refract_ray_in);
            sum_lights += 1.0f * ray_color(ray, depth - 1, 2, index);
        }


        plane_color = sum_lights;

        if (plane_color.x > 1.0f) {
            plane_color = vec3(1.0, plane_color.y, plane_color.z);
        }
        if (plane_color.y > 1.0f) {
            plane_color = vec3(plane_color.x, 1.0, plane_color.z);
        }
        if (plane_color.z > 1.0f) {
            plane_color = vec3(plane_color.x, plane_color.y, 1.0f);
        }
        return plane_color;
    }
    //background
    vec3 unit_dir = normalize(r.direction());
    t = 0.5 * unit_dir[1] + 1.0;
    float t2 = 1.0 - t;
    return t2 * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
}

vec3 ray_color(Ray& r, int depth, int planeOrSphere, int objectIndex) {

    // maximum recursion
    if (depth == 0) {
        return vec3(0, 0, 0);
    }

    float diffuseK = 1;
    float specularK = 0.7;
    std::pair<int, int> hit = minHit(r, planeOrSphere, objectIndex);

    if (hit.first == -1) { //background
        //        return vec3(0, 0, 0);
        vec3 unit_dir = normalize(r.direction());
        float t = 0.5 * unit_dir[1] + 1.0;
        float t2 = 1.0 - t;
        return t2 * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
    }
    if (hit.second == 1) { //sphere
        Sphere sphere = spheres[hit.first];
        float t_sphere = hit_sphere(sphere.center, sphere.radius, r);
        vec3 normal = normalize(r.at(t_sphere) - sphere.center);
        vec3 sum_lights = ambient * sphere.color;
        vec3 coordinate = r.at(t_sphere);

        for (int i = 0; i < dirLights.size(); i++) {
            DirLight light = dirLights[i];
            Ray ray = Ray(coordinate, -1.0f * light.direction);
            std::pair<int, int> sphereHit = minHit(ray, 1, hit.first);
            if (sphereHit.second == 1) {
                sum_lights += vec3(0, 0, 0);
                continue;
            }
            vec3 R = 1.0f * light.direction - 2 * dot(light.direction, normal) * normal;
            float p = pow(max(dot(normalize(-1.0f * r.direction()), normalize(R)), 0.0f), sphere.shininess);

            vec3 specular = vec3(1, 1, 1);
            specular = specularK * p * light.intensity;
            vec3 diffuse =
                (float)std::max((double)dot(normalize(-light.direction), normal), 0.0) * sphere.color * light.intensity *
                diffuseK;
            sum_lights += (diffuse + specular);
        }

        for (int i = 0; i < spotLights.size(); i++) {
            SpotLight light = spotLights[i];

            vec3 spotLightDir = coordinate - light.position; // spotlight dir
            Ray ray = Ray(coordinate, -1.0f * spotLightDir);
            std::pair<int, int> sphereHit = minHit(ray, 1, hit.first);
            if (sphereHit.second == 1) {
                sum_lights += vec3(0,0,0);
                continue;
            }

            vec3 spot_light_dir = normalize(light.position - coordinate);
            vec3 L = normalize(light.direction);
            float spotCosine = -1.0f * dot(L, spot_light_dir);

            if (spotCosine >= light.angle) {
                vec3 normal = normalize(coordinate - sphere.center);
                float diff = (dot(normal, spot_light_dir));
                vec3 diffuse = light.intensity * diff * diffuseK * sphere.color;
                sum_lights += diffuse;

                vec3 R = -1.0f * spot_light_dir - 2 * dot(-1.0f * spot_light_dir, normal) * normal;
                float p = pow(max(dot(normalize(-1.0f * r.direction()), normalize(R)), 0.0f), sphere.shininess);

                vec3 specular = vec3(1, 1, 1);
                specular = specularK * p * light.intensity;
                sum_lights += specular;
            }
        }
        vec3 color = sum_lights;
        if (sphere.matType == Reflective) {
            Ray reflect_ray = Ray(r.at(t_sphere),
                normalize(r.direction()) - 2 * dot(normalize(r.direction()), normal) * normal);
            color += specularK * ray_color(reflect_ray, depth - 1, 1, hit.first);
        }

        if (sphere.matType == Transparent) {
            vec3 refract_ray_in = refract(normalize(r.direction()), normal, (1.0f / 1.5f));
            Ray ray = Ray(r.at(t_sphere), refract_ray_in);
            float outPoint = hit_sphere(sphere.center, sphere.radius, ray);
            vec3 newNormal = normalize(sphere.center - ray.at(outPoint));
            vec3 refract_ray_out = refract(normalize(ray.direction()), newNormal, (1.5f / 1.0f));
            Ray refract_ray;
            refract_ray = Ray(ray.at(outPoint), refract_ray_out);
            color += 1.0f * ray_color(refract_ray, depth - 1, 1, hit.first);
        }

        if (color.x > 1.0f) {
            color = vec3(1.0, color.y, color.z);
        }
        if (color.y > 1.0f) {
            color = vec3(color.x, 1.0, color.z);
        }
        if (color.z > 1.0f) {
            color = vec3(color.x, color.y, 1.0f);
        }


        return color;


    }
    else { //plane
        //        return vec3(0,0,0);
        return ray_color_plane(r, hit.first, depth);
    }
}

unsigned char* part1() {

    float mat[256][256][3];
    vec3 origin = vec3(0.0, 0.0, 3.0);
    ambient = vec3(0.2, 0.2, 0.3);
    samples_per_pixel = 15;

    // Scene from PDF
//    spheres.push_back(Sphere(vec3(-0.7, -0.7, -2.0), 0.5, vec3(1, 0, 0), 10, Normal));
//    spheres.push_back(Sphere(vec3(0.6, -0.5, -1.0), 0.5, vec3(0.6, 0.0, 0.8), 10, Normal));
//    planes.push_back(Plane(vec4(0, -0.5, -1.0, -3.5), vec3(0.0, 1.0, 1.0), 10, Normal));
//    dirLights.push_back(DirLight(vec3(0.0, 0.5, -1.0), vec3(0.7, 0.5, 0.0)));
//    spotLights.push_back(SpotLight(vec3(2.0, 1.0, 3.0), vec3(0.5, 0.0, -1), vec3(0.2, 0.5, 0.7)));


    // Scene 3

//    planes.push_back(Plane(vec4(1, 0, -0.1, -3.0), vec3(0.6, 0.0, 0.8), 20, Normal));
//    planes.push_back(Plane(vec4(0, 0, -1.0, -3.5), vec3(0.7, 0.7, 0), 10, Normal));
//    planes.push_back(Plane(vec4(-1, 0, -0.1, -3.0), vec3(0.0, 0.9, 0.5), 15, Normal));
//    planes.push_back(Plane(vec4(0, 1, -0.1, -3.0), vec3(0.0, 0.8, 0.8), 10, Normal));
//    planes.push_back(Plane(vec4(0, -1, -0.1, -3.0), vec3(0.9, 0.0, 0.1), 15, Normal));
//    spotLights.push_back(SpotLight(vec3(0.0, 0.0, 0.0), vec3(0.0, 0.5, -1), vec3(0.3, 0.9, 0.2), 0.8));
//    spotLights.push_back(SpotLight(vec3(0.0, 0.0, 0.0), vec3(0.5, 0.0, -1), vec3(0.9, 0.5, 0.5), 0.9));
//    spotLights.push_back(SpotLight(vec3(-0.2, 0.0, 0.0), vec3(-0.4, -0.3, -1), vec3(0.3, 0.5, 0.9), 0.7));

    // Scene 4
    spheres.push_back(Sphere(vec3(-0.0, -0.5, -1.0), 0.3, vec3(0, 1, 1), 10, Normal));
    spheres.push_back(Sphere(vec3(-1.7, -0.7, -3.0), 0.7, vec3(1, 0.0, 0), 20, Normal));
    spheres.push_back(Sphere(vec3(2.6, -1.5, -10.0), 2.5, vec3(0.6, 0, 0.8), 15, Normal));
    spheres.push_back(Sphere(vec3(1.3, 1.5, -7.0), 1.5, vec3(0.9, 0.0, 0.0), 10, Normal));
    spheres.push_back(Sphere(vec3(-0.6, -0.5, -5.0), 1.0, vec3(0.0, 0.0, 0.8), 10, Normal));
    planes.push_back(Plane(vec4(0, -1, -1, -8.5), vec3(0.7, 0.7, 0.0), 10, Normal));
    dirLights.push_back(DirLight(vec3(0.0, -0.7, -1.0), vec3(0.9, 0.5, 0.0)));
    spotLights.push_back(SpotLight(vec3(-2.0, -1.0, 3.0), vec3(1.5, 0.9, -1), vec3(0.2, 0.5, 0.7), 0.6));


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
                vec3 color = ray_color(r, 2, -1, -1);
                pixelColor += color;
            }
            //            vec3 color = ray_color_plane(r);
            mat[y][x][0] = pixelColor[0] / samples_per_pixel;
            mat[y][x][1] = pixelColor[1] / samples_per_pixel;
            mat[y][x][2] = pixelColor[2] / samples_per_pixel;

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

    //    part1();

    //	AddTexture("../res/textures/box0.bmp",false);

    unsigned char* zibi = part1();
    AddTexture(256, 256, zibi);

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
