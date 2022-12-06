#include "game.h"
#include "Ray.h"
#include <iostream>
#include <algorithm>
#include <glm/gtc/matrix_transform.hpp>

using namespace glm;

std::vector<glm::vec4> spheres;
std::vector<glm::vec3> spheres_colors;

std::vector<glm::vec4> planes;
std::vector<glm::vec3> planes_colors;

std::vector<glm::vec3> directional_dirs;
std::vector<glm::vec3> spotLights_dirs;
std::vector<glm::vec3> spotLights_positions;
std::vector<glm::vec3> spotLights_intensity;
std::vector<glm::vec3> directional_intensity;
glm::vec3 ambient;


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

float hit_sphere(vec3 center, double radius, Ray r) {
    vec3 oc = r.origin() - center;
    float a = dot(r.direction(), r.direction());
    float b = 2.0 * dot(oc, r.direction());
    float c = dot(oc, oc) - radius * radius;
    float discriminant = b * b - 4 * a * c;
    if (discriminant < 0) //no hit
        return -1;
    return (-b - sqrt(discriminant)) / (2.0 * a);
}


float hit_plane(float a, float b, float c, float d, Ray r) {
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

std::pair<int, int> minHit(Ray &r, int planeOrSphere, int objectIndex) {
    float min_t_sphere = 100000;
    float min_t_plane = 100000;
    int min_sphere_index = -1;
    int min_plane_index = -1;
    int i = 0;

    for (vec4 sphere: spheres) {
        if (planeOrSphere == 1) {
            if (objectIndex == i) {
                continue;
            }
        }
        float curr_t = hit_sphere(vec3(sphere.x, sphere.y, sphere.z), sphere.w, r);
        if (curr_t > -1) {
            if (min_t_sphere > curr_t) {
                min_sphere_index = i;
                min_t_sphere = curr_t;
            }
        }
        i++;
    }

    i = 0;
    for (vec4 plane: planes) {
        if (planeOrSphere == 2) {
            if (objectIndex == i) {
                continue;
            }
        }
        float curr_t = hit_plane(plane.x, plane.y, plane.z, plane.w, r);
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
    } else if (min_t_plane < min_t_sphere) {
        p.first = min_plane_index;
        p.second = 2; // plane vector
    } else {
        p.first = min_sphere_index;
        p.second = 1; // sphere vector
    }
    return p;
}

vec3 ray_color_plane(Ray &r, int index) {
    vec4 plane = planes[index];
    vec3 plane_color = planes_colors[index];
    float diffuseK = 1;
    float specularK = 0.7;
    float spot_factor = 1.0;
    float t = hit_plane(plane[0], plane[1], plane[2], plane[3], r);

    if (t > 0) {
        vec3 normal = normalize(vec3(plane[0], plane[1], plane[2]));
        vec3 sum_lights = ambient;
        vec3 diffuse_and_specular = vec3(0, 0, 0);
        vec3 coordinate = r.at(t);
        for (int i = 0; i < directional_dirs.size(); i++) {
            vec3 light_dir = directional_dirs[i];
            Ray ray = Ray(coordinate, -1.0f * light_dir);
            std::pair<int, int> hit = minHit(ray, 2, index);
            if (hit.second == 1) {
                sum_lights += vec3(0, 0, 0);
//                return vec3(0, 0, 0);
                continue;
            }
            vec3 light_color = directional_intensity[i];
            vec3 normal = normalize(vec3(plane[0], plane[1], plane[2]));
            vec3 R = light_dir - 2 * dot(light_dir, normal) * normal;
            float p = pow(max(dot(normalize(-1.0f * r.direction()), normalize(R)), 0.0f), 256);
            vec3 specular = vec3(0, 0, 0);
            specular = specularK * p * light_color;
            vec3 diffuse = (float) std::max((double) dot(normalize(light_dir), normal), 0.0) * plane_color * diffuseK;
            diffuse_and_specular += diffuse + specular;
            sum_lights += diffuse_and_specular;
        }
        for (int i = 0; i < spotLights_dirs.size(); i++) {
            vec3 spot_light_dir = spotLights_dirs[i]; // spotlight dir
            Ray ray = Ray(coordinate, -1.0f * spot_light_dir);
            std::pair<int, int> hit = minHit(ray, 2, index);
            if (hit.second == 1) {
                sum_lights += vec3(0, 0, 0);
//                return vec3(0,x 0, 0);
                continue;
            }
            vec3 spot_light_color = spotLights_intensity[i]; // spotlight
            vec3 spot_light_origin = spotLights_positions[i]; // spotlight
            vec3 L = normalize(spot_light_dir);
            vec3 D = normalize(r.at(t) - spot_light_origin);
            float spotCosine = dot(D, L);
            if (spotCosine >= 0.6) { // ?
                spot_factor = pow(spotCosine, 5);//?
            } else {
                spot_factor = 0.0;
            }
            vec3 spot = spot_factor * vec3(1, 1, 1);
            sum_lights += spot;
            vec3 light_color = directional_intensity[i];
            vec3 normal = normalize(vec3(plane[0], plane[1], plane[2]));
            vec3 R = spot_light_dir - 2 * dot(spot_light_dir, normal) * normal;
            float p = pow(max(dot(normalize(-1.0f * r.direction()), normalize(R)), 0.0f), 256);
            vec3 specular = vec3(0, 0, 0);
            specular = specularK * p * light_color;
//            vec3 diffuse = (float) std::max((double) dot(normalize(light_dir), normal), 0.0) * plane_color * diffuseK;
//            diffuse_and_specular += diffuse + specular;
            sum_lights += specular;
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

vec3 ray_color(Ray &r) {
    float diffuseK = 1;
    float specularK = 0.7;
    float spot_factor = 1.0;
    std::pair<int, int> hit = minHit(r, -1, -1);

    if (hit.first == -1) { //background
        vec3 unit_dir = normalize(r.direction());
        float t = 0.5 * unit_dir[1] + 1.0;
        float t2 = 1.0 - t;
        return t2 * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
    }
    if (hit.second == 1) { //sphere
        vec3 sphere_center = vec3(spheres[hit.first].x, spheres[hit.first].y, spheres[hit.first].z);
        float t_sphere = hit_sphere(sphere_center, spheres[hit.first].w, r);
        vec3 sum_lights = ambient;
        for (int i = 0; i < directional_dirs.size(); i++) {
            vec3 light_dir = directional_dirs[i];
            vec3 light_color = directional_intensity[i];

            vec3 normal = normalize(r.at(t_sphere) - sphere_center);
            vec3 R = 1.0f * light_dir - 2 * dot(light_dir, normal) * normal;
            float p = pow(max(dot(normalize(-1.0f * r.direction()), normalize(R)), 0.0f), 8);

            vec3 specular = vec3(1, 1, 1);
            specular = specularK * p * light_color;
            vec3 sphere_color = spheres_colors[hit.first];
            vec3 diffuse =
                    (float) std::max((double) dot(normalize(-light_dir), normal), 0.0) * sphere_color * light_color *
                    diffuseK;
            sum_lights += (diffuse + specular);
        }

        for (int i = 0; i < spotLights_dirs.size(); i++) {
            vec3 light_color = spotLights_intensity[i];

            vec3 spot_light_dir = normalize(spotLights_positions[i] - r.at(t_sphere));
            vec3 L = normalize(spotLights_dirs[i]);
            float spotCosine = -1.0f * dot(L, spot_light_dir);

            if (spotCosine >= 0.6) {
//                spot_factor = pow(spotCosine, 1);
//                vec3 spot = spot_factor * light_color;
                vec3 normal = normalize(r.at(t_sphere) - sphere_center);
                float diff = (dot(normal, spot_light_dir));
                vec3 product = vec3(light_color.x + spheres_colors[i].x, light_color.y + spheres_colors[i].y,
                                    light_color.z + spheres_colors[i].z);
                vec3 diffuse = product * diff * diffuseK;
                sum_lights += diffuse;

                vec3 R = -1.0f * spot_light_dir - 2 * dot(-1.0f * spot_light_dir, normal) * normal;
                float p = pow(max(dot(normalize(-1.0f * r.direction()), normalize(R)), 0.0f), 8);

                vec3 specular = vec3(1, 1, 1);
                specular = specularK * p * light_color;
                sum_lights += specular;
            }
        }
        vec3 color = sum_lights;

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


    } else { //plane
//        return vec3(0,0,0);
        return ray_color_plane(r, hit.first);
    }
}


unsigned char *part1() {
    vec3 center = vec3(-0.7, -0.7, -2.0); // sphere center
    vec3 center2 = vec3(0.6, -0.5, -1.0); // sphere2 center
    vec3 sphere_color = vec3(1, 0, 0);
    vec3 sphere_color_2 = vec3(0.6, 0.0, 0.8);
    vec3 plane_color = vec3(0, 1, 1);

    float mat[256][256][3];
    vec3 origin = vec3(0.0, 0.0, 4.0);
    ambient = vec3(0.1, 0.2, 0.3);
    spheres.push_back(vec4(center.x, center.y, center.z, 0.5));
    spheres.push_back(vec4(center2.x, center2.y, center2.z, 0.5));
    planes.push_back(vec4(0, -0.5, -1.0, -3.5));
    spheres_colors.push_back(sphere_color);
    planes_colors.push_back(plane_color);
    spheres_colors.push_back(sphere_color_2);

    vec3 spot_light_dir = normalize(vec3(0.5, 0, -1)); // spotlight dir
    vec3 spot_light_color = vec3(0.2, 0.5, 0.7); // spotlight
    vec3 spot_light_origin = vec3(2.0, 1.0, 3.0); // spotlight
    vec3 light_dir = normalize(vec3(0, 0.5, -1)); // light dir
    vec3 light_color = vec3(0.7, 0.5, 0.0);

    spotLights_dirs.push_back(spot_light_dir);
    spotLights_intensity.push_back(spot_light_color);
    spotLights_positions.push_back(spot_light_origin);
    directional_dirs.push_back(light_dir);
    directional_intensity.push_back(light_color);

    for (int y = 0; y < 256; y++) {
        for (int x = 0; x < 256; x++) {
            double u = double(x) / 255;
            double v = double(y) / 255;
            vec3 xDir = vec3(2 * u, 0, 0);
            vec3 yDir = vec3(0, 2 * v, 0);
            vec3 left_corner = vec3(-1, -1, 0);

            vec3 dir = left_corner + xDir + yDir - origin;
            Ray r = Ray(origin, dir);
            vec3 color = ray_color(r);
//            vec3 color = ray_color_plane(r);
            mat[y][x][0] = color[0];
            mat[y][x][1] = color[1];
            mat[y][x][2] = color[2];

        }
    }

    int index = 0;
    unsigned char *data = new unsigned char[256 * 256 * 4];

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

    unsigned char *zibi = part1();
    AddTexture(256, 256, zibi);

    AddShape(Plane, -1, TRIANGLES);

    pickedShape = 0;

    SetShapeTex(0, 0);
    MoveCamera(0, zTranslate, 10);
    pickedShape = -1;

    //ReadPixel(); //uncomment when you are reading from the z-buffer
}

void Game::Update(const glm::mat4 &MVP, const glm::mat4 &Model, const int shaderIndx) {
    Shader *s = shaders[shaderIndx];
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
