#include "game.h"
#include "Ray.h"
#include <iostream>
#include <algorithm>
#include <glm/gtc/matrix_transform.hpp>

using namespace glm;


enum MaterialType
{
    Normal, Reflective, Transparent
};


// TODO combine all vectors to Vector<Sphere>
std::vector<glm::vec4> spheres;
std::vector<glm::vec3> spheres_colors;
std::vector<MaterialType> spheres_materials;

std::vector<glm::vec4> planes;
std::vector<glm::vec3> planes_colors;
std::vector<MaterialType> planes_materials;

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

vec3 ray_color(Ray &r, int depth , int planeOrSphere, int objectIndex) ;// forward deceleration

vec3 refract(vec3 ray, vec3 normal, float refract_ratio) {
    float cos_theta = min(dot(-ray, normal), 1.0f);
    vec3 r_out_perp = refract_ratio * (ray + cos_theta * normal);
    vec3 r_out_parallel = normal * (float) (-sqrt(abs(1.0 - dot(r_out_perp, r_out_perp))));
    return r_out_perp + r_out_parallel;
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
                i++;
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
                i++;
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

vec3 ray_color_plane(Ray &r, int index, int depth) {
    vec4 plane = planes[index];
    vec3 plane_color = planes_colors[index];
    float diffuseK = 1;
    float specularK = 0.7;
    float spot_factor = 1.0;
    float t = hit_plane(plane[0], plane[1], plane[2], plane[3], r);

    if (t > 0) {
        vec3 normal = normalize(vec3(plane[0], plane[1], plane[2]));
        vec3 sum_lights = ambient * plane_color;
        vec3 diffuse_and_specular = vec3(0, 0, 0);
        vec3 coordinate = r.at(t);
        float chess = floor(coordinate.x * 1.5) + floor(coordinate.y * 1.5);
        chess = 2 * ((chess * 0.5) - int(chess * 0.5));
        diffuseK = chess == 0 ? 1 : 0.5;
        for (int i = 0; i < directional_dirs.size(); i++) {
//            continue;
            vec3 light_dir = directional_dirs[i];
            Ray ray = Ray(coordinate, -1.0f * light_dir);
            std::pair<int, int> hit = minHit(ray, 2, index);
            if (hit.second == 1) {
                sum_lights += vec3(0, 0, 0);
                continue;
            }
            vec3 light_color = directional_intensity[i];
            vec3 normal = normalize(vec3(plane[0], plane[1], plane[2]));
            vec3 R = light_dir - 2 * dot(light_dir, normal) * normal;
            float p = pow(max(dot(normalize(-1.0f * r.direction()), normalize(R)), 0.0f), 256);
            vec3 specular = vec3(0, 0, 0);
            specular = specularK * p * light_color;
            vec3 diffuse = (float) std::max((double) dot(normalize(light_dir), normal), 0.0)  * plane_color * light_color * diffuseK;
            diffuse_and_specular += diffuse + specular;
            sum_lights += diffuse_and_specular;
        }
        for (int i = 0; i < spotLights_dirs.size(); i++) {
            vec3 spot_light_dir = spotLights_positions[i] - coordinate ; // spotlight dir
            Ray ray = Ray(coordinate,  spot_light_dir);
            std::pair<int, int> hit = minHit(ray, 2, index);
            if (hit.second == 1) {
                sum_lights += vec3(0, 0, 0);
                continue;
            }
            vec3 spot_light_color = spotLights_intensity[i]; // spotlight
            vec3 spot_light_origin = spotLights_positions[i]; // spotlight
            vec3 L = normalize(spotLights_dirs[i]);

            vec3 D = normalize(r.at(t) - spot_light_origin);
            float spotCosine = dot(D, L);

            if (spotCosine  >= 0.6) { // ?
                vec3 light_color = spotLights_intensity[i];
                vec3 normal = normalize(vec3(plane[0], plane[1], plane[2]));
                vec3 R = 1.0f * spot_light_dir - 2 * dot(-1.0f * spot_light_dir, normal) * normal;//?
                float p = pow(max(dot(normalize(1.0f * r.direction()), normalize(R)), 0.0f), 8);
                vec3 specular = vec3(0, 0, 0);
                specular = specularK * p * light_color;
                vec3 diffuse =spot_light_color * plane_color * (float) std::max((double) dot(normalize(1.0f*spot_light_dir), normal), 0.0) *
                               diffuseK;
                sum_lights += specular + diffuse;
            }
        }
        if (planes_materials[index] == Reflective)
        {
            Ray reflect_ray = Ray(r.at(t), normalize(r.direction()) - 2 * dot(normalize(r.direction()), -normal) * -normal);
            sum_lights += specularK * ray_color(reflect_ray, depth-1, 1, index);
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

vec3 ray_color(Ray &r, int depth , int planeOrSphere, int objectIndex) {

    // maximum recursion
    if (depth == 0) {
        return vec3(0,0,0);
    }

    float diffuseK = 1;
    float specularK = 0.7;
    std::pair<int, int> hit = minHit(r, planeOrSphere, objectIndex);

    if (hit.first == -1) { //background
        return vec3(0,0,0);
        vec3 unit_dir = normalize(r.direction());
        float t = 0.5 * unit_dir[1] + 1.0;
        float t2 = 1.0 - t;
        return t2 * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
    }
    if (hit.second == 1) { //sphere
        vec3 sphere_center = vec3(spheres[hit.first].x, spheres[hit.first].y, spheres[hit.first].z);
        float t_sphere = hit_sphere(sphere_center, spheres[hit.first].w, r);
        vec3 normal = normalize(r.at(t_sphere) - sphere_center);

        vec3 sum_lights = ambient * spheres_colors[hit.first];
        for (int i = 0; i < directional_dirs.size(); i++) {
            vec3 light_dir = directional_dirs[i];
            vec3 light_color = directional_intensity[i];
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
                vec3 normal = normalize(r.at(t_sphere) - sphere_center);
                float diff = (dot(normal, spot_light_dir));
                vec3 diffuse = light_color * diff * diffuseK * spheres_colors[hit.first];
                sum_lights += diffuse;

                vec3 R = -1.0f * spot_light_dir - 2 * dot(-1.0f * spot_light_dir, normal) * normal;
                float p = pow(max(dot(normalize(-1.0f * r.direction()), normalize(R)), 0.0f), 8);

                vec3 specular = vec3(1, 1, 1);
                specular = specularK * p * light_color;
                sum_lights += specular;
            }
        }
        vec3 color = sum_lights;
        if (spheres_materials[hit.first] == Reflective)
        {
            Ray reflect_ray = Ray(r.at(t_sphere), normalize(r.direction()) - 2 * dot(normalize(r.direction()), normal) * normal);
            color += specularK * ray_color(reflect_ray, depth-1, 1, hit.first);
        }

        if (spheres_materials[hit.first] == Transparent)
        {
            vec3 refract_ray_in = refract(normalize(r.direction()), normal, (1.0f/1.5f));
            Ray ray = Ray(r.at(t_sphere), refract_ray_in);
            float outPoint = hit_sphere(sphere_center, spheres[hit.first].w, ray);
            vec3 newNormal = normalize(sphere_center - ray.at(outPoint));
            vec3 refract_ray_out = refract(normalize(ray.direction()), newNormal, (1.5f/1.0f));
            Ray refract_ray = Ray(ray.at(outPoint), refract_ray_out);
//            Ray refract_ray = Ray(r.at(t_sphere + 0.2), r.direction());
            color = 1.0f * ray_color(refract_ray, depth-1, 1, hit.first);
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


    } else { //plane
//        return vec3(0,0,0);
        return ray_color_plane(r, hit.first, depth);
    }
}

unsigned char *part1() {
    // PDF Test scene
//    vec3 center = vec3(-0.7, -0.7, -2.0); // sphere center
//    vec3 center2 = vec3(0.6, -0.5, -1.0); // sphere2 center

    float mat[256][256][3];
    vec3 origin = vec3(0.0, 0.0, 4.0);
    ambient = vec3(0.1, 0.2, 0.3);

    vec3 sphere_color = vec3(0, 1, 0);
    vec3 center = vec3(-0.4, 1, -1.0); // sphere center
    spheres.push_back(vec4(center.x, center.y, center.z, 0.3));
    spheres_colors.push_back(sphere_color);
    spheres_materials.push_back(Transparent);

    vec3 sphere_color_2 = vec3(0.6, 0.0, 0.8);
    vec3 center2 = vec3(0.7, 0, -1.0); // sphere2 center
    spheres.push_back(vec4(center2.x, center2.y, center2.z, 0.5));
    spheres_colors.push_back(sphere_color_2);
    spheres_materials.push_back(Transparent);

    vec3 sphere_color_3 = vec3(0.2, 0.3, 0.4);
    vec3 center3 = vec3(-0.4, 0, -1.0); // sphere2 center
    spheres.push_back(vec4(center3.x, center3.y, center3.z, 0.5));
    spheres_colors.push_back(sphere_color_3);
    spheres_materials.push_back(Transparent);

    planes.push_back(vec4(0, -0.5, -1.0, -3.5));
    vec3 plane_color = vec3(1, 0, 0);
    planes_colors.push_back(plane_color);
    planes_materials.push_back(Normal);

    vec3 spot_light_dir = normalize(vec3(0.5, 0, -1)); // spotlight dir
    vec3 spot_light_color = vec3(0.2, 0.5, 0.7); // spotlight
    vec3 spot_light_origin = vec3(2.0, 1.0, 3.0); // spotlight
    vec3 light_dir = normalize(vec3(0.5, -0.5, -1)); // light dir
//    vec3 light_dir = normalize(normalize(1.0f * vec3(0, -0.5, -1.0))); // light dir
    vec3 light_color = vec3(0.7, 0.5, 0.0);
//    vec3 light_color = vec3(1, 1, 1);

//    spotLights_dirs.push_back(spot_light_dir);
//    spotLights_intensity.push_back(spot_light_color);
//    spotLights_positions.push_back(spot_light_origin);
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
            vec3 color = ray_color(r, 2, -1, -1);
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
