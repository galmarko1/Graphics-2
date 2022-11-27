//
// Created by gal on 26/11/2022.
//

#ifndef GAME_RAY_H
#define GAME_RAY_H


#include "glm/vec3.hpp"

class Ray {
    public:
        Ray(){}
        Ray(glm::vec3 origin, glm::vec3 direction) : orig(origin), dir(direction)
        {}

        glm::vec3 origin() {return orig;}
        glm::vec3 direction() {return dir;}
        glm::vec3 at (float t) {
            glm::vec3 vec_t = glm::vec3(t,t,t);
            return orig + vec_t * dir;
        }
    public:
    glm::vec3 orig;
    glm::vec3 dir;

};


#endif //GAME_RAY_H
