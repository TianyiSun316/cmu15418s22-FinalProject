#ifndef AABB_H
#define AABB_H

#include "ray.h"
#include "vec3.h"


class aabb {
    public:
        __device__ aabb() {}
        __device__ aabb(const vec3& a, const vec3& b) { minimum = a; maximum = b;}

        __device__ virtual bool hit(const ray& r, float t_min, float t_max) const;
        
        vec3 minimum;
        vec3 maximum;
};

__device__ bool aabb::hit(const ray& r, float t_min, float t_max) const {
    for (int a = 0; a < 3; a++) {
        float invD = 1.0f / r.direction()[a];
        float t0 = (minimum[a] - r.origin()[a]) * invD;
        float t1 = (maximum[a] - r.origin()[a]) * invD;
        if (invD < 0.0f) {
            float temp = t0;
            t0 = t1;
            t1 = temp;
        }
        t_min = t0 > t_min ? t0 : t_min;
        t_max = t1 < t_max ? t1 : t_max;
        if (t_max <= t_min)
            return false;
    }
    return true;
}

__device__ aabb surrounding_box(aabb box0, aabb box1) {
    vec3 small(fmin(box0.minimum.x(), box1.minimum.x()),
               fmin(box0.minimum.y(), box1.minimum.y()),
               fmin(box0.minimum.z(), box1.minimum.z()));

    vec3 big  (fmax(box0.maximum.x(), box1.maximum.x()),
               fmax(box0.maximum.y(), box1.maximum.y()),
               fmax(box0.maximum.z(), box1.maximum.z()));

    return aabb(small, big);
}

#endif