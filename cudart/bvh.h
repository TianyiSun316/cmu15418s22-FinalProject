#ifndef BVH_H
#define BVH_H

#include "hitable.h"
#include "hitable_list.h"
#include "bbox.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>

using std::shared_ptr;
using std::make_shared;
using std::sqrt;

class bvh_node : public hitable  {
    public:
        __device__ bvh_node();
        __device__ bvh_node(hitable** src_objects, int start, int end, int axis);
        __device__ virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const override;
        __device__ virtual bool bounding_box(aabb& output_box) const override;

    public:
        hitable* left;
        hitable* right;
        aabb box;
        curandState local_rand_state;
};


__device__ bool box_compare(hitable* a, hitable* b, int axis) {
    aabb box_a;
    aabb box_b;
    a->bounding_box(box_a);
    b->bounding_box(box_b);

    bool ret = box_a.minimum.e[axis] < box_b.minimum.e[axis];
    return ret;
}


__device__ bool box_x_compare (hitable* a, hitable* b) {
    return box_compare(a, b, 0);
}

__device__ bool box_y_compare (hitable* a, hitable* b) {
    return box_compare(a, b, 1);
}

__device__ bool box_z_compare (hitable* a, hitable* b) {
    return box_compare(a, b, 2);
}

__device__ void sortobj(hitable** arr, int start, int end, int axis) {
    int N = end - start - 1;
    hitable** newarr = arr + start;
    for (int i = 0; i < N-1; i++) {
        bool swapped = false;
        for (int j = 0; j < N-1-i; j++) {
            if (box_compare(newarr[j], newarr[j + 1], axis)) {
                hitable* temp = newarr[i];
                newarr[i] = newarr[j + 1];
                newarr[j + 1] = temp;
                swapped = true;
            }
        }
        if (!swapped) {
            return;
        }
    }
}

__device__ bvh_node::bvh_node(hitable** objects, int start, int end, int axis) {
    // int axis = random_int(0,2);
    auto comparator = (axis == 0) ? box_x_compare
                    : (axis == 1) ? box_y_compare
                                  : box_z_compare;
    // printf("Axis=%d\n", axis);
    int object_span = end - start;
    if (object_span == 1) {
        left = right = objects[start];
    } else if (object_span == 2) {
        if (comparator(objects[start], objects[start+1])) {
            left = objects[start];
            right = objects[start+1];
        } else {
            left = objects[start+1];
            right = objects[start];
        }
    } else {
        if (start == 0) {
            left = objects[0];
            right = new bvh_node(objects, 1, end, axis);
        } else {
            thrust::sort(objects + start, objects + end, comparator);
            // printf("%d, %d, axis=%d\n", start, end, axis);
            // sortobj(objects, start, end, axis);
            int mid = start + object_span / 2;
            // left = new bvh_node(objects, start, mid, axis == 0 ? 2 : 0);
            // right = new bvh_node(objects, mid, end, axis == 0 ? 2 : 0);
            left = new bvh_node(objects, start, mid, (axis+1) % 3);
            right = new bvh_node(objects, mid, end, (axis+1) % 3);
        }
    }
    aabb box_left, box_right;
    left->bounding_box(box_left);
    right->bounding_box(box_right);
    box = surrounding_box(box_left, box_right);
}


__device__ bool bvh_node::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
    if (!box.hit(r, t_min, t_max))
        return false;

    bool hit_left = left->hit(r, t_min, t_max, rec);
    bool hit_right = right->hit(r, t_min, hit_left ? rec.t : t_max, rec);

    return hit_left || hit_right;
}


__device__ bool bvh_node::bounding_box(aabb& output_box) const {
    output_box = box;
    return true;
}


#endif