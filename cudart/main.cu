#include <iostream>
#include <time.h>
#include <float.h>
#include <curand_kernel.h>
#include "vec3.h"
#include "ray.h"
#include "sphere.h"
#include "hitable_list.h"
#include "camera.h"
#include "material.h"

// limited version of checkCudaErrors from helper_cuda.h in CUDA examples
#define checkCudaErrors(val) check_cuda( (val), #val, __FILE__, __LINE__ )
#define NUM_SPHERES 488
#define NUM_LAM_SPHERES 2
#define NUM_MET_SPHERES 1
#define NUM_DIE_SPHERES 485

void check_cuda(cudaError_t result, char const *const func, const char *const file, int const line) {
    if (result) {
        std::cerr << "CUDA error = " << static_cast<unsigned int>(result) << " at " <<
            file << ":" << line << " '" << func << "' \n";
        // Make sure we call CUDA Device Reset before exiting
        cudaDeviceReset();
        exit(99);
    }
}

__device__ bool sphere_hit(const sphere &s, const ray& r, float t_min, float t_max, hit_record& rec) {
    vec3 oc = r.origin() - s.center;
    float a = dot(r.direction(), r.direction());
    float b = dot(oc, r.direction());
    float c = dot(oc, oc) - s.radius*s.radius;
    float discriminant = b*b - a*c;
    if (discriminant > 0) {
        float temp = (-b - sqrt(discriminant))/a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - s.center) / s.radius;
            rec.mat_ptr = s.mat_ptr;
            return true;
        }
        temp = (-b + sqrt(discriminant)) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - s.center) / s.radius;
            rec.mat_ptr = s.mat_ptr;
            return true;
        }
    }
    return false;
}

__device__ bool world_hit(const ray& r, const sphere** list, int list_size, float t_min, float t_max, hit_record& rec) {
        hit_record temp_rec;
        bool hit_anything = false;
        float closest_so_far = t_max;
        for (int i = 0; i < list_size; i++) {
            if (sphere_hit((*list)[i], r, t_min, closest_so_far, temp_rec)) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
                if (i < NUM_LAM_SPHERES) {
                    rec.type = 0;
                }
                else if (i >= NUM_LAM_SPHERES && i < NUM_LAM_SPHERES + NUM_MET_SPHERES) {
                    rec.type = 1;
                }
                else if (i >= NUM_LAM_SPHERES + NUM_MET_SPHERES && i < NUM_LAM_SPHERES + NUM_MET_SPHERES + NUM_DIE_SPHERES) {
                    rec.type = 2;
                }
            }
        }
        return hit_anything;
}

__device__ bool lambertian_scatter(const lambertian* l, const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered, curandState *local_rand_state) {
    vec3 target = rec.p + rec.normal + random_in_unit_sphere(local_rand_state);
    scattered = ray(rec.p, target-rec.p);
    attenuation = l->albedo;
    return true;
}

__device__ bool metal_scatter(const metal* m, const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered, curandState *local_rand_state) {
    vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
    scattered = ray(rec.p, reflected + m->fuzz*random_in_unit_sphere(local_rand_state));
    attenuation = m->albedo;
    return (dot(scattered.direction(), rec.normal) > 0.0f);
}

__device__ bool dielectric_scatter(
                        const dielectric* d,
                        const ray& r_in,
                        const hit_record& rec,
                        vec3& attenuation,
                        ray& scattered,
                        curandState *local_rand_state) {
    vec3 outward_normal;
    vec3 reflected = reflect(r_in.direction(), rec.normal);
    float ni_over_nt;
    attenuation = vec3(1.0, 1.0, 1.0);
    vec3 refracted;
    float reflect_prob;
    float cosine;
    if (dot(r_in.direction(), rec.normal) > 0.0f) {
        outward_normal = -rec.normal;
        ni_over_nt = d->ref_idx;
        cosine = dot(r_in.direction(), rec.normal) / r_in.direction().length();
        cosine = sqrt(1.0f - d->ref_idx*d->ref_idx*(1-cosine*cosine));
    }
    else {
        outward_normal = rec.normal;
        ni_over_nt = 1.0f / d->ref_idx;
        cosine = -dot(r_in.direction(), rec.normal) / r_in.direction().length();
    }
    if (refract(r_in.direction(), outward_normal, ni_over_nt, refracted))
        reflect_prob = schlick(cosine, d->ref_idx);
    else
        reflect_prob = 1.0f;
    if (curand_uniform(local_rand_state) < reflect_prob)
        scattered = ray(rec.p, reflected);
    else
        scattered = ray(rec.p, refracted);
    return true;
}

// Matching the C++ code would recurse enough into color() calls that
// it was blowing up the stack, so we have to turn this into a
// limited-depth loop instead.  Later code in the book limits to a max
// depth of 50, so we adapt this a few chapters early on the GPU.
__device__ vec3 color(const ray& r, const sphere** list, int list_size, curandState *local_rand_state) {
    ray cur_ray = r;
    vec3 cur_attenuation = vec3(1.0,1.0,1.0);
    for(int i = 0; i < 50; i++) {
        hit_record rec;
        if (world_hit(cur_ray, list, list_size, 0.001f, FLT_MAX, rec)) {
            ray scattered;
            vec3 attenuation;
            if (rec.type == 0 && lambertian_scatter((lambertian*)rec.mat_ptr, cur_ray, rec, attenuation, scattered, local_rand_state)) {
                cur_attenuation *= attenuation;
                cur_ray = scattered;
            }
            else if (rec.type == 1 && metal_scatter((metal*)rec.mat_ptr, cur_ray, rec, attenuation, scattered, local_rand_state)) {
                cur_attenuation *= attenuation;
                cur_ray = scattered;
            }
            else if (rec.type == 2 && dielectric_scatter((dielectric*)rec.mat_ptr, cur_ray, rec, attenuation, scattered, local_rand_state)) {
                cur_attenuation *= attenuation;
                cur_ray = scattered;
            }
            else {
                return vec3(0.0,0.0,0.0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(cur_ray.direction());
            float t = 0.5f*(unit_direction.y() + 1.0f);
            vec3 c = (1.0f-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
            return cur_attenuation * c;
        }
    }
    return vec3(0.0,0.0,0.0); // exceeded recursion
}

__global__ void rand_init(curandState *rand_state) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        curand_init(1984, 0, 0, rand_state);
    }
}

__global__ void render_init(int max_x, int max_y, curandState *rand_state) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    if((i >= max_x) || (j >= max_y)) return;
    int pixel_index = j*max_x + i;
    // Original: Each thread gets same seed, a different sequence number, no offset
    // curand_init(1984, pixel_index, 0, &rand_state[pixel_index]);
    // BUGFIX, see Issue#2: Each thread gets different seed, same sequence for
    // performance improvement of about 2x!
    curand_init(1984+pixel_index, 0, 0, &rand_state[pixel_index]);
}

__global__ void render(vec3 *fb, int max_x, int max_y, int ns, camera **cam, hitable **world, curandState *rand_state) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    if((i >= max_x) || (j >= max_y)) return;

    hitable_list **hitable_list_word = (hitable_list**)(world);

    extern __shared__ char array[];
    int index = 0;
    sphere *sphere_list = (sphere *)array;
    index += NUM_SPHERES*sizeof(sphere);
    lambertian* lam_list = (lambertian *)(array + index);
    index += NUM_LAM_SPHERES*sizeof(lambertian);
    metal* met_list = (metal *)(array + index);
    index += NUM_MET_SPHERES*sizeof(metal);
    dielectric* die_list = (dielectric *)(array + index);
    index += NUM_DIE_SPHERES*sizeof(dielectric);

    for (int i = 0; i < NUM_LAM_SPHERES; ++i) {
        sphere* s = (sphere*)(*((*hitable_list_word)->list + i));
        lam_list[i].albedo.e[0] = ((lambertian*)(s->mat_ptr))->albedo.e[0];
        lam_list[i].albedo.e[1] = ((lambertian*)(s->mat_ptr))->albedo.e[1];
        lam_list[i].albedo.e[2] = ((lambertian*)(s->mat_ptr))->albedo.e[2];
    }

    for (int i = NUM_LAM_SPHERES; i < NUM_LAM_SPHERES+NUM_MET_SPHERES; ++i) {
        int met_i = i - NUM_LAM_SPHERES;
        sphere* s = (sphere*)(*((*hitable_list_word)->list + i));
        met_list[met_i].albedo.e[0] = ((metal*)(s->mat_ptr))->albedo.e[0];
        met_list[met_i].albedo.e[1] = ((metal*)(s->mat_ptr))->albedo.e[1];
        met_list[met_i].albedo.e[2] = ((metal*)(s->mat_ptr))->albedo.e[2];
        met_list[met_i].fuzz = ((metal*)(s->mat_ptr))->fuzz;
    }

    for (int i = NUM_LAM_SPHERES+NUM_MET_SPHERES; i < NUM_LAM_SPHERES+NUM_MET_SPHERES+NUM_DIE_SPHERES; ++i) {
        int die_i = i - NUM_LAM_SPHERES-NUM_MET_SPHERES;
        sphere* s = (sphere*)(*((*hitable_list_word)->list + i));
        die_list[die_i].ref_idx = ((dielectric*)(s->mat_ptr))->ref_idx;
    }

    for (int i = 0; i < (*hitable_list_word)->list_size; ++i) {
        sphere* s = (sphere*)(*((*hitable_list_word)->list + i));
        sphere_list[i].center.e[0] = s->center.e[0];
        sphere_list[i].center.e[1] = s->center.e[1];
        sphere_list[i].center.e[2] = s->center.e[2];
        sphere_list[i].radius = s->radius;
        if (i < NUM_LAM_SPHERES) {
            sphere_list[i].mat_ptr = &lam_list[i];
        }
        else if (i >= NUM_LAM_SPHERES && i < NUM_LAM_SPHERES + NUM_MET_SPHERES) {
            sphere_list[i].mat_ptr = &met_list[i-NUM_LAM_SPHERES];
        }
        else if (i >= NUM_LAM_SPHERES + NUM_MET_SPHERES && i < NUM_LAM_SPHERES + NUM_MET_SPHERES + NUM_DIE_SPHERES) {
            sphere_list[i].mat_ptr = &die_list[i-NUM_LAM_SPHERES-NUM_MET_SPHERES];
        }
    }

    int pixel_index = j*max_x + i;
    curandState local_rand_state = rand_state[pixel_index];
    vec3 col(0,0,0);
    for(int s=0; s < ns; s++) {
        float u = float(i + curand_uniform(&local_rand_state)) / float(max_x);
        float v = float(j + curand_uniform(&local_rand_state)) / float(max_y);
        ray r = (*cam)->get_ray(u, v, &local_rand_state);
        col += color(r, &sphere_list, (*hitable_list_word)->list_size, &local_rand_state);
    }
    rand_state[pixel_index] = local_rand_state;
    col /= float(ns);
    col[0] = sqrt(col[0]);
    col[1] = sqrt(col[1]);
    col[2] = sqrt(col[2]);
    fb[pixel_index] = col;
}

#define RND (curand_uniform(&local_rand_state))

__global__ void create_world(hitable **d_list, hitable **d_world, camera **d_camera, int nx, int ny, curandState *rand_state) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        curandState local_rand_state = *rand_state;
        d_list[0] = new sphere(vec3(0,-1000.0,-1), 1000,
                               new lambertian(vec3(0.5, 0.5, 0.5)));
        int i = 1;
        for (; i < NUM_LAM_SPHERES-1; ++i) {
            int a = i / 22 - 11;
            int b = i % 22 - 11;
            vec3 center(a+RND,0.2,b+RND);
            d_list[i] = new sphere(center, 0.2,
                                     new lambertian(vec3(RND*RND, RND*RND, RND*RND)));
        }
        d_list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));
        for (; i < NUM_LAM_SPHERES+NUM_MET_SPHERES-1; ++i) {
            int a = i / 22 - 11;
            int b = i % 22 - 11;
            vec3 center(a+RND,0.2,b+RND);
            d_list[i] = new sphere(center, 0.2,
                                        new metal(vec3(0.5f*(1.0f+RND), 0.5f*(1.0f+RND), 0.5f*(1.0f+RND)), 0.5f*RND));
        }
        d_list[i++] = new sphere(vec3(4, 1, 0),  1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));
        for (; i < NUM_LAM_SPHERES+NUM_MET_SPHERES+NUM_DIE_SPHERES-1; ++i) {
            int a = i / 22 - 11;
            int b = i % 22 - 11;
            vec3 center(a+RND,0.2,b+RND);
            d_list[i] = new sphere(center, 0.2, new dielectric(1.5));
        }
        d_list[i++] = new sphere(vec3(0, 1,0),  1.0, new dielectric(1.5));
        *rand_state = local_rand_state;
        *d_world  = new hitable_list(d_list, NUM_SPHERES);

        vec3 lookfrom(13,2,3);
        vec3 lookat(0,0,0);
        float dist_to_focus = 10.0; (lookfrom-lookat).length();
        float aperture = 0.1;
        *d_camera   = new camera(lookfrom,
                                 lookat,
                                 vec3(0,1,0),
                                 30.0,
                                 float(nx)/float(ny),
                                 aperture,
                                 dist_to_focus);
    }
}

__global__ void free_world(hitable **d_list, hitable **d_world, camera **d_camera) {
    for(int i=0; i < NUM_SPHERES; i++) {
        delete ((sphere *)d_list[i])->mat_ptr;
        delete d_list[i];
    }
    delete *d_world;
    delete *d_camera;
}

int main() {
    int nx = 1200;
    int ny = 800;
    int ns = 10;
    int tx = 8;
    int ty = 8;

    std::cerr << "Rendering a " << nx << "x" << ny << " image with " << ns << " samples per pixel ";
    std::cerr << "in " << tx << "x" << ty << " blocks.\n";

    int num_pixels = nx*ny;
    size_t fb_size = num_pixels*sizeof(vec3);

    // allocate FB
    vec3 *fb;
    checkCudaErrors(cudaMallocManaged((void **)&fb, fb_size));

    // allocate random state
    curandState *d_rand_state;
    checkCudaErrors(cudaMalloc((void **)&d_rand_state, num_pixels*sizeof(curandState)));
    curandState *d_rand_state2;
    checkCudaErrors(cudaMalloc((void **)&d_rand_state2, 1*sizeof(curandState)));

    // we need that 2nd random state to be initialized for the world creation
    rand_init<<<1,1>>>(d_rand_state2);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());

    // make our world of hitables & the camera
    hitable **d_list;
    int num_hitables = NUM_SPHERES;
    checkCudaErrors(cudaMalloc((void **)&d_list, num_hitables*sizeof(hitable *)));
    hitable **d_world;
    checkCudaErrors(cudaMalloc((void **)&d_world, sizeof(hitable *)));
    camera **d_camera;
    checkCudaErrors(cudaMalloc((void **)&d_camera, sizeof(camera *)));
    create_world<<<1,1>>>(d_list, d_world, d_camera, nx, ny, d_rand_state2);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());

    clock_t start, stop;
    start = clock();
    // Render our buffer
    dim3 blocks(nx/tx+1,ny/ty+1);
    dim3 threads(tx,ty);
    render_init<<<blocks, threads>>>(nx, ny, d_rand_state);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());
    render<<<blocks, threads, NUM_SPHERES*sizeof(sphere)+NUM_LAM_SPHERES*sizeof(lambertian)+NUM_MET_SPHERES*sizeof(metal)+NUM_DIE_SPHERES*sizeof(dielectric)>>>(fb, nx, ny,  ns, d_camera, d_world, d_rand_state);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());
    stop = clock();
    double timer_seconds = ((double)(stop - start)) / CLOCKS_PER_SEC;
    std::cerr << "took " << timer_seconds << " seconds.\n";

    // Output FB as Image
    std::cout << "P3\n" << nx << " " << ny << "\n255\n";
    for (int j = ny-1; j >= 0; j--) {
        for (int i = 0; i < nx; i++) {
            size_t pixel_index = j*nx + i;
            int ir = int(255.99*fb[pixel_index].r());
            int ig = int(255.99*fb[pixel_index].g());
            int ib = int(255.99*fb[pixel_index].b());
            std::cout << ir << " " << ig << " " << ib << "\n";
        }
    }

    // clean up
    checkCudaErrors(cudaDeviceSynchronize());
    free_world<<<1,1>>>(d_list,d_world,d_camera);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaFree(d_camera));
    checkCudaErrors(cudaFree(d_world));
    checkCudaErrors(cudaFree(d_list));
    checkCudaErrors(cudaFree(d_rand_state));
    checkCudaErrors(cudaFree(d_rand_state2));
    checkCudaErrors(cudaFree(fb));

    cudaDeviceReset();
}
