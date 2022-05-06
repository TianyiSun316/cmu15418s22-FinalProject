# Presentation Video
https://drive.google.com/file/d/1f5Ro3GWk_dcy1DaiIg2fLZNUFoMxV7mA/view?usp=sharing
# Final Report
https://github.com/TianyiSun316/cmu15418s22-FinalProject/blob/main/Parallel_Final_Project.pdf

# Project Proposal

## SUMMARY
We are going to implement ray tracing algorithm of both sequential and parallel version on GPU using CUDA. We will make use of GPU computing units to accelerate the process as quick as possible. We will also conduct detailed performance profiling for our implementation considering characteristics including memory bandwidth, shared memory usage and cache behavior. Besides, a comparision with the sequential implementation will also be made.

## BACKGROUND
The idea of ray tracing comes from as early as the 16th century when it was described by Albrecht Dürer, who is credited for its invention. In Four Books on Measurement, he described an apparatus called a Dürer's door using a thread attached to the end of a stylus that an assistant moves along the contours of the object to draw. The thread passes through the door's frame and then through a hook on the wall. The thread forms a ray and the hook acts as the center of projection and corresponds to the camera position in ray tracing. For decades, global illumination in major films using computer generated imagery was faked with additional lights. Since 2018, however, hardware acceleration for real-time ray tracing has become standard on new commercial graphics cards. Nowadays, ray tracing is becoming more and more accessible to game players and it's now the new trend of game image quality evolution.

## THE CHALLENGE
- The key challenge here is to determine a efficient way so that the computation is paralleled with balanced workload and lost cost. As there is a lot of memory access involoved in the algorithm, we also need to think very carefully about the way of leveraging the shared memory so that most of the memory access can be fast.
- There is a lot of physics of optics that we should implement such as reflection calculation, intersection calculation, etc. Besides, all these calculations should be implemented in a efficient way.

## RESOURCES
We plan to use ghc server with Nvidia GPU and CUDA environment. We will start from scratch for the ray tracing algorithm. We may still need to use OpenGL rasterizer for hybrid realtime rendering. Here are some reference papers that introduce ray tracing algorithm and implementation on GPU.

- [Purcell, Timothy J., et al. "Ray Tracing on Programmable Graphics Hardware."](https://graphics.stanford.edu/papers/rtongfx/rtongfx.pdf)
- [Glassner, Andrew S., ed. An introduction to ray tracing. Morgan Kaufmann, 1989.](https://www.realtimerendering.com/raytracing/An-Introduction-to-Ray-Tracing-The-Morgan-Kaufmann-Series-in-Computer-Graphics-.pdf)

## GOALS AND DELIVERABLES
### 75% goal
- Implement static ray tracing algorithm of both sequential version and GPU-parallel version. 
- Support light reflection and refraction on a single texture.
- Supports simple 3D geometry(sphere, cube, cone) rendering.
- Conduct performance analysis on our ray tracing algorithm compared with sequential version and on other type of GPUs.

### 100% goal
- Support 1080p, 30FPS real-time rendering.
- Support multiple object textures and multiple light sources.
- Conduct detailed performance analysis and find potential bottleneck and optimize our render.

### 125% goal
- Implement hybrid real-time rendering (ray tracing + rasterization).
- Optimize our render to support 4K, 100FPS real-time rendering.

## PLATFORM CHOICE
Language: C++

Platform: GHC machine with NVIDIA RTX 2080 GPU. (Our rendering system runs on single machine and normally requires no more than 8GB GPU memory.)

## SCHEDULE
- 4.1 Complete sequential ray tracing algorithm
- 4.8 Complete GPU version static ray tracing
- 4.12 Support light reflection and refraction
- 4.17 Support real-time rendering
- 4.25 Optimize our algorithm and hopefully achieve fps requirements.
- 4.29 Finish detailed performance analysis and assignment report.


## Milestone
In the past few weeks, we made progress in the following directions:
- Study the theoretical knowledge of ray tracing and optical physics.
- Find reference sequential implementation and test its performance on GHC machine.
- Implement the CUDA version of the static ray tracing algorithm and conduct experiments to evaluate its performance and speedup.
- Exploring more optimization techniques on accelerating GPU ray tracing.

We believe we can accomplish our goals. Since we found aboundant information on this topic, and got preliminary results on cuda static ray tracing. The next steps will be how to speedup the ray tracer to support real-time rendering. We have found some potential optimization opportunities to achieve our performance goal:
- **Replacing recursion with iterative loop.** This can avoid the complexity of the recursive stack during ray tracing, which can cause unnecessary pressure on registers and local memory
- **Using shared memory to cache surface data.** If an object is intersected with multiple rays, the surface data of the object will be requested multiple times by each ray. Using shared memory can avoid the cost of fetching a large amount of data from global memory.
- **Quantization.** As the speed of GPU units is highly influenced by the precision of the operands, we will try to accelerate the speed by using single-precision float point.
- **Intersection Optimization.** Using spatial partitioning data structures include BSP trees, quadtrees, octrees, and k-D trees to reduce the number of surfaces that we need to check for each ray.
- **Minimize data movement cost.** The image data can be compressed to 8-bit components before sending it back to the host memory. This can minimize the amount of data being moved between the GPU and the host.

In our poster session, we will have both a real-time ray tracing demo and performance analysis graphs. Idealy, we will impelment an interactive ray tracer which users can design a scene and assign textures to objects by themselves. 

### Preliminary Results
Sequential, Static CPU version:

Platform: 2.6 GHz 6-Core Intel Core i7, 16 GB 2400 MHz DDR4, MAC OSX

| num of objects      | Sequential Time(s) | GPU Time(s)|
| ----------- | ----------- | ----------- |
| 576      | 400.1       | 3.13 |
| 220   | 175.2        | - |
| 100   | 83.8        | - |
| 10   | 17.0        | - |

### Concerns
- Consider current static sequential Ray tracer takes 400s to render a 960P image. We are not sure if we can optimize the cuda rendering version that support 30FPS realtime rendering. 
- Not sure about the workload and time to build an interactive interface.
- Still have to investigate collision detection and reflection algorithm on other 3D shapes other than sphere.