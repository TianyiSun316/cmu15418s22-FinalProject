# cmu15418s22-FinalProject



## SUMMARY


## BACKGROUND


## THE CHALLENGE


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
