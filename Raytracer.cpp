// one include path:
// C:\Program Files (x86)\Windows Kits\10\Include\10.0.18362.0\shared

//->  g++ Raytracer.cpp -std=c++17 -fopenmp -O1
// g++ arg : https://linux.die.net/man/1/g++,  -ftime-report for compile time

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "Ray.h"
#include "Scene.h"
#include "Vector.h"
#include "ProgressBar.h"
#include "Triangle.h"

#include <algorithm> 
#include <chrono> 
#include <iostream>

using namespace std;
using namespace std::chrono; 

#define M_PI 3.14159265359
#define M_2PI 6.28318530718

int main(int argc, char **argv) {
    int max_bounce = 10;
    int n_rays = 1;

    if (argc < 3)
    {
        fprintf(stderr, "Error: you must provide at least width and height (n_rays or recursive depth are options, default:1 & 10)\n");
        return -1;
    }
    if (argc == 4) n_rays = atoi(argv[3]);
    if (argc == 5) {
        n_rays = atoi(argv[3]);
        max_bounce = atoi(argv[4]);
    }

    auto start = high_resolution_clock::now();        

    const int W = atoi(argv[1]);
    const int H = atoi(argv[2]);
    const int WH = W*H;
    
    Geometry test42 = Geometry("cube.off", 2, Vector3d(0,0,0));

    const double fov = M_PI/3; 
    const double z = - W / (2*tan(fov/2.0));
  
    const Vector3d C = Vector3d(0.0, 0.0, 60); // camera origin

    Scene scene = Scene(n_rays, max_bounce);
        
    std::vector<unsigned char> image(W*H * 3, 0);

    ProgressBar progBar = ProgressBar(H, 10);
    #pragma omp parallel for schedule(static,4)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {                   

            Vector3f intersec_col = Vector3f(0,0,0);   
            for (int k = 0; k < scene.n_rays; ++k) {   
                // antialiasing                
                double r1 = uniform(engine);
			    double r2 = uniform(engine); 
                double R = 0.25*sqrt(-2*log(r1));
                double dx = R*cos(M_2PI*r2);
                double dy = R*sin(M_2PI*r2);
                //dx = dy = 0;
                Vector3d u = Vector3d(j-W/2.0+0.5+dx, -i+H/2.0+0.5+dy, z).getNormalized(); // DIRECTION DU RAYON            

                // casting ray
                intersec_col += scene.intersect(Ray(C, u));
            }  

            intersec_col /= scene.n_rays;

            float index = (i*W + j) * 3;
            image[index + 0] = std::clamp(std::pow(intersec_col[0], 0.45), 0.0, 255.0);
            image[index + 1] = std::clamp(std::pow(intersec_col[1], 0.45), 0.0, 255.0);
            image[index + 2] = std::clamp(std::pow(intersec_col[2], 0.45), 0.0, 255.0);                
        }
        if( i%progBar.modulo == 0) progBar.cout_progress(float(i));                            
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    auto cduration = duration_cast<milliseconds>(high_resolution_clock::now() - start);
    cout << endl;
    cout << "compute time: " << cduration.count() << " milliseconds" << endl; 
    cout << "total ray casted: " << scene.total_ray_casted / 1E6 << " million" <<endl; 

    return 0;
}
