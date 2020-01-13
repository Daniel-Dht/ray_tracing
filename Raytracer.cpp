
// one include path:
// C:\Program Files (x86)\Windows Kits\10\Include\10.0.18362.0\shared

//->  g++ Raytracer.cpp -std=c++17 -fopenmp -O1
// g++ arg : https://linux.die.net/man/1/g++,  -ftime-report for cimpile time
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "Ray.h"
#include "Sphere.h"
#include "Scene.h"
#include "Light.h"

#include <chrono> 

#include "iostream" 
#include <Eigen/Dense>

using namespace std;
using namespace std::chrono; 

using Eigen::Vector3d; // norme / squaredNorm / normalize / dot / +,-,/,* /
using Eigen::Vector3f;


int main() {

    auto start = high_resolution_clock::now();    

    const int W = 800;
    const int H = 600;
    const int WH = W*H;

    const double PI = 3.14159265359;
    const double fov = PI/3; 
    const double z = - W / (2*tan(fov/2.0));
  
    const Vector3d C = Vector3d(0.0, 0.0, 60); // camera origin

    Scene scene = Scene();
        
    std::vector<unsigned char> image(W*H * 3, 0);

    int i,j;
    #pragma omp parallel for schedule(dynamic,8)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
        // for (int index = 0; index < WH; ++index)
        // {
        //     i = int(index/H);
        //     j = int(index-W*i);

            Vector3d u = Vector3d(j-W/2.0+0.5, -i+H/2.0+0.5, z); // DIRECTION DU RAYON
            u.normalize();        

            Vector3f intersec_col = scene.intersect(Ray(C, u));
            //if (intersec_col[0] != -1.0) {
                float index = (i*W + j) * 3;
                image[index + 0] = intersec_col[0];
                image[index + 1] = intersec_col[1];
                image[index + 2] = intersec_col[2];                
            //}             
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    auto stop = high_resolution_clock::now();
    auto cduration = duration_cast<milliseconds>(stop - start);
    cout << "compute time: "
    << cduration.count() << " milliseconds" << endl; 
    return 0;
}

// l'intensité lumineuse n'est pas perçu par l'oeil comme multplié par 2 lorque on la multiplié par 2. A la place on multipli par i^0.45.