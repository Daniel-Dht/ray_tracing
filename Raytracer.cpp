
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "Ray.h"
#include "Sphere.h"
#include "Scene.h"
#include "Light.h"

#include "iostream" 
#include <Eigen/Dense>
 
using Eigen::Vector3d;
using Eigen::Vector3f;

int main() {

    const int W = 512;
    const int H = 512;

    const double PI = 3.14159265359;
    const double fov = PI/3; 
    const double z = - W / (2*tan(fov/2.0));
  
    const Vector3d C = Vector3d(0.0, 0.0, 60); // camera origin

    Scene scene = Scene();
        
    std::vector<unsigned char> image(W*H * 3, 0);

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
   
            Vector3d u = Vector3d(j-W/2.0+0.5, -i+H/2.0+0.5, z); // DIRECTION DU RAYON
            u.normalize();

            Ray ray = Ray(C, u);

            //float intersec = scene.spheres[0]->intersect(ray);
            Vector3f intersec_col = scene.intersect(ray);
            if (intersec_col[0] != -1.0) {
                float index = (i*W + j) * 3;
                image[index + 0] = intersec_col[0];
                image[index + 1] = intersec_col[1];
                image[index + 2] = intersec_col[2];                
            }             
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    return 0;
}

// l'intensité lumineuse n'est pas perçu par l'oeil comme multplié par 2 lorque on la multiplié par 2. A la place on multipli par i^0.45.