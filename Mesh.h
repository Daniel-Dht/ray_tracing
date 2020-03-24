#ifndef MESH_H
#define MESH_H

#include <string>
#include <map>
#include <limits>
#include <vector>
#include <iostream>
#include <list>

#include "Vector.h"
#include "Object.h"
#include "Triangle.h"
#include "Material.h"
#include "BBox.h"

using namespace std;



class Face
{
public:

    int _v[3] = {-1,-1,-1};

    //Vector3d N;
    Face(){};
    Face(int v1, int v2, int v3) {
        _v[0]=v1;_v[1]=v2;_v[2]=v3;
    }
    // get
    double v1() const { return _v[0]; }
    double v2() const { return _v[1]; }
    double v3() const { return _v[2]; }
  
    friend ostream& operator<<(ostream& os, const Face& dt);
};

class Mesh : public Object {
public:

    BBox bbox;
	BVH  bvh;
    int nb_box = 0;
    std::vector<Face> faces;
	std::vector<Vector3d> vertices;
	std::vector<Vector3d> normals;
    std::vector<Vector3d> normalsVert;
    int Nv, Nf;

	Mesh() {};
	Mesh(const char* obj, double scaling, const Vector3d& offset, Material mat, Vector3f rot=Vector3f(0,0,0)) : Object(mat) {

		readOBJ(obj);	
        

        
        std::vector<Vector3d> tts;
        tts.reserve(Nv);
        // for (int i = 0; i < Nv; i++) {
        //     tts.push_back( vertices[i] - Vector3d(0));
        //     vertices[i] = tts[i]*(-1);
        // }
        rotate(rot);
		for (int i = 0; i < Nv; i++) {            
			vertices[i] = vertices[i] * scaling + offset ;
		}	

        computeFacesNormal();
        computeVertNormal();

        create_BBox();
        create_BVH(&bvh, 0, Nf);
        cout << "nombre de box dans le BVH : " << nb_box << endl;
	} 
    
	double intersect0(const Ray& r, int& total_ray_casted, InterStruct &interStruct) const {
		
        if (!bbox.intersect(r)) { 
        // if (!bvh.bbox.intersect(r))
             return -1;
        }
        
        //return intersectBVH(r, total_ray_casted, interStruct, &bvh);        

		double min_t = 1E99;
        InterStruct interStruct_tri;
		for(int i =0 ; i < faces.size(); i++) // faces.size()
		{
			const Face& face = faces[i];
			const Vector3d &v1 = vertices[face.v1()];
			const Vector3d &v2 = vertices[face.v2()];
			const Vector3d &v3 = vertices[face.v3()];

			Triangle tri = Triangle(v1, v2, v3, mat, false, normals[i]);
			const double t = tri.intersect(r, total_ray_casted, interStruct_tri);

    		if (t>0 && t<min_t){
	        	min_t = t;
                interStruct = interStruct_tri;							
	        }
		}

        // InterStruct interStructBVH;
        // double minBVH =  intersectBVH(r, total_ray_casted, interStructBVH, &bvh);

        
        // if(min_t != minBVH){
        //     cout << "** N "<< interStruct.N << " , et BVH   " << interStructBVH.N <<endl;
        // }
        // if(min_t != minBVH){
        //     cout << "** P "<< interStruct.P << " , et BVH   " << interStructBVH.P <<endl;
        // }
        // //if(min_t <1E99){
        // if(min_t != minBVH){
        //     cout << min_t << " , et BVH   " << minBVH<<endl;
        // }
		return min_t;		
	}

    double intersectBVH(const Ray& r, int& total_ray_casted, InterStruct &interStruct, const BVH *current) const {
        if (!current->bbox.intersect(r))
            return 1E99;

        if (!current->bl){
      
            double min_t = 1E99;
            InterStruct interStruct_tri;
            for(int i =current->i0 ; i < current->i1; i++) // faces.size()
            {
                const Face& face = faces[i];
                const Vector3d &v1 = vertices[face.v1()];
                const Vector3d &v2 = vertices[face.v2()];
                const Vector3d &v3 = vertices[face.v3()];

                Triangle tri = Triangle(v1, v2, v3, mat, false, normals[i]);
                const double t = tri.intersect(r, total_ray_casted, interStruct_tri);

                if (t>0 && t<min_t){
                    min_t = t;
                    interStruct = interStruct_tri;							
                }
            }
            return min_t;
        } 

        double t_bl =  intersectBVH(r, total_ray_casted, interStruct, current->bl);
        double t_br =  intersectBVH(r, total_ray_casted, interStruct, current->br);
        return fmin(t_bl,t_br);

    }
    double intersect(const Ray& r, int& total_ray_casted, InterStruct &interStruct) const {
		
        double min_t = 1E99;

        //if (!bbox.intersect(r))
        if (!bvh.bbox.intersect(r))
            return min_t;
        
        std::list<const BVH*> l;
        l.push_front(&bvh);
        
        while(!l.empty()) {

            const BVH *current = l.front();
            l.pop_front();
            if (current->bl && current->bl->bbox.intersect(r)) {
                l.push_back(current->bl);
            }
            if (current->br && current->br->bbox.intersect(r)) {
                l.push_back(current->br);
            }
            if (!current->bl){
                InterStruct interStruct_tri;
                for(int i =current->i0 ; i < current->i1; i++) // faces.size()
                {
                    const Face& face = faces[i];
                    const Vector3d &v1 = vertices[face.v1()];
                    const Vector3d &v2 = vertices[face.v2()];
                    const Vector3d &v3 = vertices[face.v3()];
                    
                    //Triangle tri = Triangle(v1, v2, v3, mat, false, normals[i]);
                    Triangle tri = Triangle(v1, v2, v3, normals[i]);

                    double t;
                    if(interStruct.computeN && 1 ){
                        const Vector3d &n1 = normalsVert[face.v1()];
                        const Vector3d &n2 = normalsVert[face.v2()];
                        const Vector3d &n3 = normalsVert[face.v3()];
                        t = tri.intersect_interpolateN(r, total_ray_casted, interStruct_tri, n2, n3, n1, min_t);
                    } else {
                        t = tri.intersect(r, total_ray_casted, interStruct_tri);
                    }            

                    if (t>0 && t<min_t){
                        min_t = t;
                        interStruct = interStruct_tri;							
                    }
                }
            }
        }
		return min_t;		
	}
    
    void computeFacesNormal() {
        normals.clear();
        normals.reserve(Nf);
        for (int i=0 ; i < Nf; i++)
        {
            Face &f = faces[i];

            Vector3d A = vertices[f.v1()]-vertices[f.v2()];
            Vector3d B = vertices[f.v1()]-vertices[f.v3()];
            normals.push_back( B.cross(A).getNormalized()*(1.) );
        }
    }
    void computeVertNormal() {
        normalsVert.clear();
        normalsVert.reserve(Nv);

        // certains vertex ont leur normal nulle après calcul, 
        // si c'est le cas on leur attribut la normal d'une de leur face incicente
        std::vector<int> faceofVert;
        faceofVert.reserve(Nv);

        for (int i=0 ; i < Nv; i++){ 
            normalsVert.push_back(Vector3d(0));
            faceofVert.push_back(-1);
        }
        
        for (int i=0 ; i < Nf; i++)
        {
            const Face &f = faces[i];
            normalsVert[f.v1()] += normals[i];
            normalsVert[f.v2()] += normals[i];
            normalsVert[f.v3()] += normals[i];    
            if( faceofVert[f.v1()]==-1) faceofVert[f.v1()] = f.v1();
            if( faceofVert[f.v2()]==-1) faceofVert[f.v2()] = f.v2();
            if( faceofVert[f.v3()]==-1) faceofVert[f.v3()] = f.v3();
        }
        for (int i=0 ; i < Nv; i++){ 
            if( normalsVert[i].norm()==0) {
                // cout << i << " CACA" << endl;
                normalsVert[i] = normals[faceofVert[i]];
            }
            else normalsVert[i].normalize();
        }

        // DEBUG
        // for (int i=0 ; i < Nf; i++)
        // {
        //     const Face &f = faces[i];
        //     normals[i] += normalsVert[f.v1()] ;
        //     normals[i] += normalsVert[f.v2()] ;
        //     normals[i] += normalsVert[f.v3()] ;   
        //     normals[i].normalize();         
        // }
    } 

    void readOBJ(const char file_name[]) {
        FILE *fp;
        fp = fopen(file_name, "r");  //Ouverture d'un fichier en lecture
        if(fp == NULL)
        {
            cout <<"Error opening file: "<< file_name << endl;
            //exit(1);
            return;
        }
        cout <<"Opening "<< file_name << endl;

        vertices.clear();
        faces.clear();
        int nb_edge = 0;
        fscanf(fp, "%d %d %d\n", &Nv, &Nf, &nb_edge);
        cout <<"nb vertices: "<< Nv<< ", nb faces: " << Nf  << endl;

        vertices.reserve(Nv);
        faces.reserve(Nf);

        double x, y, z;
        for(int i_vertex = 0; i_vertex < Nv; i_vertex++)
        {
            fscanf(fp, "%lf %lf %lf\n", &x, &y, &z);
            vertices.push_back(Vector3d(x,z,y));
        }

        int n_face, v1, v2, v3;
        for(int i_triangle = 0; i_triangle < Nf; i_triangle++)
        {
            fscanf(fp, "%d %d %d %d\n", &n_face, &v1, &v2, &v3);
            faces.push_back(Face(v1, v2, v3));
        }
        fclose(fp);
    }

    void create_BBox(){

        Vector3d c_min = vertices[0];
        Vector3d c_max = vertices[0];
        for (int i = 1; i < vertices.size(); ++i)
        {
            c_min.x = fmin(c_min.x, vertices[i].x);
            c_max.x = fmax(c_max.x, vertices[i].x);

            c_min.y = fmin(c_min.y, vertices[i].y);
            c_max.y = fmax(c_max.y, vertices[i].y);

            c_min.z = fmin(c_min.z, vertices[i].z);
            c_max.z = fmax(c_max.z, vertices[i].z);
        }
        bbox = BBox(c_min, c_max);
    };

    BBox create_BBox(int i0, int i1){ 
        // i0 and i1 are indices in the faces vector

        Vector3d c_min = vertices[faces[i1].v1()];
        Vector3d c_max = vertices[faces[i1].v1()];
        
        for (int i = i0; i < i1; ++i) {

            const Face &f = faces[i];

            for(const int id_vert : f._v){ // the 3 vertice indices of the face f

                c_min.x = fmin(c_min.x, vertices[id_vert].x);
                c_max.x = fmax(c_max.x, vertices[id_vert].x);

                c_min.y = fmin(c_min.y, vertices[id_vert].y);
                c_max.y = fmax(c_max.y, vertices[id_vert].y);

                c_min.z = fmin(c_min.z, vertices[id_vert].z);
                c_max.z = fmax(c_max.z, vertices[id_vert].z);           
            }           
        }
        return BBox(c_min, c_max);
    }

    void create_BVH(BVH *node, int i0, int i1) {
        // i0 and i1 are indices in the faces vector

        nb_box++;
        node->bbox = create_BBox(i0, i1);
        node->i0 = i0;
        node->i1 = i1;

        const Vector3d diag = node->bbox.c_max - node->bbox.c_min;
        int split_dim; // dimension sur laquelle la bbox est la plus longue
        if      (diag.x>diag.y && diag.x>diag.z) split_dim = 0; // x
        else if (diag.y>diag.x && diag.y>diag.z) split_dim = 1; // z
        else split_dim = 2;

        // valeur de la longueur à laquelle la BBox est coupée en 2
        const double split_val = node->bbox.c_min[split_dim] + diag[split_dim] / 2;
        
        // index on which the face will be swap
        int pivot = i0-1;
        // take all face under the split_value, and place them at the beggining of the list.
        for(int i=i0 ; i< i1 ;i++) {

            double center_face_split_dim = (vertices[faces[i].v1()][split_dim] + vertices[faces[i].v2()][split_dim] + vertices[faces[i].v3()][split_dim]) /3. ;

            //cout << center_face_split_dim << " " <<split_val<<endl;

            if ( center_face_split_dim < split_val){
                pivot++;                
                std::swap(faces[pivot], faces[i]);
                std::swap(normals[pivot], normals[i]);
            }
        }
        
        // si tous les triangles restant ne sont que d'un côté:
        if (pivot<=i0 || pivot>=i1) return;
        // si il ya trop peu de triangles dans uen des prochaines sous boite:
        if (fmin(pivot-i0, i1-pivot) < 2) return;

        node->bl = new BVH();
        node->br = new BVH();
        create_BVH(node->bl, i0, pivot);
        create_BVH(node->br, pivot, i1);
      
    }
    void rotate(Vector3f rot) {
        double cosa = cos(rot.x);
        double sina = sin(rot.x);

        double cosb = cos(rot.y);
        double sinb = sin(rot.y);

        double cosc = cos(rot.z);
        double sinc = sin(rot.z);

        double Axx = cosa*cosb;
        double Axy = cosa*sinb*sinc - sina*cosc;
        double Axz = cosa*sinb*cosc + sina*sinc;

        double Ayx = sina*cosb;
        double Ayy = sina*sinb*sinc + cosa*cosc;
        double Ayz = sina*sinb*cosc - cosa*sinc;

        double Azx = -sinb;
        double Azy = cosb*sinc;
        double Azz = cosb*cosc;

        for (int i = 0; i < Nv; i++) {
            double px = vertices[i].x;
            double py = vertices[i].y;
            double pz = vertices[i].z;

            vertices[i].x = Axx*px + Axy*py + Axz*pz;
            vertices[i].y = Ayx*px + Ayy*py + Ayz*pz;
            vertices[i].z = Azx*px + Azy*py + Azz*pz;
        }
    }
}; 


#endif // MESH_H