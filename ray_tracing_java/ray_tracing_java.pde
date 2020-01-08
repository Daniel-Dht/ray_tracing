
float fov = PI/3;
PVector u = null;

void setup(){
  size(400,400);
}

void draw() {
  noLoop(); 
  
  float d = width / (2*tan(fov/2.0));  
  PVector C = new PVector(width/2, height/2, fov);
  Sphere s = new Sphere(new PVector(width/2, height/2,0) , 10.0);
  
  for(int i =0; i<width; i++){
    for(int j =0; j<height; j++){
      
         
      u = new PVector(i-width/2, j-height/2,  -d);
      u.normalize();
      
      Ray ray = new Ray(C, u);
      
      boolean intersec = s.intersect(ray);
      if (intersec) {
        stroke(0,255,0);
        point(i,j);
      } 
    }
  }
}
