#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#define MAX_SPHERES 5
#define MAX_TRIANGLES 20
#define MAX_LIGHTS 5
#define MAX_TEXTURE 11
#define MAX_NAME_LENGTH 256
#define MAX_TEXTURE_SIZE 2000


typedef struct{
    float x, y, z;
}VectorType;

typedef struct{
    float r, g, b, eta; 
}ColorType;

typedef struct{
    float n1,n2,n3;
}Normal;

typedef struct{
    float x,y,z;
}VertexType;

typedef struct{
    float u,v;
}TextureCoordType;

typedef struct{ //struct to store each individual texture image with its own 2d array
    ColorType **pixelColor;
    int text_idx;
    char *file_name;
    int height,width;
}TexturePic;

typedef struct{
    int v[3], vn[3],vt[3];
    float normal[3];
    float smoothF;
    int c;
    int textF;
}TriangleType;

typedef struct{
    VectorType diffuse;
    VectorType specular;
    float ka;
    float kd;
    float ks;
    float n;
    float alpha;
    float eta;
}MtlType;

typedef struct{
    float x,y,z,w;
}LightDir;

typedef struct{
    LightDir dir;
    float i;
}LightType;

typedef struct{
    float x, y, z;
    float r;
    int c;
    int textF;
    
} SphereType;

typedef struct{
    float x, y, z;
    float dx, dy, dz; 
}RayType;

typedef struct{
    int w;
    int h;
}ImsizeType;
ColorType Trace_ray(RayType ray, SphereType spheres[MAX_SPHERES], MtlType color_list[10],
    ColorType back_color, LightType light[MAX_LIGHTS], int light_count,
    TriangleType triangles[MAX_TRIANGLES], VertexType vertices[20],
    VertexType verticesN[20], TexturePic textureFiles[MAX_TEXTURE],
    TextureCoordType coords[MAX_TEXTURE], int depth);


// Functions to help compute operations for vectors 
float vector_length(VectorType v1){ 
    float length = sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
    return length;
}

VectorType cross_product(VectorType v1, VectorType v2){
    // Computes the cross product of two vectors
    VectorType cross;
    cross.x = (v1.y *v2.z) - (v1.z * v2.y);
    cross.y = (v1.z * v2.x) - (v1.x * v2.z);
    cross.z = (v1.x * v2.y) - (v1.y * v2.x);
    return cross;                             

}


VectorType normalize_vector(VectorType v1){
    float length = vector_length(v1);
    if(length > 0.0){ 
        // Only want to do if length > 0 to avoid errors
        v1.x = v1.x/length;
        v1.y = v1.y/length;
        v1.z = v1.z/length;
        return v1;
    }
    VectorType zero;
    zero.x = 0.0;
    zero.y = 0.0;
    zero.z = 0.0;
    return zero ;
}

VectorType scale_vector(VectorType v1, float s){
    // Method to scale a vector with given scalar
    v1.x = v1.x * s;
    v1.y = v1.y * s;
    v1.z = v1.z * s;
    return v1;
}

float dot_product(VectorType v1, VectorType v2){
    float x = v1.x * v2.x;
    float y = v1.y * v2.y;
    float z = v1.z * v2.z;
    float w = x + y + z;
    return w; 
}



RayType make_ray(VectorType origin, VectorType point){
    // Takes in origin and point to make ray 
    RayType ray;
    ray.x = origin.x;
    ray.y = origin.y;
    ray.z = origin.z;
    ray.dx = point.x;
    ray.dy = point.y;
    ray.dz = point.z;
    return ray;
}
TexturePic load_ppm(char *texturFileName){
    // Helper Function for storing texture image into a 2D array define by TexturePic Struct
    if(texturFileName == NULL){
        printf("invalid file name");
    }
    TexturePic texture; // Where image will be stored
    texture.file_name = texturFileName;
    
    FILE *fp = fopen(texture.file_name, "r");
    if (fp == NULL) {
        perror("Error opening texture file");
        exit(1);
    }
    // Reading ppm file header, checking for valid format
    char header[3];
    int max_color;
    if (fscanf(fp, "%2s %d %d %d", header, &texture.width, &texture.height, &max_color) != 4) {
        // Check if correct header
        fclose(fp);
        exit(1);
    }
    
    texture.pixelColor = malloc(texture.height * sizeof(ColorType *));
    if (texture.pixelColor==NULL) {
        fclose(fp);
        exit(1);
    }
    
    for (int i = 0; i < texture.height; i++) {
        texture.pixelColor[i] = malloc(texture.width * sizeof(ColorType));
        if (texture.pixelColor[i] == NULL) {
            fclose(fp);
            exit(1);
        }
    }
    // Loops for storing image r,g,b values into texture struct
    for (int i = 0; i < texture.height; i++) {
        for (int j = 0; j < texture.width; j++) {
           if (fscanf(fp, "%f %f %f", 
               &texture.pixelColor[i][j].r, 
               &texture.pixelColor[i][j].g, 
               &texture.pixelColor[i][j].b) != 3) {
                fclose(fp);
                exit(1);
            }
        }
    }
    fclose(fp);
    return texture; 
}

void free_textures(TexturePic textures[MAX_TEXTURE]){
    for (int i = 0; i < MAX_TEXTURE; i++) {
        for (int j = 0; j < textures[i].height; j++) {
            free(textures[i].pixelColor[j]);
        }
        free(textures[i].pixelColor); 
    }
    free(textures);
}

ColorType texture_coord_sphere(SphereType sphere, TexturePic textureFiles[MAX_TEXTURE], VectorType intersection){
    float u; // pheta/ 2pi
    float v; // o / pi 
    float pheta;//0
    pheta = atan2(intersection.y - sphere.y, intersection.x - sphere.x);
    float phi;//o
    phi = acos((intersection.z - sphere.z)/ sphere.r);
   
    if(pheta > 0){
       // then u = pheta/2pi 
       u = (pheta/ (2 * M_PI));
    }
    else if(pheta< 0){
        u = (pheta + (2 * M_PI))/ (2 * M_PI);
       //u = (pheta + 2pi) / 2pi
    }
    v = phi / M_PI;
    TexturePic texture1 = textureFiles[sphere.textF];
    int height = texture1.height;
    int width = texture1.width;


    ColorType text_color;
    int text_Y = (int)(v * height);
    int text_X = (int)(u * width);
    text_color.r = texture1.pixelColor[text_Y][text_X].r;
    text_color.g = texture1.pixelColor[text_Y][text_X].g;
    text_color.b = texture1.pixelColor[text_Y][text_X].b;
    text_color.r = text_color.r/255;
    text_color.b = text_color.b/255;
    text_color.g = text_color.g/255;
    return text_color;
}
ColorType texture_coord_triangle(TriangleType triangle, TexturePic textureFiles[MAX_TEXTURE], TextureCoordType coords[10], float beta, float alpha, float gamma){
    // Takes in texture files array and correctly indexes into assignmed ppm file for Triangle
    // useing textF, then computes texture coordinate and gets coordinat form texture file image at 
    // that location
    float u1 = coords[triangle.vt[0]].u;
    float u2 = coords[triangle.vt[1]].u;
    float u3 = coords[triangle.vt[2]].u;
    float v1 = coords[triangle.vt[0]].v;
    float v2 = coords[triangle.vt[1]].v;
    float v3 = coords[triangle.vt[2]].v;
    
    float u = (alpha * u1) + (beta * u2) + (gamma * u3);
    float v = (alpha * v1) + (beta * v2) + (gamma * v3);
    TexturePic texture1 = textureFiles[triangle.textF];
    int height = texture1.height;
    int width = texture1.width;
    int text_X = (int)(u * (width - 1));
    int text_Y = (int)(v * (height - 1));
   
    ColorType text_color;
    text_color.r = texture1.pixelColor[text_Y][text_X].r;
    text_color.g = texture1.pixelColor[text_Y][text_X].g;
    text_color.b = texture1.pixelColor[text_Y][text_X].b;
    text_color.r = text_color.r/255;
    text_color.b = text_color.b/255;
    text_color.g = text_color.g/255;
   
    return text_color;
}
ColorType shade_ray(SphereType spheres[MAX_SPHERES],RayType ray,SphereType sphere,VectorType intersection, MtlType color_list[10],
    LightType light[MAX_LIGHTS], int light_count, TriangleType triangles[MAX_TRIANGLES], TriangleType triangle, VertexType vertices[20],
     VertexType verticesN[20], TexturePic textureFiles[MAX_TEXTURE], TextureCoordType coords[20], int depth, ColorType back_color){
    //computes the color at each point of the intersected object using blinn phong model
    //computes the ambient term, and diffuse and specular term and assigns each component to the color
    // calls trace ray when handling refraction and reflection after computing their rays
    int idx; //index of the object mtl color 
    float Fr; //for Frensel formula 
    float F0;
    ColorType col  = back_color;
    ColorType reflected_col = {0,0,0};
    if(depth>= 10){
        // Reached max recursion depth
      return back_color;
    }
    VectorType N; //initializ N vector
    if(sphere.c != -1){
        //if we hit a sphere its index wont be -1
        idx = sphere.c; 
        N.x = (intersection.x - sphere.x) / sphere.r;
        N.y = (intersection.y - sphere.y) / sphere.r;
        N.z = (intersection.z - sphere.z) / sphere.r;    
        N = normalize_vector(N);
        if(sphere.textF != -1 ){
            ColorType text_color;
            text_color = texture_coord_sphere(sphere, textureFiles, intersection);
            color_list[idx].diffuse.x  = text_color.r;
            color_list[idx].diffuse.y  = text_color.g;
            color_list[idx].diffuse.z  = text_color.b;
        }
    }
   
    if(triangle.c != -1){
        //if we hit a triangle its index wont be -1
        //Computing N vector for triangle
        idx = triangle.c;
        VertexType p0 = vertices[triangle.v[0]]; //v1
        VertexType p1 = vertices[triangle.v[1]]; //v2
        VertexType p2 = vertices[triangle.v[2]]; //v3
        VectorType e1; 
        e1.x = p1.x - p0.x;
        e1.y = p1.y - p0.y;
        e1.z = p1.z - p0.z;
        VectorType e2; 
        e2.x = p2.x - p0.x;
        e2.y = p2.y - p0.y;
        e2.z = p2.z - p0.z;
        N = cross_product(e1,e2);
        if(triangle.smoothF == 1){ 
            // Check if triangle has smooth shading
            // To get smooth shading use normal vertices to 
            // compute N vector
            // Reusing math from triangle intersection to get alpha, beta, and gamma
            VertexType n0 = verticesN[triangle.vn[0]]; 
            VertexType n1 = verticesN[triangle.vn[1]]; 
            VertexType n2 = verticesN[triangle.vn[2]]; 
            float A = N.x;
            float B = N.y;
            float C = N.z; 
            float D =  -((A * p0.x) + (B * p0.y) + (C * p0.z));
            int t;
            t = - ((A * ray.x) +( B * ray.y) + (C * ray.z) + D);
            int denominator = (A * ray.dx + B * ray.dy + C * ray.dz);
            t = t/(A * ray.dx + B * ray.dy + C * ray.dz);
            if(t > 0 ){
                VectorType pt;
                pt.x = ray.x + (t * ray.dx);
                pt.y = ray.y + (t * ray.dy);
                pt.z = ray.z + (t * ray.dz);
                VectorType ep;
                ep.x = pt.x - p0.x;
                ep.y = pt.y - p0.y;
                ep.z = pt.z - p0.z;
                float d11 = dot_product(e1,e1);
                float d12 = dot_product(e1,e2);
                float d22 = dot_product(e2,e2);
                float d1p = dot_product(e1,ep);
                float d2p = dot_product(e2,ep);
                float det = ((d11 * d22) - (d12 * d12));
            
                float beta = ((d22 * d1p) - (d12 * d2p))/ det;
                float gamma = ((d11 * d2p) - (d12 * d1p))/ det;
                float alpha = 1 - (beta + gamma); 
        
                N.x = n0.x * alpha + n1.x * beta + n2.x * gamma;
                N.y = n0.y * alpha + n1.y * beta + n2.y * gamma;
                N.z = n0.z * alpha + n1.z * beta + n2.z * gamma;
               //N = normalize_vector(N);
            } 
        }
       N = normalize_vector(N);
       if(triangle.textF != -1 ){
        ColorType text_color;
        VertexType p0 = vertices[triangle.v[0]]; //v1
        VertexType p1 = vertices[triangle.v[1]]; //v2
        VertexType p2 = vertices[triangle.v[2]]; //v3
        VectorType e1; 
        e1.x = p1.x - p0.x;
        e1.y = p1.y - p0.y;
        e1.z = p1.z - p0.z;
        VectorType e2;
        VectorType W;
        e2.x = p2.x - p0.x;
        e2.y = p2.y - p0.y;
        e2.z = p2.z - p0.z;
        
                VectorType pt;
                pt.x = intersection.x;
                pt.y = intersection.y;
                pt.z = intersection.z;
                VectorType ep;
                ep.x = pt.x - p0.x;
                ep.y = pt.y - p0.y;
                ep.z = pt.z - p0.z;
                float d11 = dot_product(e1,e1);
                float d12 = dot_product(e1,e2);
                float d22 = dot_product(e2,e2);
                float d1p = dot_product(e1,ep);
                float d2p = dot_product(e2,ep);
                float det = ((d11 * d22) - (d12 * d12));
            
                float beta = ((d22 * d1p) - (d12 * d2p))/ det;
                float gamma = ((d11 * d2p) - (d12 * d1p))/ det;
                float alpha = 1 - (beta + gamma); 
                    text_color = texture_coord_triangle(triangle, textureFiles,coords,beta,alpha,gamma);
                    color_list[triangle.c].diffuse.x  = text_color.r;
                    color_list[triangle.c].diffuse.y  = text_color.g;
                    color_list[triangle.c].diffuse.z  = text_color.b;
               
            }
        
    }
    
    ColorType color; //r is ambient, g is diffuse term, b is specular term
    VectorType ambient_term; //ka * Od
    
    ambient_term.x = color_list[idx].diffuse.x * color_list[idx].ka;
    ambient_term.y = color_list[idx].diffuse.y * color_list[idx].ka;
    ambient_term.z = color_list[idx].diffuse.z * color_list[idx].ka;

    VectorType diffuse_term = {0.0,0.0,0.0}; 
    VectorType specular_term = {0.0,0.0,0.0}; 
    VectorType L;
    VectorType H;
    VectorType V; 
    V.x = ray.x - intersection.x;
    V.y = ray.y - intersection.y;
    V.z = ray.z - intersection.z;

    V = normalize_vector(V);
    
    //Where work for assignment 1d starts 
    float IDotN = (ray.dx * N.x) + (ray.dy * N.y) + (ray.dz* N.z);
    VectorType Inc; // the incident ray
    Inc.x = ray.dx;
    Inc.y = ray.dy;
    Inc.z = ray.dz;
    float eta = color_list[idx].eta; //refraction index 
    float cos_theta_i =  IDotN;
    if(color_list[idx].ks > 0){
        // Only want to make reflection if ks >0
        VectorType reflected; 
        // Creating the reflected ray direction 
        reflected.x = 2*(IDotN) * N.x - ray.dx;
        reflected.y = 2*(IDotN) * N.y - ray.dy;
        reflected.z = 2*(IDotN) * N.z - ray.dz;
        reflected = normalize_vector(reflected);
        VectorType origin = intersection; 
        origin.x += Inc.x + reflected.x * 0.001;  
        origin.y += Inc.y + reflected.y * 0.001;  
        origin.z += Inc.z + reflected.z * 0.001; 
                
        RayType reflected_ray = make_ray(origin, reflected);
        // Making ray that is reflected to call traceray on  
        reflected_col = Trace_ray(reflected_ray,spheres, color_list,back_color,light,light_count,
            triangles,vertices,verticesN,textureFiles,coords,depth + 1);
        // Trace ray will return the color of the reflection
      
    }
    ColorType refracted_color;

    if(color_list[idx].alpha < 1) { 
        //dont do refraction on things that are fully opaque
        if (cos_theta_i < 0) { 
            // ray is entering object  
           // have to change refraction index to apply snells law 
            eta = 1.0/ color_list[idx].eta;
        } 
        // if ray is exiting eta stays the same 
            F0 = pow((eta - 1) / (eta + 1), 2); //Fresnel reflectance
            Fr = F0 + (1 - F0) * pow(1 - fabs(cos_theta_i), 5);
           
            float d = 1 - (eta * eta) * (1 - (cos_theta_i * cos_theta_i));
        if(d>=0){
            //If d is less than 0 that means the transmitted ray doesn't exist 
            VectorType transmitted;
        
            float cos_theta_t = sqrt(d);
            //use snell laws to get cos theta t, and compute T
            transmitted.x = eta * Inc.x + (eta * cos_theta_i - cos_theta_t) * N.x;
            transmitted.y = eta * Inc.y + (eta * cos_theta_i - cos_theta_t) * N.y;
            transmitted.z = eta * Inc.z + (eta * cos_theta_i - cos_theta_t) * N.z;
            transmitted = normalize_vector(transmitted);
           
            VectorType origin = intersection; 
            //adding offset so transmitted ray correclty hits object
            origin.x += transmitted.x * 0.001;  
            origin.y += transmitted.y * 0.001;  
            origin.z += transmitted.z * 0.001; 
            
           
            RayType transmitted_ray = make_ray(origin, transmitted); 
            // making transmitted ray   
            refracted_color = Trace_ray(transmitted_ray, spheres, color_list, back_color,
                 light, light_count, triangles, vertices, verticesN, textureFiles, coords, depth + 1);
        
           }
           else{
            // otherwise have total internal reflection so return reflected color computes earlier 
            refracted_color = reflected_col;
           }
            
             
    }
    // if it is directional light then w would be 0 if point light w is 1
    // Because need to handle multiple light sources and H depends on L as well
    // both are computed in light loop
    for(int i =0; i <light_count; i++){
        LightType cur_light = light[i];
        if(cur_light.dir.w ==0){
            //then L vector is a drectional light
            L.x = -cur_light.dir.x ;
            L.y = -cur_light.dir.y ;
            L.z = -cur_light.dir.z ;
            L = normalize_vector(L);
        }
        if(cur_light.dir.w ==1){
            //then L vector is a positional light
            L.x = cur_light.dir.x - intersection.x;
            L.y = cur_light.dir.y - intersection.y;
            L.z = cur_light.dir.z - intersection.z;
            L = normalize_vector(L); 
        }

        RayType shadow; 
        //ray that is casted to check for shadows 
        shadow.dx = L.x;
        shadow.dy = L.y;
        shadow.dz = L.z;
        shadow.x = intersection.x + (L.x * 0.001);
        shadow.y = intersection.y + (L.y * 0.001);
        shadow.z = intersection.z + (L.z * 0.001);


       int has_shadow = 0; //flag to track if there is a shadow
        float shadow_opt = 1.0; //shadow opacity flag added
        for (int j = 0; j < MAX_SPHERES; j++) { 
            //for loop to check if any spheres lie in the plane of light
            if (spheres[j].c == sphere.c) {
                //checking not the same sphere
                continue;
            }
            float a = (shadow.dx * shadow.dx) + (shadow.dy * shadow.dy) + (shadow.dz * shadow.dz);

            float b = 2 * ((shadow.dx * (shadow.x - spheres[j].x)) + (shadow.dy * (shadow.y - spheres[j].y)) 
            + (shadow.dz * (shadow.z - spheres[j].z)));

            float c = ((shadow.x - spheres[j].x) * (shadow.x - spheres[j].x)) + 
            ((shadow.y - spheres[j].y) * (shadow.y - spheres[j].y)) + ((shadow.z - spheres[j].z) * (shadow.z - spheres[j].z)) 
            - (spheres[j].r * spheres[j].r);

            float discriminant = (b * b) - (4 * a * c);
            
         
            
    
            if (discriminant >0) {  
                //only less than 0 if we hit background 
                if (color_list[spheres[j].c].alpha < 1.0) {
                    // transparent use shadow opt to allow some light
                    shadow_opt *= color_list[spheres[j].c].alpha;
                    
                } else {
                    has_shadow =1;
                    shadow_opt = 0;
                    break;
                }
               
            }
    
          
        }   
            int tri_has_shadow = 0; // shadow flag for triangles 
             for (int j = 0; j < MAX_TRIANGLES; j++) { 
                float smaller_t = INFINITY;
                //for loop to check if any spheres lie in the plane of light
                VertexType p0 = vertices[triangles[j].v[0]]; //v1
                VertexType p1 = vertices[triangles[j].v[1]]; //v2
                VertexType p2 = vertices[triangles[j].v[2]]; //v3
                VectorType e1; 
                e1.x = p1.x - p0.x;
                e1.y = p1.y - p0.y;
                e1.z = p1.z - p0.z;
                VectorType e2; 
                e2.x = p2.x - p0.x;
                e2.y = p2.y - p0.y;
                e2.z = p2.z - p0.z;
                VectorType n;
                n = cross_product(e1,e2);
                //use componeents of n to determmine A,B,C,D
                float A = n.x;
                float B = n.y;
                float C = n.z; 
                float D =  -((A * p0.x) + (B * p0.y) + (C * p0.z));
                //after finding D we have the eq x + y + z + D = 0 
                //substitute ray eq into plane equation 
                smaller_t = - ((A * shadow.x) +( B * shadow.y) + (C * shadow.z) + D);
                float denominator = (A * shadow.dx + B * shadow.dy + C * shadow.dz);
                if(denominator <=0){
                   //has_shadow = 0;
                    continue;
                }
               
                smaller_t = smaller_t/denominator;
                float max_t = sqrt(L.x * L.x + L.y * L.y + L.z * L.z); // Distance to the light

                if (smaller_t > 0 && (cur_light.dir.w == 1 && smaller_t < max_t)) {
                    has_shadow = 1;
                    //break;
                }
                
                if(smaller_t <0 ){
                    continue;
                }
                if(smaller_t > 0){
                VectorType pt;
                pt.x = shadow.x + (smaller_t * shadow.dx);
                pt.y = shadow.y + (smaller_t * shadow.dy);
                pt.z = shadow.z + (smaller_t * shadow.dz);
                //make sure the pt satisfies the plane equation
                //by finding baycentric coordinates gamma and beta
                VectorType ep;
                ep.x = pt.x - p0.x;
                ep.y = pt.y - p0.y;
                ep.z = pt.z - p0.z;
                float d11 = dot_product(e1,e1);
                float d12 = dot_product(e1,e2);
                float d22 = dot_product(e2,e2);
                float d1p = dot_product(e1,ep);
                float d2p = dot_product(e2,ep);
                float det = ((d11 * d22) - (d12 * d12));
                if(det == 0 ){
                    printf("no solution");
                    continue;
                }
                 
                float beta = ((d22 * d1p) - (d12 * d2p))/ det;
                float gamma = ((d11 * d2p) - (d12 * d1p))/ det;
                float alpha = 1 - (beta + gamma); 
                if (alpha > 0 && beta > 0 && gamma > 0 && alpha < 1 && beta < 1 && gamma < 1){
                    tri_has_shadow = 1;
                     break;
                }
                
            }
        }
        if(tri_has_shadow == 0  && has_shadow == 0){
            float NdotL = fmax(dot_product(N, L), 0.0);
            // Now with recursive ray tracing need to make H =(L + I / || L + I)
            //BEFORE:
            //  H.x = L.x + V.x;
            //  H.y = L.y + V.y;
            //  H.z = L.z + V.z;    
            //NOW:
            // H is halfway betwen L and the opposit of 
            H.x = L.x + (-ray.dx); 
            H.y = L.y + (-ray.dy);
            H.z = L.z + (-ray.dz);       
            H = normalize_vector(H);
           
            float NdotH = fmin(fmax(dot_product(N, H), 0.0),1.0);
            float specular =  pow(NdotH, color_list[idx].n);
            

            //after getting all information needed computing the diffuse and specular colors of objects
            // shadow opt added for tranperent object to allow light to come through
            diffuse_term.x += fmin(color_list[idx].diffuse.x * color_list[idx].kd * (NdotL)*shadow_opt, 1.0);
            diffuse_term.y += fmin(color_list[idx].diffuse.y * color_list[idx].kd * (NdotL)*shadow_opt, 1.0);
            diffuse_term.z += fmin(color_list[idx].diffuse.z * color_list[idx].kd * (NdotL)*shadow_opt,1.0);
            specular_term.x += fmin(color_list[idx].specular.x * color_list[idx].ks * specular* shadow_opt, 1.0);
            specular_term.y += fmin(color_list[idx].specular.y * color_list[idx].ks * specular*shadow_opt,1.0);
            specular_term.z += fmin(color_list[idx].specular.z * color_list[idx].ks * specular*shadow_opt,1.0);
        }
    }
    
            
    VectorType I;
    // Adding Fr * reflected  and refraction to illumination equation 
    //refraction isn't added if alpha =1 becuase object is fully opaq
    I.x = ambient_term.x + diffuse_term.x + specular_term.x + (Fr * reflected_col.r + (1 - Fr) * (1 - color_list[idx].alpha) * refracted_color.r);
    I.y = ambient_term.y + diffuse_term.y + specular_term.y + (Fr * reflected_col.g + (1 - Fr)  * (1 - color_list[idx].alpha)* refracted_color.g);
    I.z = ambient_term.z + diffuse_term.z + specular_term.z +  (Fr * reflected_col.b + (1 - Fr) * (1 - color_list[idx].alpha)* refracted_color.b);
   
  
  
    color.r = I.x;
    color.g = I.y;
    color.b = I.z;
    
    return color;
       
}





ColorType Trace_ray(RayType ray, SphereType spheres[MAX_SPHERES], MtlType color_list[10],
     ColorType back, LightType light[10], int light_count, TriangleType traingles[10], 
     VertexType vertices[20], VertexType verticesN[20], TexturePic textureFiles[20], TextureCoordType coords[20], int depth){
    // Returns color of object ray hits, if doesnt hit anything will return background color
    
    ColorType back_color;
    back_color.r = back.r; 
    back_color.g = back.g; 
    back_color.b = back.b; 
    if(depth >=10){
        //reached max recursion depth so just return back
        return back_color;
    }  
    float smallest_t = INFINITY;  // To keep track of the smallest t value
    SphereType closest; 
    TriangleType phony; 
    closest.c = -1; // If no sphere is hit will return index -1 
    phony.c = -1; // To track if closest object is Triangle

    for (int i = 0; i < MAX_SPHERES; i++) {
        //Looping through alll possible spheres in the image
        // To get the closest one the ray hits
       
        float a = (ray.dx * ray.dx) + (ray.dy * ray.dy) + (ray.dz * ray.dz);

        float b = 2 * ((ray.dx * (ray.x - spheres[i].x)) + (ray.dy * (ray.y - spheres[i].y)) 
        + (ray.dz * (ray.z - spheres[i].z)));

        float c = ((ray.x - spheres[i].x) * (ray.x - spheres[i].x)) + 
        ((ray.y - spheres[i].y) * (ray.y - spheres[i].y)) + ((ray.z - spheres[i].z) * (ray.z - spheres[i].z)) 
        - (spheres[i].r * spheres[i].r);

        float discriminant = (b * b) - (4 * a * c);
       

        if (discriminant < 0) { 
            // If less than zero we hit the backgorund
    
            continue;
        }

        if(discriminant > 0){
            // Two possibilities either hit sphere or grazed it 
            // If not zero we hit the sphere 
            float small_t1 = (-b - sqrt(discriminant)) / (2 * a);
            float small_t2 = (-b + sqrt(discriminant)) / (2 * a);
            // Check if t values are valid (positive) and if they are the smallest
            // check if smallest because we want the closest sphere
            if (small_t1 >0 && small_t1 < smallest_t) {
                smallest_t = small_t1;
                closest = spheres[i];
            }
            if (small_t2 > 0 && small_t2 < smallest_t) {
                smallest_t = small_t2;
                closest = spheres[i];
            }
        }
        if(discriminant ==0){
            // Discriminant is 0 meaning we grazed the sphere
            float small_t1 = (-b) / (2 * a);
            if (small_t1 > 0 && small_t1 < smallest_t) {

                smallest_t = small_t1;
                closest = spheres[i];
            }
        }
    }
    float smaller_t = INFINITY; 
    for(int i =0; i < MAX_TRIANGLES;  i++){
       VertexType p0 = vertices[traingles[i].v[0]]; //v1
       VertexType p1 = vertices[traingles[i].v[1]]; //v2
       VertexType p2 = vertices[traingles[i].v[2]]; //v3
       VectorType e1; 
       e1.x = p1.x - p0.x;
       e1.y = p1.y - p0.y;
       e1.z = p1.z - p0.z;
       VectorType e2; 
       e2.x = p2.x - p0.x;
       e2.y = p2.y - p0.y;
       e2.z = p2.z - p0.z;
       VectorType n;
       n = cross_product(e1,e2);
       //use componeents of n to determmine A,B,C,D
       float A = n.x;
       float B = n.y;
       float C = n.z; 
       float D =  -((A * p0.x) + (B * p0.y) + (C * p0.z));
       //after finding D we have the eq x + y + z + D = 0 
       //substitute ray eq into plane equation 
       smaller_t = - ((A * ray.x) +( B * ray.y) + (C * ray.z) + D);
       int denominator = (A * ray.dx + B * ray.dy + C * ray.dz);
       if(denominator == 0){
           continue;
       }
      
        smaller_t = smaller_t/(A * ray.dx + B * ray.dy + C * ray.dz);
        
        
       if(smaller_t > 0 ){
        VectorType pt;
        pt.x = ray.x + (smaller_t * ray.dx);
        pt.y = ray.y + (smaller_t * ray.dy);
        pt.z = ray.z + (smaller_t * ray.dz);
        //make sure the pt satisfies the plane equation
        //by finding baycentric coordinates gamma and beta
        VectorType ep;
        ep.x = pt.x - p0.x;
        ep.y = pt.y - p0.y;
        ep.z = pt.z - p0.z;
        float d11 = dot_product(e1,e1);
        float d12 = dot_product(e1,e2);
        float d22 = dot_product(e2,e2);
        float d1p = dot_product(e1,ep);
        float d2p = dot_product(e2,ep);
        float det = ((d11 * d22) - (d12 * d12));
        if(det == 0 ){
            printf("no solution");
            continue;
        }
        
        float beta = ((d22 * d1p) - (d12 * d2p))/ det;
        float gamma = ((d11 * d2p) - (d12 * d1p))/ det;
        float alpha = 1 - (beta + gamma); 
        if(alpha >= 0 && beta >= 0 && gamma >= 0 && alpha < 1 && beta <1 && gamma < 1){
            if(smaller_t < smallest_t){
                smallest_t = smaller_t;
                phony = traingles[i];
            }
        }
       }
    }
    // Closest sphere funtion will return the index of the spheres that is the closest
    SphereType hit = closest;
    float t = smallest_t;
    VectorType intersection_p;
    intersection_p.x = ray.x + (t * ray.dx); 
    intersection_p.y = ray.y + (t * ray.dy);
    intersection_p.z = ray.z + (t * ray.dz);

    if(hit.c != -1){ 
        phony.c = -1;
        // If closest sphere hits a sphere we get a real index, not -1
        // In case we dont hit a sphere but triangle hit.c will still be -1
        return shade_ray(spheres, ray, hit, intersection_p, color_list, light, light_count, traingles,phony,vertices,verticesN,textureFiles,coords,depth+1, back_color); // Get color form color list
    }
    if(phony.c != -1){
       
        hit.c =-1;
        // If we hit a sphere and not triangle phony will sitll be -1
        return shade_ray(spheres, ray,hit,intersection_p,color_list,light,light_count,traingles,phony, vertices,verticesN,textureFiles,coords,depth+1,back_color);
    }
   
    return back_color; 
}


 
int main(int argc, char *argv[]){ 
    int textCoord = 0;
    int text_file_num = 0;
    int vertex_c = 0;
    int vertexN_c = 0; //index into normals vertex
    int mtl_index = 0; //to keep track of index into material color array
    VertexType *vertices = (VertexType *)malloc(20 * sizeof(VertexType)); 
    VertexType *verticeN = (VertexType *)malloc(20 * sizeof(VertexType)); 
    MtlType *mtl_colors = malloc(10 * sizeof(MtlType));  // List to keep track of colors of objects DOESN'T hold background color 
    LightType *lights = malloc(10 * sizeof(LightType));//to have more than 1 light
    //char textureFileNames[MAX_TEXTURE][MAX_NAME_LENGTH]; //array to store texture file names  
    TexturePic *textureImages = (TexturePic*)malloc(MAX_TEXTURE * sizeof(TexturePic)); //array that stores texture images
    TextureCoordType *coords = (TextureCoordType*)malloc(20 * sizeof(TextureCoordType)); //to store texture coord vt
    if (textureImages == NULL) {
        return -1;
    }
    if (mtl_colors == NULL) {
        return -1;
    }
    if (lights == NULL) {
        return -1;
    }
    if(vertices == NULL){
        return -1;
    }
    if(verticeN == NULL){
        return -1;
    }
    if(coords == NULL){
        return -1;
    }
    char *file1 = argv[1];
    char *file2 = argv[1];
    FILE *fp = fopen(file1, "r");
    if(fp == NULL){
        printf("ERROR CANT OPEN FILE");
        return -1;

    }
    int text_flag = -1; //determines what item is using texture if -1 triangle if 0 sphere
    int sphere_count = 0; // To track number of spheres generated to not Exceed Max
    int light_count = 0; //Track number of lights on scene
    int triangle_count = 0;
    VectorType eye ;
    VectorType view_dir;
    VectorType up_dir;
    ImsizeType pic_size;
    ColorType back_ground_color;
    TriangleType *triangles = (TriangleType *)malloc(MAX_TRIANGLES * sizeof(TriangleType));
    SphereType *spheres = (SphereType *)malloc(MAX_SPHERES * sizeof(SphereType));
    // Making a list of spheres so scene can have more than one sphere
    if (spheres == NULL) {
        return -1;
    }
    if(triangles == NULL){
        return -1;
    }

    float v_fov; // Vertical Field of View
    float width;
    float height;
    
    char buffer[256]; 
    // To help read input File
    while(fgets(buffer, 256, fp)!=NULL){
        // While loop goes through input file
        // Using buffer and name to check for key words and storing information

       char *name = strtok(buffer, " "); 
       // Use buffer and name to break up input file 
       
       if(strcmp(name, "eye") == 0){

        eye.x = atof(strtok(NULL, " "));
        eye.y = atof(strtok(NULL, " "));
        eye.z = atof(strtok(NULL, " ")); 
       }
       if(strcmp(name, "viewdir") == 0){

        view_dir.x = atof(strtok(NULL, " "));
        view_dir.y = atof(strtok(NULL, " "));
        view_dir.z = atof(strtok(NULL, " ")); 
        
       }
       if(strcmp(name, "updir") == 0){

        up_dir.x = atof(strtok(NULL, " "));
        up_dir.y = atof(strtok(NULL, " "));
        up_dir.z = atof(strtok(NULL, " ")); 
        
       }
       if(strcmp(name, "vfov") == 0){

        v_fov = atof(strtok(NULL, " ")); 
        
       }
       if(strcmp(name, "imsize") == 0){
      
        width = atof(strtok(NULL, " ")) ;
        height = atof(strtok(NULL, " ")) ;
        pic_size.w = width ;
        pic_size.h = height ;
       
        
       }
       if(strcmp(name, "bkgcolor") == 0){

        float r = atof(strtok(NULL, " "));
        float g = atof(strtok(NULL, " "));
        float b = atof(strtok(NULL, " "));
        float e = atof(strtok(NULL, " "));
        back_ground_color.r = r;
        back_ground_color.g = g;
        back_ground_color.b = b;
        back_ground_color.eta = e;
        //should also include a value N
        //printf("%2f, this is the new value added\n", back_ground_color.eta);
        
        
       }
       if(strcmp(name, "mtlcolor") == 0){

        float odr = atof(strtok(NULL, " "));
        float odg = atof(strtok(NULL, " "));
        float odb = atof(strtok(NULL, " "));
        float osr = atof(strtok(NULL, " "));
        float osg = atof(strtok(NULL, " "));
        float osb = atof(strtok(NULL, " "));
        float ka = atof(strtok(NULL, " "));
        float kd = atof(strtok(NULL, " "));
        float ks = atof(strtok(NULL, " "));
        float highlight = atof(strtok(NULL, " "));
        float a = atof(strtok(NULL, " "));
        float e = atof(strtok(NULL, " "));
        mtl_colors[mtl_index].diffuse = (VectorType) {odr,odg,odb};
        mtl_colors[mtl_index].specular = (VectorType) {osr,osg,osb};
        mtl_colors[mtl_index].ka = ka;
        mtl_colors[mtl_index].kd = kd;
        mtl_colors[mtl_index].ks = ks;
        mtl_colors[mtl_index].n = highlight;
        mtl_colors[mtl_index].alpha = a;
        mtl_colors[mtl_index].eta = e;
        //Add the materials opacity a and its refractions N
        //printf("alpha : %2f, eta: %2f, this is the new value added\n", mtl_colors[mtl_index].alpha,mtl_colors[mtl_index].eta );

        //printf("%d",mtl_index);
        mtl_index++;
        } 
        
       if(strcmp(name, "sphere") == 0){
        sphere_count ++;
        
        if(text_file_num>= mtl_index && text_file_num != 0){
            spheres[sphere_count-1].textF = 2;
     
        }
        else{
            spheres[sphere_count-1].textF = -1;
        }
        
        //mtl_index ++; 
        if(sphere_count > MAX_SPHERES){
            printf("%s\n", "too many spheres, max number of spheres are 5");
            return -1;
        }
        float x = atof(strtok(NULL, " ")) ;
        float y = atof(strtok(NULL, " ")) ;
        float z = atof(strtok(NULL, " ")) ;
        float r = atof(strtok(NULL, " ")) ;
        if(text_file_num>= mtl_index && text_file_num != 0){
           spheres[sphere_count-1] = (SphereType) {x,y,z,r,mtl_index - 1, text_file_num-1};
            spheres[sphere_count-1].c = mtl_index-1;
            
        }
        else{
        spheres[sphere_count-1] = (SphereType) {x,y,z,r,mtl_index - 1, -1};
        spheres[sphere_count-1].c = mtl_index-1;
        }
        //soheres[sphere_count-1].textF
    
       }
       if(strcmp(name, "light") == 0){
            //multiple lights 
            float x = atof(strtok(NULL, " ")) ;
            float y = atof(strtok(NULL, " ")) ;
            float z = atof(strtok(NULL, " ")) ;
            float w = atof(strtok(NULL, " ")) ;
            float i = atof(strtok(NULL, " ")) ;
            light_count ++;
            if(light_count > 10){
                printf("%s\n", "too many lights, max number of lights are 10");
                return -1;
            }
            lights[light_count-1].dir = (LightDir) {x,y,z,w};
            lights[light_count-1].i = i;
        }
        if(strcmp(name, "v") == 0){
            vertex_c ++; //starts zero after first input should be 1;
            float x = atof(strtok(NULL, " "));
            vertices[vertex_c ].x = x;
            float y = atof(strtok(NULL, " "));
            vertices[vertex_c].y = y;
            float z = atof(strtok(NULL, " "));
            vertices[vertex_c].z = z;
        }
        if(strcmp(name, "vn") == 0){
            vertexN_c ++;
            float x = atof(strtok(NULL, " "));
            verticeN[vertexN_c ].x = x;
            float y = atof(strtok(NULL, " "));
            verticeN[vertexN_c].y = y;
            float z = atof(strtok(NULL, " "));
            verticeN[vertexN_c].z = z;
        }
        if(strcmp(name, "vt") == 0){
            textCoord++;
            float u = atof(strtok(NULL, " "));
            coords[textCoord].u = u;
            float v = atof(strtok(NULL, " "));
            coords[textCoord].v = v;
        }
        if(strcmp(name, "texture") == 0){
           char *text_file = strtok(NULL,"\n");    
            textureImages[text_file_num].text_idx = text_file_num;
            textureImages[text_file_num] = load_ppm(text_file); 
            text_file_num++;
          
        }
        
      
        if (strcmp(name, "f") == 0) {
            int v1, v2, v3;
            int vt1, vt2, vt3 ;
            int vn1, vn2, vn3;
            triangles[triangle_count].textF = -1;
            triangles[triangle_count].c = mtl_index-1;
            if(text_file_num >= mtl_index && text_file_num != 0 ){
                // Only if the number of mtl colors and texture are equal 
                // than that means this object has a texture, if multiple have same texture text num would stay same
                triangles[triangle_count].textF = text_file_num -1 ;
    
                
            }
            else{
                triangles[triangle_count].textF = -1;
            }
            buffer[strcspn(buffer, "\n")] = 0;
         
            char *token = strtok(NULL, " ");  // The first token should be "f"
            //token = strtok(NULL, " ");  // Skip "f", now token points to the first face component
            if (sscanf(token, "%d/%d/%d", &v1, &vt1, &vn1) == 3) {
                triangles[triangle_count].v[0] = v1;
                triangles[triangle_count].vt[0] = vt1;  
                triangles[triangle_count].vn[0] = vn1; 
                token = strtok(NULL, " ");
                sscanf(token, "%d/%d/%d", &v2, &vt2, &vn2);
                triangles[triangle_count].v[1] = v2;
                triangles[triangle_count].vt[1] = vt2;
                triangles[triangle_count].vn[1] = vn2;
                token = strtok(NULL, " ");
                sscanf(token, "%d/%d/%d", &v3, &vt3, &vn3);
                triangles[triangle_count].v[2] = v3;
               triangles[triangle_count].vt[2] = vt3;
                triangles[triangle_count].vn[2] = vn3;
             
                triangle_count++;
            } 
            else if (sscanf(token, "%d//%d", &v1, &vn1) == 2) {
                token = strtok(NULL, " ");
                sscanf(token, "%d//%d", &v2, &vn2);
                token = strtok(NULL, " ");
                sscanf(token, "%d//%d", &v3, &vn3);
                triangles[triangle_count].v[0] = v1;
                triangles[triangle_count].v[1] = v2;
                triangles[triangle_count].v[2] = v3;
                triangles[triangle_count].vn[0] = vn1;
                triangles[triangle_count].vn[1] = vn2;
                triangles[triangle_count].vn[2] = vn3;
                triangles[triangle_count].smoothF = 1;
                //mtl_index ++;
                triangle_count++;
            } 
            else if (sscanf(token, "%d/%d", &v1, &vt1) == 2) {
                // Success reading a face in v/t format
                token = strtok(NULL, " ");
                sscanf(token, "%d/%d", &v2, &vt2);
                token = strtok(NULL, " ");
                sscanf(token, "%d/%d", &v3, &vt3);
                triangles[triangle_count].v[0] = v1;
                triangles[triangle_count].v[1] = v2;
                triangles[triangle_count].v[2] = v3;
                triangles[triangle_count].vt[0] = vt1;
                triangles[triangle_count].vt[1] = vt2;
                triangles[triangle_count].vt[2] = vt3;
             
                triangle_count++;
            } 
            else if (sscanf(token, "%d", &v1) == 1) {
               
                token = strtok(NULL, " ");
                sscanf(token, "%d", &v2);
                token = strtok(NULL, " ");
                sscanf(token, "%d", &v3);
                triangles[triangle_count].v[0] = v1;
                triangles[triangle_count].v[1] = v2;
                triangles[triangle_count].v[2] = v3;
              
           
                triangle_count++;
            } 
            else {
                printf("Failed to read face data\n");
                return -1;
            }
        }
    }
    
    fclose(fp);
    // Changing file name to ppm
    char *fname = strchr(file2, '.');
    fname[1] = 'p';
    fname[2] = 'p';
    fname[3] = 'm';

    // Setting up Viewing Coordinat system 
    VectorType w = normalize_vector(view_dir);
    VectorType u = cross_product(w, up_dir);
    normalize_vector(u);
    VectorType v = cross_product(u,w);
    normalize_vector(v);
   

    float d =1;//= (tan((v_fov *( M_PI / 180.0))/ 2)) * ( pic_size.h ); 
    // d is DISTANCE FORM VIEWING WINDOW 
   

    float aspect_ratio = pic_size.h / pic_size.w; 
    // Aspect ratio for real image size
    
    float h_fov = 2 * atan(aspect_ratio * tan((v_fov * (M_PI / 180.0) )/ 2.0) );
    // Horizontal field of View

    h_fov = ( 180.0/ M_PI) *h_fov; 
    // Converting to degrees
    
    float window_height = 2 * d * tan(v_fov * (M_PI / 360.0));
    float window_width = 2 * d * tan(h_fov * (M_PI / 360.0));
    
   
    VectorType ul;
    // Vector for Defining upper left corner of viewing window 
    ul.x = eye.x + (w.x  * d) - (u.x * (window_width/2)) +  (v.x * (window_height/2));
    ul.y = eye.y +( w.y  * d) - (u.y * (window_width/2)) +  (v.y * (window_height/2));
    ul.z = eye.z +( w.z  * d) - (u.z * (window_width/2)) +  (v.z * (window_height/2));

    VectorType ur;
    // Vector for Defining upper right corner of viewing window 
    ur.x = eye.x + (w.x  * d) + (u.x *  window_width/2) +  (v.x * window_height/2);
    ur.y = eye.y +( w.y  * d) + (u.y *  window_width/2) +  (v.y * window_height/2);
    ur.z = eye.z +( w.z  * d) + (u.z *  window_width/2) +  (v.z * window_height/2);

    VectorType ll;
    // Vector for Defining lower left corner of viewing window 
    ll.x = eye.x + (w.x  * d) - (u.x *  window_width/2) -  (v.x * window_height/2);
    ll.y = eye.y +( w.y  * d) - (u.y *  window_width/2) -  (v.y * window_height/2);
    ll.z = eye.z +( w.z  * d) - (u.z *  window_width/2) -  (v.z * window_height/2);

    VectorType lr; 
    // Vector for Defining lower right corner of viewing window 
    lr.x = eye.x + (w.x  * d) + (u.x *  window_width/2) -  (v.x * window_height/2);
    lr.y = eye.y +( w.y  * d) + (u.y *  window_width/2) -  (v.y * window_height/2);
    lr.z = eye.z +( w.z  * d) + (u.z *  window_width/2) -  (v.z * window_height/2);


    VectorType h_offset;
    // Offset to move through viewing window in horizontal direction 
    h_offset.x = (ur.x - ul.x ) / (pic_size.w - 1);
    h_offset.y = (ur.y - ul.y ) / (pic_size.w - 1);
    h_offset.z = (ur.z - ul.z ) / (pic_size.w - 1);
    VectorType v_offset; 
    // Offset to move through viewing window in vertical direction  
    v_offset.x = (ll.x - ul.x) / (pic_size.h -1);
    v_offset.y = (ll.y - ul.y) / (pic_size.h -1);
    v_offset.z = (ll.z - ul.z) / (pic_size.h -1);

   
    printf("%s", file1); //check file name is correct
    FILE *fp2 = fopen(file2, "w");
     if(fp2 == NULL){
        printf("ERROR CANT OPEN FILE TO WRITE IN");
        return -1;

    }
    fputs("P3\n", fp2);
    fprintf(fp2, "%d %d\n", pic_size.w, pic_size.h); 
    fputs("255\n", fp2);
    
    for(int i = 0; i < pic_size.h; i++){
        for(int j = 0; j < pic_size.w; j++){
            VectorType curr_pixel;
            // Curr pixel to know which pixel the ray is pointing towards
            curr_pixel.x = ul.x + (j * h_offset.x) + (i * v_offset.x);
            curr_pixel.y = ul.y + (j * h_offset.y) + (i * v_offset.y);
            curr_pixel.z = ul.z + (j * h_offset.z) + (i * v_offset.z);
          
            ColorType pixel_color = back_ground_color; 
            

            VectorType ray_dir ;
            // Ray direction 
            ray_dir.x = curr_pixel.x - eye.x;
            ray_dir.y = curr_pixel.y - eye.y;
            ray_dir.z = curr_pixel.z - eye.z;
            ray_dir = normalize_vector(ray_dir); //ray direction (p-e)/||p-e||
           
            RayType ray = make_ray(eye,ray_dir);
            // Making Vector ray dir into a rayType

            //pixel_color = Trace_ray(ray,spheres, mtl_colors, back_ground_color);
            pixel_color = Trace_ray(ray,spheres, mtl_colors, back_ground_color, lights, light_count, triangles,vertices,verticeN, textureImages, coords,0);
            // Calls Trace Ray on current pixels color, if hits a sphere colors the sphere color 
            // Other wise returns the background color 
            
            fprintf(fp2, "%d %d %d\n", (int)(255.0 * fmin(fmax(pixel_color.r, 0.0),1.0)), 
            (int)(255.0 * fmin(fmax(pixel_color.g,0.0),1.0)), (int)(255.0 * fmin(fmax(pixel_color.b,0.0),1.0)));
            // Writing pixel colors to ppm file multiplies by 255 since current values in 0-1

      }
    }
   
    fclose(fp2);
    free(mtl_colors);
    free(lights);
    free(spheres);
    free(triangles);
    free(vertices);
    free(verticeN);
    free(coords);
    free_textures(textureImages);
    return 0;
 
}