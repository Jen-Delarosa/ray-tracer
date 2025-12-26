# ray-tracer

## Project Overview:
	
Implemented a Ray tracer program in C that recursively casts rays into a 3D scene to determine the color stored for each pixel. The program outputs a PPM image that displays the scene the user created. The program creates a scene from an input file that specifies the objects, materials, and lighting. 

The ray tracer supports input files containing multiple objects, multiple light sources, and different material properties, including texture mapping, smooth shading, flat shading, reflectance, and transparency. 

## Main Features: 

- File handling and scene parsing: 
The program parses input files to render the scene and image files to perform texture mapping. The program also writes to a PPM file, storing the RGB values for each pixel in the scene.
 
- Object creation and geometry:
Defines and manages objects used to create the scene, such as spheres, triangles, and vectors.

- Shadow detection:
The program handles multiple objects by detecting shadows to determine which object appears closer to the camera first.

- Lighting:
Used the Phong Illumination Model to calculate the diffuse, specular, and ambient terms, and determined the lighting of the scene.

- Reflectance:
Schlick's approximation is used to reflect rays off the object’s surfaces to determine which objects the ray hits and are reflected onto the current object in the scene.

- Transparency: 
Handles transparency by casting refracted rays through objects and blending the transmitted color with the surface color, so that objects behind transparent objects are still visible.

- Material types: 
Different material types, including reflective, transparent, solid, smooth shading, and flat shading. 

To run the program, use the following commands:
		1. Make ./raytracer 
		2. ./raytracer input.txt


## Directory structure: 

- The images/ folder contains images converted to PNG from the PPM files produced by the ray tracer.

- The src/ folder contains the source code, Makefile, and an example input file (input.txt).

- The scenes/ folder contains different input files that produce the images in the images/ folder. 

- The scenes/output/ folder contains the original PPM files produced by the ray tracer.

