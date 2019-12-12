/*
	Raytracer
	- Avishkar Kolahalu
*/

#include <windows.h>
#include <stdlib.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <stdio.h>
#include <string>

#include "opencv2/core/core.hpp"
#include "opencv2/core/utility.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/imgcodecs/imgcodecs.hpp"

using namespace std;

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename = 0;

// Different Display Modes:
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode = MODE_DISPLAY;

// Window Size:
#define WIDTH 320
#define HEIGHT 240

// Anti-aliasing by supersampling - Change to test - (Higher value takes longer) - Change to 1 to remove.
#define SAMPLER_VALUE 3

double sampled_height = HEIGHT * SAMPLER_VALUE;	// Height bound.
double sampled_width = WIDTH * SAMPLER_VALUE;	// Width bound.

#define fov 60.0	// Field of view of camera.

unsigned char buffer[HEIGHT][WIDTH][3];

double screen_space_width = 0.0, screen_space_height = 0.0;

struct Vertex
{
	double position[3];
	double color_diffuse[3];
	double color_specular[3];
	double normal[3];
	double shininess;
};

// For x, y, z values
struct point
{
	double x, y, z;
};

// For r, g, b values
struct colour
{
	double r, g, b;
};

typedef struct _Triangle
{
	struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
	double position[3];
	double color_diffuse[3];
	double color_specular[3];
	double shininess;
	double radius;
} Sphere;

typedef struct _Light
{
	double position[3];
	double color[3];
} Light;

// Ray structure
struct ray_t
{
	point d;	// direction.
	point o;	// origin.
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];

double ambient_light[3];

int num_triangles = 0, num_spheres = 0, num_lights = 0;

// Corners of the screen space:

point left_top;		// Left-top of screen
point left_bottom;	// Left-bottom of screen

point right_top;	// Right-top of screen
point right_bottom;	// Right-bottom of screen

point cam_point = { 0.0, 0.0, 0.0 };

colour** pixels;	// Screen pixel matrix that stores rgb values

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
	glColor3f(((double)r) / 256.f, ((double)g) / 256.f, ((double)b) / 256.f);
	glVertex2i(x, y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
	buffer[HEIGHT - y - 1][x][0] = r;
	buffer[HEIGHT - y - 1][x][1] = g;
	buffer[HEIGHT - y - 1][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
	plot_pixel_display(x, y, r, g, b);
	if (mode == MODE_JPEG)
		plot_pixel_jpeg(x, y, r, g, b);
}

// Calculates r(t) = o + td
point rayPosition(ray_t shot_ray, double t)
{
	point ray;

	ray.x = shot_ray.o.x + t * shot_ray.d.x;
	ray.y = shot_ray.o.y + t * shot_ray.d.y;
	ray.z = shot_ray.o.z + t * shot_ray.d.z;

	return ray;
}

// Dot product
double dotProduct(point t1, point t2)
{
	return (t1.x * t2.x) + (t1.y * t2.y) + (t1.z * t2.z);
}

// Cross product
point xProduct(point t1, point t2)
{
	point temp;

	temp.x = (t1.y * t2.z) - (t1.z * t2.y);
	temp.y = (t1.z * t2.x) - (t1.x * t2.z);
	temp.z = (t1.x * t2.y) - (t1.y * t2.x);

	return temp;
}

// Distance between 2 points
double calcDist(point t1, point t2)
{
	double x_d, y_d, z_d;

	x_d = t1.x - t2.x;
	y_d = t1.y - t2.y;
	z_d = t1.z - t2.z;

	return sqrt((x_d * x_d) + (y_d * y_d) + (z_d * z_d));
}

// Normaliser
point normalise(point t1)	
{
	double n_d = sqrt((t1.x * t1.x) + (t1.y * t1.y) + (t1.z * t1.z));

	t1.x = t1.x / n_d;
	t1.y = t1.y / n_d;
	t1.z = t1.z / n_d;

	return t1;
}

// Vector Subtraction
point coordSubtract(point t1, point t2)
{
	point temp;

	temp.x = t1.x - t2.x;
	temp.y = t1.y - t2.y;
	temp.z = t1.z - t2.z;

	return temp;
}

// Find Reflection, given the normal
point Reflector(point incident, point normal)
{
	double cos = dotProduct(incident, normal);

	point reflection;

	reflection.x = (cos * normal.x * 2) - incident.x;
	reflection.y = (cos * normal.y * 2) - incident.y;
	reflection.z = (cos * normal.z * 2) - incident.z;

	return reflection;
}

// Intersection check - Ray x Light source
bool xLightSrc(int light_no, ray_t shot_ray, double &x_dist)
{
	double y_dist, z_dist;

	x_dist = (lights[light_no].position[0] - shot_ray.o.x) / shot_ray.d.x;
	y_dist = (lights[light_no].position[1] - shot_ray.o.y) / shot_ray.d.y;
	z_dist = (lights[light_no].position[2] - shot_ray.o.z) / shot_ray.d.z;

	if ((x_dist != y_dist) || (x_dist != z_dist))
		return false;
	else
		return true;
}

// Intersection check - Ray x Sphere
bool xSphere(int i_sph, ray_t shot_ray, point &N, double &dist)
{
	double b, c, calc, t0, t1, x_diff, y_diff, z_diff, s_x, s_y, s_z;
	point temp;

	// Calculate Solution: 

	s_x = spheres[i_sph].position[0];
	s_y = spheres[i_sph].position[1];
	s_z = spheres[i_sph].position[2];

	x_diff = shot_ray.o.x - s_x;
	y_diff = shot_ray.o.y - s_y;
	z_diff = shot_ray.o.z - s_z;

	// Note: a = 1

	b = 2.0 * (shot_ray.d.x * (shot_ray.o.x - s_x) + shot_ray.d.y * (shot_ray.o.y - s_y) + shot_ray.d.z * (shot_ray.o.z - s_z));

	c = (x_diff * x_diff) + (y_diff * y_diff) + (z_diff * z_diff) - (spheres[i_sph].radius * spheres[i_sph].radius);

	calc = (b * b) - (4.0 * c);

	if (calc >= 0.0)
	{
		t0 = (-b + sqrt(calc)) / 2;
		t1 = (-b - sqrt(calc)) / 2;
	}

	else
		return false;

	// Ensure both aren't less than zero
	if (t1 <= 0 && t0 <= 0)
		return false;

	// Find minimal of the solutions:
	(t0 > 0 && t1 > 0) ? dist = MIN(t0, t1) : dist = MAX(t0, t1);

	// limiter
	if (dist <= 0.000001)
		return false;	

	temp = rayPosition(shot_ray, dist);

	// Normalise:

	N.x = temp.x - s_x;
	N.y = temp.y - s_y;
	N.z = temp.z - s_z;

	N = normalise(N);

	return true;
}

// Intersection check - Ray and Triangle
bool xTriangle(int i_tri, ray_t shot_ray, point &N, double &alpha, double &beta, double &gamma, double &t)
{
	point p, p1, p2, p3, temp;	// The points on the triangle
	double n_dot_d, area;

	p1.x = triangles[i_tri].v[0].position[0];
	p2.x = triangles[i_tri].v[1].position[0];
	p3.x = triangles[i_tri].v[2].position[0];

	p1.y = triangles[i_tri].v[0].position[1];
	p2.y = triangles[i_tri].v[1].position[1];
	p3.y = triangles[i_tri].v[2].position[1];

	p1.z = triangles[i_tri].v[0].position[2];
	p2.z = triangles[i_tri].v[1].position[2];
	p3.z = triangles[i_tri].v[2].position[2];

	temp = xProduct(coordSubtract(p2, p1), coordSubtract(p3, p1));

	N = normalise(temp);	// Normal

	n_dot_d = dotProduct(N, shot_ray.d);

	// n.d = 0; therefore, no intersection
	if (n_dot_d == 0.0)
		return false;

	t = -(dotProduct(N, coordSubtract(shot_ray.o, p1))) / n_dot_d;	//	-(n.p1 + d)/n.d

	// Intersection is behind ray origin
	if (round(t) <= 0)
		return false;

	// Intersection Test:

	p = rayPosition(shot_ray, t);

	temp = xProduct(coordSubtract(p3, p2), coordSubtract(p, p2));
	alpha = dotProduct(N, temp);		// Area(p2p3p)

	temp = xProduct(coordSubtract(p1, p3), coordSubtract(p, p3));
	beta = dotProduct(N, temp);			// Area(p3p1p)

	temp = xProduct(coordSubtract(p2, p1), coordSubtract(p, p1));
	gamma = dotProduct(N, temp);		// Area(p1p2p)

	// Intersection within Triangle
	if (alpha >= 0 && beta >= 0 && gamma >= 0)
	{
		area = dotProduct(N, xProduct(coordSubtract(p2, p1), coordSubtract(p3, p1)));	// Triangle Area

		alpha = alpha / area;
		beta = beta / area;
		gamma = gamma / area;

		return true;
	}

	return false;
}

// Phong Model function
void phongShader(point int_point, point N, colour &Colour, point kd, point ks, double shiny, point viewer)
{
	int i_light, i_sph, i_tri;	// Indices

	bool shadow_checker;

	point light_xyz, int_dir, t1_normal, t2_normal;
	ray_t shade;		// Shadow Ray

	double light_dist, sph_dist;
	double tri_dist, a, b, c;
	double cos_theta, cos_phi, Lr, Lg, Lb;

	// Send shadow ray to every light
	for (i_light = 0; i_light < num_lights; i_light++)
	{
		shadow_checker = false;

		light_xyz.x = lights[i_light].position[0];
		light_xyz.y = lights[i_light].position[1];
		light_xyz.z = lights[i_light].position[2];

		light_dist = calcDist(light_xyz, int_point);	// Distance b/w origin point and light source

		int_dir.x = lights[i_light].position[0] - int_point.x;
		int_dir.y = lights[i_light].position[1] - int_point.y;
		int_dir.z = lights[i_light].position[2] - int_point.z;

		int_dir = normalise(int_dir);

		// Ray shot towards light source
		shade.o = int_point;
		shade.d = int_dir;

		// Check Shadow Ray Intersection with Spheres.
		for (i_sph = 0; i_sph < num_spheres; i_sph++)
			if (xSphere(i_sph, shade, t1_normal, sph_dist))
				if (calcDist(rayPosition(shade, sph_dist), int_point) <= light_dist)
					shadow_checker = true;		// Reached sphere before reaching a light source. So, it's in shadow

		// Check Shadow Ray Intersection with Triangles.
		for (i_tri = 0; i_tri < num_triangles; i_tri++)
			if (xTriangle(i_tri, shade, t2_normal, a, b, c, tri_dist))
				if (calcDist(rayPosition(shade, tri_dist), int_point) <= light_dist)
					shadow_checker = true;		// Reached triangle before reaching a light source. So, it's in shadow.

		if (shadow_checker == false)
		{
			// Colour calcuation:

			cos_theta = dotProduct(int_dir, N);
			cos_phi = dotProduct(normalise(Reflector(int_dir, N)), viewer);

			// Clamp: 
			if (cos_theta < 0.0)
				cos_theta = 0.0;

			if (cos_phi < 0.0)
				cos_phi = 0.0;

			// Light contribution:

			Lr = lights[i_light].color[0];
			Lg = lights[i_light].color[1];
			Lb = lights[i_light].color[2];

			// Phong model equation:

			Colour.r += Lr * (kd.x*cos_theta + ks.x*pow(cos_phi, shiny));

			Colour.g += Lg * (kd.y*cos_theta + ks.y*pow(cos_phi, shiny));

			Colour.b += Lb * (kd.z*cos_theta + ks.z*pow(cos_phi, shiny));
		}
	}
}

// Raytracer Function
void raytracer(ray_t shot_ray, colour &RGB)
{
	point N_S, N_T;					// Normals
	point x_S, x_T;					// Intersection points
	point kd_S, kd_T, ks_S, ks_T;	// Coefficients
	point v_S, v_T;					// To Viewer

	colour Colour_S, Colour_T;

	int i;

	bool x_checker = false;	// Intersection checker

	double light_dist, sph_dist, tri_dist, shiny, shiny_T, a, b, g;

	double min_obj_dist = 1000000000.0;	// Tracks closer objects

	// Triangle Intersection check with every triangle
	for (int index = 0; index < num_triangles; index++)
		if (xTriangle(index, shot_ray, N_T, a, b, g, tri_dist))
			if (min_obj_dist > tri_dist)
			{
				min_obj_dist = tri_dist;

				// Viewer:

				v_T.x = -1.0 * shot_ray.d.x;
				v_T.y = -1.0 * shot_ray.d.y;
				v_T.z = -1.0 * shot_ray.d.z;

				v_T = normalise(v_T);

				// Normal:

				N_T.x = triangles[index].v[0].normal[0] * a;
				N_T.x += triangles[index].v[1].normal[0] * b;
				N_T.x += triangles[index].v[2].normal[0] * g;

				N_T.y = triangles[index].v[0].normal[1] * a;
				N_T.y += triangles[index].v[1].normal[1] * b;
				N_T.y += triangles[index].v[2].normal[1] * g;

				N_T.z = triangles[index].v[0].normal[2] * a;
				N_T.z += triangles[index].v[1].normal[2] * b;
				N_T.z += triangles[index].v[2].normal[2] * g;

				// Diffuse coefficient:

				kd_T.x = triangles[index].v[0].color_diffuse[0] * a;
				kd_T.x += triangles[index].v[1].color_diffuse[0] * b;
				kd_T.x += triangles[index].v[2].color_diffuse[0] * g;

				kd_T.y = triangles[index].v[0].color_diffuse[1] * a;
				kd_T.y += triangles[index].v[1].color_diffuse[1] * b;
				kd_T.y += triangles[index].v[2].color_diffuse[1] * g;

				kd_T.z = triangles[index].v[0].color_diffuse[2] * a;
				kd_T.z += triangles[index].v[1].color_diffuse[2] * b;
				kd_T.z += triangles[index].v[2].color_diffuse[2] * g;


				// Specular coefficient:

				ks_T.x = triangles[index].v[0].color_specular[0] * a;
				ks_T.x += triangles[index].v[1].color_specular[0] * b;
				ks_T.x += triangles[index].v[2].color_specular[0] * g;

				ks_T.y = triangles[index].v[0].color_specular[1] * a;
				ks_T.y += triangles[index].v[1].color_specular[1] * b;
				ks_T.y += triangles[index].v[2].color_specular[1] * g;


				ks_T.z = triangles[index].v[0].color_specular[2] * a;
				ks_T.z += triangles[index].v[1].color_specular[2] * b;
				ks_T.z += triangles[index].v[2].color_specular[2] * g;


				// Shininess Factor:

				shiny_T = triangles[index].v[0].shininess * a;
				shiny_T += triangles[index].v[0].shininess * b;
				shiny_T += triangles[index].v[0].shininess * g;


				Colour_T = { 0.0,0.0,0.0 };
				x_T = rayPosition(shot_ray, tri_dist);	// Intersection point

				phongShader(x_T, N_T, Colour_T, kd_T, ks_T, shiny_T, v_T);	// Calc. the colour

				RGB.r = 255.0 * Colour_T.r;
				RGB.g = 255.0 * Colour_T.g;
				RGB.b = 255.0 * Colour_T.b;

				x_checker = true;
			}

	// Sphere Intersection check with every sphere
	for (i = 0; i < num_spheres; i++)
		if (xSphere(i, shot_ray, N_S, sph_dist))
			if (min_obj_dist > sph_dist)
			{
				min_obj_dist = sph_dist;

				// Viewer:

				v_S.x = -1.0 * shot_ray.d.x;
				v_S.y = -1.0 * shot_ray.d.y;
				v_S.z = -1.0 * shot_ray.d.z;

				v_S = normalise(v_S);

				// Diffuse coefficient:

				kd_S.x = spheres[i].color_diffuse[0];
				kd_S.y = spheres[i].color_diffuse[1];
				kd_S.z = spheres[i].color_diffuse[2];

				// Specular coefficient:

				ks_S.x = spheres[i].color_specular[0];
				ks_S.y = spheres[i].color_specular[1];
				ks_S.z = spheres[i].color_specular[2];

				// Shininess Factor:

				shiny = spheres[i].shininess;


				Colour_S = { 0.0, 0.0, 0.0 };
				x_S = rayPosition(shot_ray, sph_dist);	// Intersection point

				phongShader(x_S, N_S, Colour_S, kd_S, ks_S, shiny, v_S);	// Calc. the colour

				RGB.r = 255.0 * Colour_S.r;
				RGB.g = 255.0 * Colour_S.g;
				RGB.b = 255.0 * Colour_S.b;

				x_checker = true;
			}

	// Light Source Intersection check with every light source
	for (i = 0; i < num_lights; i++)
		if (xLightSrc(i, shot_ray, light_dist))
			if (min_obj_dist > light_dist)
			{
				min_obj_dist = light_dist;

				// Light src colour:

				RGB.r = 255.0 * lights[i].color[0];
				RGB.g = 255.0 * lights[i].color[1];
				RGB.b = 255.0 * lights[i].color[2];

				x_checker = true;
			}

	if (x_checker == true)
	{
		// Add Ambient:

		RGB.r += 255.0 * ambient_light[0];
		RGB.g += 255.0 * ambient_light[1];
		RGB.b += 255.0 * ambient_light[2];
	}

	else
		RGB = { 255.0, 255.0, 255.0 };
}

void draw_scene()
{
	double width_jump;	// To skip to the next column
	double height_jump;	// To skip to the next row

	double x, y;
	int row, col;

	point pixel_xyz;	// Pixel co-ordinates
	ray_t shot_ray;		// Shot ray
	colour RGB;

	pixel_xyz.z = -1.0;	// z is always -1.0

	width_jump = screen_space_width / (sampled_width);
	height_jump = screen_space_height / (sampled_height);

	// Raytracing Loop:

	for (row = 0, y = left_bottom.y; row < sampled_height; row++, y += height_jump)
		for (col = 0, x = left_bottom.x; col < sampled_width; col++, x += width_jump)
		{
			RGB = { 0.0, 0.0, 0.0 };		// Init r, g, b values

			// Iterative pixel traversal to shoot rays
			pixel_xyz.x = x;
			pixel_xyz.y = y;

			// Ray shot from camera origin
			shot_ray.o = cam_point;

			// Direction of ray to be shot through pixel position
			shot_ray.d = normalise(coordSubtract(pixel_xyz, cam_point));

			// Trace the shot ray
			raytracer(shot_ray, RGB);

			// Clamp:

			RGB.r = MIN(255.0, RGB.r);
			RGB.r = MAX(0.0, RGB.r);

			RGB.g = MIN(255.0, RGB.g);
			RGB.g = MAX(0.0, RGB.g);

			RGB.b = MIN(255.0, RGB.b);
			RGB.b = MAX(0.0, RGB.b);

			// Assign calculated RGB values to pixel matrix:

			pixels[row][col].r = RGB.r;
			pixels[row][col].g = RGB.g;
			pixels[row][col].b = RGB.b;
		}

	// Render Loop:

	int c_width, c_height, w, h;
	double red, green, blue;
	double avger = SAMPLER_VALUE * SAMPLER_VALUE;	// Averager

	for (c_width = 0, col = 0; c_width < sampled_width; c_width += SAMPLER_VALUE, col++)
	{
		glBegin(GL_POINTS);
		for (c_height = 0, row = 0; c_height < sampled_height; c_height += SAMPLER_VALUE, row++)
		{
			// Supersampling:

			// Summing up r, g, b values of the 'subpixels':
			for (w = 0, red = 0.0, green = 0.0, blue = 0.0; w < SAMPLER_VALUE; w++)
				for (h = 0; h < SAMPLER_VALUE; h++)
				{
					red += pixels[c_height + h][c_width + w].r;
					green += pixels[c_height + h][c_width + w].g;
					blue += pixels[c_height + h][c_width + w].b;
				}

			// Average them, and render:
			plot_pixel(col, row, red / avger, green / avger, blue / avger);
		}
		glEnd();
	}

	glFlush();

	printf("Done!\n");
}

// Write a jpg image from buffer
void save_jpg()
{
	if (filename == NULL)
		return;

	// Allocate a picture buffer
	cv::Mat3b bufferBGR = cv::Mat::zeros(HEIGHT, WIDTH, CV_8UC3); // rows, cols, 3-channel 8-bit
	printf("File to save to: %s\n", filename);

	// unsigned char buffer[HEIGHT][WIDTH][3];
	for (int r = 0; r < HEIGHT; r++)
		for (int c = 0; c < WIDTH; c++)
			for (int chan = 0; chan < 3; chan++)
			{
				unsigned char red = buffer[r][c][0];
				unsigned char green = buffer[r][c][1];
				unsigned char blue = buffer[r][c][2];
				bufferBGR.at<cv::Vec3b>(r, c) = cv::Vec3b(blue, green, red);
			}

	if (cv::imwrite(filename, bufferBGR)) 
		printf("File saved Successfully\n");

	else 
		printf("Error in Saving\n");

}

void parse_check(char *expected, char *found)
{
	if (_stricmp(expected, found))
	{
		printf("Expected '%s ' found '%s '\n", expected, found);
		printf("Parse error, abnormal abortion\n");
		exit(0);
	}
}

void parse_doubles(FILE *file, char *check, double p[3])
{
	char str[100];
	fscanf(file, "%s", str);
	parse_check(check, str);
	fscanf(file, "%lf %lf %lf", &p[0], &p[1], &p[2]);
	printf("%s %lf %lf %lf\n", check, p[0], p[1], p[2]);
}

void parse_rad(FILE *file, double *r)
{
	char str[100];
	fscanf(file, "%s", str);
	parse_check((char *)"rad:", str);
	fscanf(file, "%lf", r);
	printf("rad: %f\n", *r);
}

void parse_shi(FILE*file, double *shi)
{
	char s[100];
	fscanf(file, "%s", s);
	parse_check((char *)"shi:", s);
	fscanf(file, "%lf", shi);
	printf("shi: %f\n", *shi);
}

// Load scene from cmd line specified text file
int loadScene(char *argv)
{
	FILE *file;
	fopen_s(&file, argv, "r");
	int number_of_objects;
	char type[50];
	int i;
	Triangle t;
	Sphere s;
	Light l;
	fscanf(file, "%i", &number_of_objects);

	printf("number of objects: %i\n", number_of_objects);
	char str[200];

	parse_doubles(file, (char *)"amb:", ambient_light);

	for (i = 0; i < number_of_objects; i++)
	{
		fscanf(file, "%s\n", type);
		printf("%s\n", type);
		if (_stricmp(type, "triangle") == 0)
		{

			printf("found triangle\n");
			int j;

			for (j = 0; j < 3; j++)
			{
				parse_doubles(file, (char *)"pos:", t.v[j].position);
				parse_doubles(file, (char *)"nor:", t.v[j].normal);
				parse_doubles(file, (char *)"dif:", t.v[j].color_diffuse);
				parse_doubles(file, (char *)"spe:", t.v[j].color_specular);
				parse_shi(file, &t.v[j].shininess);
			}

			if (num_triangles == MAX_TRIANGLES)
			{
				printf("too many triangles, you should increase MAX_TRIANGLES!\n");
				exit(0);
			}
			triangles[num_triangles++] = t;
		}
		else if (_stricmp(type, "sphere") == 0)
		{
			printf("found sphere\n");

			parse_doubles(file, (char *)"pos:", s.position);
			parse_rad(file, &s.radius);
			parse_doubles(file, (char *)"dif:", s.color_diffuse);
			parse_doubles(file, (char *)"spe:", s.color_specular);
			parse_shi(file, &s.shininess);

			if (num_spheres == MAX_SPHERES)
			{
				printf("too many spheres, you should increase MAX_SPHERES!\n");
				exit(0);
			}
			spheres[num_spheres++] = s;
		}
		else if (_stricmp(type, "light") == 0)
		{
			printf("found light\n");
			parse_doubles(file, (char *)"pos:", l.position);
			parse_doubles(file, (char *)"col:", l.color);

			if (num_lights == MAX_LIGHTS)
			{
				printf("too many lights, you should increase MAX_LIGHTS!\n");
				exit(0);
			}
			lights[num_lights++] = l;
		}
		else
		{
			printf("unknown type in scene description:\n%s\n", type);
			exit(0);
		}
	}
	return 0;
}

void display() { }

void init()
{
	double a, x, y;
	int t;

	glMatrixMode(GL_PROJECTION);

	glOrtho(0, WIDTH, 0, HEIGHT, 1, -1);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glClearColor(1.0, 1.0, 1.0, 0);	// White Background as suggested
	glClear(GL_COLOR_BUFFER_BIT);

	// Find the four corners:

	a = (double)WIDTH / (double)HEIGHT;
	x = a * tan(fov * 3.141592 / (2 * 180));
	y = tan(fov * 3.141592 / (2 * 180));

	screen_space_width = x * 2.0;
	screen_space_height = y * 2.0;

	left_top.x = -x;
	left_top.y = y;
	left_top.z = -1.0;

	left_bottom.x = -x;
	left_bottom.y = -y;
	left_bottom.z = -1.0;

	right_top.x = x;
	right_top.y = y;
	right_top.z = -1.0;

	right_bottom.x = x;
	right_bottom.y = -y;
	right_bottom.z = -1.0;

	// Initialise pixel matrix:

	pixels = new colour*[sampled_height];

	for (t = 0; t < sampled_height; t++)
		pixels[t] = new colour[sampled_width];
}

void idle()
{
	static int once = 0;
	if (!once)
	{
		draw_scene();
		if (mode == MODE_JPEG)
			save_jpg();
	}
	once = 1;
}

int main(int argc, char ** argv)
{
	if (argc < 2 || argc > 3)
	{
		printf("usage: %s <scenefile> [jpegname]\n", argv[0]);
		exit(0);
	}
	if (argc == 3)
	{
		mode = MODE_JPEG;
		filename = argv[2];
	}
	else if (argc == 2)
		mode = MODE_DISPLAY;

	glutInit(&argc, argv);

	loadScene(argv[1]);

	glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(WIDTH, HEIGHT);

	glutCreateWindow("Assignment 3: Ray Tracing - Avishkar Kolahalu");

	glutDisplayFunc(display);
	glutIdleFunc(idle);

	init();

	glutMainLoop();
}