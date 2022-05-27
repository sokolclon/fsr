#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <math.h>

#include <stdio.h>
#include <stdlib.h>

#define E 2.718281828459045

void rgba_to_grayscale(float *img_in,float *img_out, int iw, int ih, int n) {
	int i, j = 0;
	for (i = 0; i < ih * iw * 1; i++) {
		img_out[i] = (img_in[j] * 11 + img_in[j + 1] * 16 + img_in[j + 2] * 5) / 32;
		if (n == 4) {
			j += 4;
		}
		else if (n == 3) {
			j += 3;
		}
	}
}

void grayscale_to_rgba(float *img_in, float *img_out, int iw, int ih) {
	int i, j = 0;
	for (i = 0; i < ih * iw * 1; i++) {
		img_out[j] = img_in[i];
		img_out[j + 1] = img_in[i];
		img_out[j + 2] = img_in[i];
		img_out[j + 3] = 255;
		j += 4;
	}
}

float gaussian(float x, float mu, float sigma) {
	return pow(E,(-((x - mu) * (x - mu)) / (2 * sigma * sigma)));
}

float* gaussian_kernel(int size, float sigma, int normalize) {
	float *kernel = malloc(sizeof(float) * size);
	float sum;
	for (int i = 0; i < size; i++) {
		kernel[i] = gaussian(i, size / 2, sigma);
		sum += kernel[i];
	}
	if (normalize) {
		for (int i = 0; i < size; i++) {
			kernel[i] /= sum;
		}
	}
	return kernel;
}


void expand_array(float *array,float *expanded_array, int radius, int ih, int iw) {
	int nih = ih + 2 * radius;
	int niw = iw + 2 * radius;
	int i = 0;
	for (int k = 0; k < radius; k++) {
		for(int l = 0; l < radius; l++) {
			expanded_array[i] = array[0 * iw + 0];
			i++;
		}
		for (int l = 0; l < iw; l++) {
			expanded_array[i] = array[0 * iw + l];
			i++;
		}
		for (int l = 0; l < radius; l++) {
			expanded_array[i] = array[0 * iw + iw - 1];
			i++;
		}
	}
	for (int k = 0; k < ih; k++) {
		for (int l = 0; l < radius; l++) {
			expanded_array[i] = array[k * iw + 0];
			i++;
		}
		for (int l = 0; l < iw; l++) {
			expanded_array[i] = array[k * iw + l];
			i++;
		}
		for (int l = 0; l < radius; l++) {
			expanded_array[i] = array[k * iw + iw - 1];
			i++;
		}
	}
	for (int k = 0; k < radius; k++) {
		for (int l = 0; l < radius; l++) {
			expanded_array[i] = array[(ih - 1) * iw + 0];
			i++;
		}
		for (int l = 0; l < iw; l++) {
			expanded_array[i] = array[(ih - 1) * iw + l];
			i++;
		}
		for (int l = 0; l < radius; l++) {
			expanded_array[i] = array[(ih - 1) * iw + iw - 1];
			i++;
		}
	}
}


void mul_h(float *arr_i,float *arr_o, float *kernel, int radius, int index) {
	float pix = 0;
	for (int i = -radius; i <= radius ; i++) {
		pix += arr_i[index + i] * kernel[i + radius];
	}
	arr_o[index] = pix;
}

void mul_v(float *arr_i, float *arr_o, float *kernel, int radius, int index, int iw) {
	float pix = 0;
	for (int i = -radius; i <= radius ; i++) {
		pix += arr_i[index + i * iw] * kernel[i + radius];
	}
	arr_o[index] = pix;
}

void shrink_array(float *array, float *shrinked_array, int iw, int ih, int radius) {
	int i = radius * (2 * radius + iw) + radius;
	for (int j = 0; j < ih ; j++) {
		for (int k = 0; k < iw ; k++) {
			shrinked_array[j * iw + k] = array[i];
			i++;
		}
		i += 2 * radius;
	}

}

void gaussian_blur(float *img_in, float *img_out, int iw, int ih, int kernel_size, float sigma, int normalize) {
	float *kernel;
	int radius = (kernel_size - 1) / 2;
	int niw = iw + kernel_size - 1;
	int nih = ih + kernel_size - 1;
	kernel = gaussian_kernel(kernel_size, sigma, normalize);
	float *e_array = malloc((niw * nih) * sizeof(float));
	expand_array(img_in,e_array,radius,ih,iw);
	float *img_blur_h = malloc((niw * nih) * sizeof(float));
	float *img_blur_v = malloc((niw * nih) * sizeof(float));
	for(int i = 0; i < nih * niw; i++) {
		img_blur_h[i] = 0;
		img_blur_v[i] = 0;
	}
	for (int i = 0; i < nih; i++) {
		for (int j = 0; j < niw; j++) {
			if (j < radius || j > niw - radius) {
				img_blur_h[i * niw + j] = e_array[i * niw + j];
			}
			else {
				mul_h(e_array, img_blur_h, kernel, radius, i * niw + j);
			}
		}
	}
	for (int i = 0; i < nih; i++) {
		for (int j = 0; j < niw; j++) {
			if (i < radius || i > nih - radius) {
				img_blur_v[i * niw + j] = img_blur_h[i * niw + j];
			}
			else{
				mul_v(img_blur_h, img_blur_v, kernel, radius, i * niw + j, niw);
			}
		}
	}
	shrink_array(img_blur_v,img_out,iw,ih,radius);
}

void compressor(float *img_in,float *img_out, int iw, int ih) {
	float max = 0;
	float min = 1000;
	for (int i = 0; i < iw * ih; i++) {
		if (img_in[i] > max) {
			max=img_in[i];
		}
		if (img_in[i] < min) {
			min=img_in[i];
		}
	}
	for (int i = 0; i < iw * ih; i++) {
		img_out[i] = (img_in[i] - (max + min) / 2) / max * 127 + 127;
	}
}

unsigned char* arr_to_img(float *img, unsigned char *img_out, int iw, int ih, int n) {		
	for (int i = 0; i < iw * ih * n; i++) {
		img_out[i] = (unsigned char)img[i];
	}
	return img_out;
}

void copy_array(int *src, int *dst) {
	for(int i = 0; i < sizeof(src) / sizeof(src[0]); i++){
		dst[i] = src[i];
	}
}

int round_to_nearest(int n, int *targets, int size) {
	int min = INT_MAX;
	int closest = 0;
	for (int i = 0; i < size; i++) {
		if (abs(targets[i] - n) < min) {
			min = abs(targets[i] - n);
			closest = targets[i];
		}
	}
	return closest;
}

int ind(int *arr, float val, int size) {
	for (int i = 0; i < size; i++) {
		if (arr[i] == val) {
			return i;
		}
	}
	return -1;
}

void k_mean_for_colors(int *colors, int *color_belonging, int k, int iw, int ih) {
	int num_of_colors = 256;
	float l = 256 / k / 2;
	int *centers = (int*)malloc(k * sizeof(int));
	int *color_belonging_old = (int*)malloc(num_of_colors * sizeof(int));
	for(int i = 0; i < k; i++) {
		centers[i] = l * (l + 2 * i);
	}
	for(int i = 0; i < num_of_colors; i++) {
		color_belonging[i] = (int)round((i / (2 * l)));
		color_belonging_old[i] = (int)round((i / (2 * l)));
	}
	int  too_long = 0;
	while (1) {
		for (int affinity = 0; affinity < k; affinity++) {
			float M = 0;
			float c = 0;
			for(int i = 0; i < num_of_colors; i++) {
				if (color_belonging[i] == affinity) {
					M += colors[i];
					c += colors[i] * i;
				}
			}
			if (M == 0) {
				continue;
			}
			centers[affinity] = c / M;
		}
		for (int i = 0; i < num_of_colors; i++) {
			color_belonging[i] = ind(centers, round_to_nearest(i, centers, k), k);
		}
		if (color_belonging_old == color_belonging) {
			break;
		}
		copy_array(color_belonging, color_belonging_old);
		if (too_long > 100) {
			break;
		}
		too_long++;
	}
	free(centers);
	free(color_belonging_old);
}

void get_quantity_of_each_color(float *image, int *colors_out, int iw, int ih) {
	int n = 1;
	int num_of_colors = 256;
	for (int i = 0; i < num_of_colors; i++) {
		colors_out[i] = 0;
	}
	for(int i = 0; i < iw * ih * n; i++) {
		colors_out[(int)(image[i])] += 1;
	}
}

void paint_image(float *image_in, float *image_out, int iw, int ih, int *colors, int k) {
	int pallet [] = {
		250, 175, 196,
		205, 184, 145,
		153, 50, 204,
		178, 236, 93,
		255, 168, 175,
		0, 153, 0,
		255, 36, 0,
		153, 255, 153,
		163, 51, 255,
		237, 145, 33,
		162, 35, 29,
		216, 169, 3,
		120, 219, 226,
		0, 191, 255,
		140, 203, 94,
		106, 7, 254,
	};
	int j = 0;
	for (int i = 0; i < iw * ih; i++) {
		image_out[j++] = pallet[(int)(colors[(int)(image_in[i])]) * 3 + 0];
		image_out[j++] = pallet[(int)(colors[(int)(image_in[i])]) * 3 + 1];
		image_out[j++] = pallet[(int)(colors[(int)(image_in[i])]) * 3 + 2];
		image_out[j++] = 255;
	}
}

int main(int argc, char **argv) {
	char *inputPath;
	int number_of_areas, radius;
	float sigma;

	inputPath = argv[1];
	number_of_areas = atoi(argv[2]);
	radius = atoi(argv[3]);
	sigma = atof(argv[4]);
	
	int iw, ih, n;
	unsigned char *idata = stbi_load(inputPath, &iw, &ih, &n, 0);
	int m = 4;

	if(idata==NULL){
		printf("error loading image\n");
		return -1;
	}


	float* img = (float*)malloc((iw * ih * n) * sizeof(float));
	float* black_white = (float*)malloc((iw * ih * 1) * sizeof(float));
	float* img_blured = (float*)malloc((iw * ih * 1) * sizeof(float));
	float* img_compessed = (float*)malloc((iw * ih * 1) * sizeof(float));
	float* img_out_rgb = (float*)malloc((iw * ih * 4) * sizeof(float));
	unsigned char *img_out = (unsigned char*)malloc((iw * ih * m) * sizeof(char));

	int* colors_count = (int*)malloc((255) * sizeof(int));
	int* colors_by_clasters = (int*)malloc((255) * sizeof(int));
	float* colored_image = (float*)malloc((iw * ih * 4) * sizeof(float));

	for (int i = 0; i < iw * ih * n; i++) {
		img[i] = idata[i];
	}

	rgba_to_grayscale(img, black_white, iw, ih,n);
	gaussian_blur(black_white, img_blured, iw, ih,5,1.5,1);
	compressor(img_blured, img_compessed, iw, ih);

	get_quantity_of_each_color(img_compessed,colors_count,iw,ih);
	k_mean_for_colors(colors_count,colors_by_clasters,number_of_areas,iw,ih);
	paint_image(img_blured,colored_image,iw,ih,colors_by_clasters,number_of_areas);

	arr_to_img(colored_image,img_out,iw,ih,m);

	char *outputPath = "result_skull.png";
	stbi_write_png(outputPath, iw, ih, m, img_out, 0);
	
	stbi_image_free(img);
	stbi_image_free(idata);
	stbi_image_free(black_white);
	stbi_image_free(img_blured);
	stbi_image_free(img_out_rgb);
	stbi_image_free(img_out);

	return 0;
}
