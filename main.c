#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <math.h>

#include <strings.h>
#include <stdio.h>
#include <stdlib.h>

#define E 2.718281828459045

void rgba_to_grayscale(unsigned char *img_in, float *black_white, int iw, int ih, int n) {
	int j = 0;
	for (int i = 0; i < ih * iw * 1; i++) {
		black_white[i] = (img_in[j] * 11 + img_in[j + 1] * 16 + img_in[j + 2] * 5) / 32;
		if (n == 4) {
			j += 4;
		}
		else if (n == 3) {
			j += 3;
		}
	}
}

void grayscale_to_rgba(float *img_in, float *img_out, int iw, int ih)
{
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

void gaussian_kernel(float * kernel, int size, float sigma, int normalize) {
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
	return;
}

float* expand_array(float *array, int size, int ih, int iw) {
	int i, j;
	int nih = ih + 2 * size;
	int niw = iw + 2 * size;
	float *expanded_array = (float*)malloc((niw * nih) * sizeof(float));
	for (int k = 0; k < nih * niw; k++) {
		i = k % niw;
		j = k / niw;
		if (i < size){
			i = size;
		}
		if (i > niw - size){
			i = niw - size;
		}
		if (j < size){
			j = size;
		}
		if (j > nih - size){
			j = nih - size;
		}
		expanded_array[k] = array[j * ih + i];
	}
	return expanded_array;
}


void mul_h(float *arr_i, float *arr_o, float *kernel, int radius, int index) {
	float pix = 0;
	for (int i = -radius; i <= radius; i++) {
		pix += arr_i[index + i] * kernel[i + radius];
	}
	arr_o[index] = pix;
}

void mul_v(float *arr_i, float *arr_o, float *kernel, int radius, int index, int iw){
	float pix = 0;
	for (int i = -radius; i <= radius; i++) {
		pix += arr_i[index + i * iw] * kernel[i + radius];
	}
	arr_o[index] = pix;
}

void gaussian_blur(float *img_in, float *img_out, int iw, int ih, int kernel_size, float sigma, int normalize) {
	float* kernel = (float*)malloc((kernel_size) * sizeof(float));
	gaussian_kernel(kernel, kernel_size, sigma, normalize);

	int radius = (kernel_size - 1) / 2, niw = iw + kernel_size - 1, nih = ih + kernel_size - 1;
	float* e_array = (float*)malloc((niw * nih) * sizeof(float));
	e_array = expand_array(img_in, radius, ih, iw);

	float *img_blur_h = malloc((nih * niw) * sizeof(float));
	float *img_blur_v = malloc((nih * niw) * sizeof(float));

	for (int i = 0; i < niw; i++) {
		for (int j = 0; j < nih; j++) {
			if (j <= radius || j >= niw - radius) {
				img_blur_h[i * niw + j] = e_array[i * niw + j];
			}
			else {
				mul_h(e_array, img_blur_h, kernel, radius, i * niw + j);
			}
		}
	}
	/*
	for (int i = 0; i < niw; i++) {
		for (int j = 0; j < nih; j++) {
			if(i <= radius || i >= niw - radius){
				img_blur_v[i * niw + j] = img_blur_h[i * niw + j];
			}
			else{
				mul_v(img_blur_h, img_blur_v, kernel, radius, i * niw + j, niw);
			}
		}
	}
	*/
	img_blur_v = img_blur_h;
	for (int i = 0; i < iw; i++) {
		for (int j = 0; j < ih; j++) {
			img_out[i * iw + j] = img_blur_h[(i + radius) * niw + (j + radius)];
		}
	}
}

unsigned char* arr_to_img(int iw, int ih, float* img) {
	unsigned char *img_out = (unsigned char *)malloc((iw * ih * 4) * sizeof(char));
	for (int i = 0; i < iw * ih * 4; i++){
		img_out[i] = img[i];
	}
	return img_out;
}

int main() {
	char *inputPath = "hamster.png";
	//char *inputPath = "3.png";
	int iw, ih, n;
	unsigned char *idata = stbi_load(inputPath, &iw, &ih, &n, 0);
	if (idata == NULL) {
		printf("ERROR: can't read file %s\n", inputPath);
		return 1;
	}

	float* img = (float *)malloc((iw * ih * n) * sizeof(float));
	float* black_white = (float *)malloc((iw * ih * 1) * sizeof(float));
	float* img_blured = (float *)malloc((iw * ih * 1) * sizeof(float));
	float* img_out_rgb = (float *)malloc((iw * ih * 4) * sizeof(float));
	unsigned char *img_out = (unsigned char *)malloc((iw * ih * 4) * sizeof(char));

	for (int i = 0; i < iw * ih * n; i++) {
		img[i] = (float)idata[i];
	}

	rgba_to_grayscale(idata, black_white, iw, ih, n);
//	gaussian_blur(black_white, img_blured, iw, ih, 3, 1, 1);
//	grayscale_to_rgba(img_blured, img_out_rgb, iw, ih);

//	img_out = arr_to_img(iw, ih, img_blured);

	char *outputPath = "result.png";
	stbi_write_png(outputPath, iw, ih, 1, black_white, 0);

	stbi_image_free(img);
	stbi_image_free(idata);
	stbi_image_free(black_white);
	stbi_image_free(img_blured);
	stbi_image_free(img_out_rgb);
	stbi_image_free(img_out);

	return 0;
}
