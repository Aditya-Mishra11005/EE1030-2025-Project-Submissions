#ifndef IMAGE_H
#define IMAGE_H

int img_read_rgb(const char *path, int *h, int *w, double **R, double **G, double **B);

int jpg_read_rgb(const char *path, int *h, int *w, double **R, double **G, double **B);
int png_read_rgb(const char *path, int *h, int *w, double **R, double **G, double **B);

int jpg_write_rgb(const char *path, int h, int w, int quality, double *R, double *G, double *B);
int png_write_rgb(const char *path, int h, int w, double *R, double *G, double *B);

int ppm_write(const char *path, int h, int w, int mx, double *R, double *G, double *B);

#endif