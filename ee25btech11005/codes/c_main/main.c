#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "../c_libs/image.h"
#include "../c_libs/svd.h"
int gray(double *R, double *G, double *B, size_t pix, int m, int n) {
    int err = 0;
    for (int i = 0; i < m; i++) {
        if (R[i] != G[i] || R[i] != B[i]) {
            err++;
        }
    }
    if (err < 0.01 * m * n) return 1;
    return 0;
}
double err(double *A, double *Ak, int m, int n) {
    double num = 0.0, den = 0.0;
    int size = m * n;
    for (int i = 0; i < size; i++) {
        double diff = A[i] - Ak[i];
        num += diff * diff;
        den += A[i] * A[i];
    }
    if (den == 0.0) return 0.0;
    return sqrt(num / den);
}


int main() {
    int m, n, max = 255;
    double *Rin, *Gin, *Bin;
    char ifname[100], ofname[100];
    int k;
    int q = 90;
    printf("Enter the value of K: ");
    scanf("%d", &k);
    printf("Enter Image File name with format: ");
    scanf("%s", ifname);
    printf("Enter where to output the image: ");
    scanf("%s", ofname);
    
    printf("Opening Image: %s\n", ifname);
    int rstat = img_read_rgb(ifname, &m, &n, &Rin, &Gin, &Bin);
    if (rstat != 0) {
        printf("Error! Couldn't open image file.\n");
        return 1;
    }

    size_t pix = (size_t)m * n;
    if (gray(Rin, Gin, Bin, pix, m, n)) {
        printf("Image is Grayscale. Compressing once...\n");
        double *AhR = channel(Rin, m, n, k, max);
        double *AhG = calloc(pix, sizeof(double));
        double *AhB = calloc(pix, sizeof(double));
        cpyv(AhR, AhG, pix);
        cpyv(AhR, AhB, pix);
        double eR = err(Rin, AhR, m, n);
        double eG = err(Gin, AhG, m, n);
        double eB = err(Bin, AhB, m, n);
        double e = sqrt(eR * eR + eG * eG + eB * eB);
        printf("\nApproximation error ∥A - A_k∥_F = %.6f\n", e);
        printf("Grayscale Compression Complete. Output file %s produced\n", ofname);

        unsigned char *out = (unsigned char *)malloc(3 * pix);
        for (size_t i = 0; i < pix; i++) {
            unsigned char val = (unsigned char)(AhR[i] > 255 ? 255 : (AhR[i] < 0 ? 0 : AhR[i]));
            out[3 * i] = out[3 * i + 1] = out[3 * i + 2] = val;
        }

        if (strstr(ofname, ".jpg") || strstr(ofname, ".jpeg")) {
            jpg_write_rgb(ofname, m, n, q, AhR, AhG, AhB);
        } else {
            png_write_rgb(ofname, m, n, AhR, AhG, AhB);
        }

        
        free(out); free(AhR); free(AhG); free(AhB);
    } else {
        printf("Image is RGB coloured. Compressing channel-wise...\n");
        double *AhR = channel(Rin, m, n, k, max);
        double *AhG = channel(Gin, m, n, k, max);
        double *AhB = channel(Bin, m, n, k, max);
        double eR = err(Rin, AhR, m, n);
        double eG = err(Gin, AhG, m, n);
        double eB = err(Bin, AhB, m, n);
        double e = sqrt(eR * eR + eG * eG + eB * eB);
        printf("\nApproximation error ∥A - A_k∥_F = %.6f\n", e);
        printf("RGB Compression Complete. Output file %s produced\n", ofname);

        unsigned char *out = (unsigned char *)malloc(3 * pix);
        for (size_t i = 0; i < pix; i++) {
            out[3 * i]     = (unsigned char)(AhR[i] > 255 ? 255 : (AhR[i] < 0 ? 0 : AhR[i]));
            out[3 * i + 1] = (unsigned char)(AhG[i] > 255 ? 255 : (AhG[i] < 0 ? 0 : AhG[i]));
            out[3 * i + 2] = (unsigned char)(AhB[i] > 255 ? 255 : (AhB[i] < 0 ? 0 : AhB[i]));
        }

        if (strstr(ofname, ".jpg") || strstr(ofname, ".jpeg")) {
            jpg_write_rgb(ofname, m, n, q, AhR, AhG, AhB);
        } else {
            png_write_rgb(ofname, m, n, AhR, AhG, AhB);
        }
        
        free(out); free(AhR); free(AhG); free(AhB);
    }

    free(Rin); free(Gin); free(Bin);
    return 0;
}
