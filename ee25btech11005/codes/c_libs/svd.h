#ifndef SVD_H
#define SVD_H

double** crmat(int m,int n);
double* crvec(int n);
double vnorm(double *v,int n);
double dot(double *a,double *b,int n);
void Axv(double **a,double *x,double *y,int m,int n);
void scl(double *v,double s,int n);
void ypx(double *x,double *y,double a,int n);
void cpyv(double *x,double *y,int n);
double **trp(double **A,int m,int n);
int loop(double **A,double *Q,double *v,int k,int m,int n,double *a,double *b);
int GRot(double *a,double *b,int k,double *ev,int mit,double lim);
double* channel(double *A1,int m,int n,int k_in,int mxg);

#endif
