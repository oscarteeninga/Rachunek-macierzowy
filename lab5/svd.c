#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double eps = 0.01;

void print_vector(double *v, int n) {
    for (int i = 0; i < n; i++) {
        printf("\t%'.4f", v[i]);
    }
    printf("\n");
}

void print_matrix(double **a, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("\t%'.4f", a[i][j]);
        }
        printf("\n");
    }
}

double *vector(int n) {
    double *v = calloc(n, sizeof(double));
    for (int i = 0; i < n; i++) {
        v[i] = 100.0;
    }
    return v;
}

double **matrix(int n) {
    double **m = calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++) {
        m[i] = calloc(n, sizeof(double));
    }
    return m;
}

double **read_matrix(char *name, int n) {
    double **m = matrix(n);
    FILE *file = fopen(name, "r");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fscanf(file, "%lf", &m[i][j]);
        }
    }
    return m;
}

void vector_div_numeric(double *v, double a, double *r, int n) {
    for (int i = 0; i < n; i++) {
        r[i] = v[i] / a;
    }
}

void vector_mul_vector(double *v1, double *v2, double **m, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            m[i][j] = v1[i]*v2[j];
        }
    }
}

void matrix_mul_vector(double **m, double *v, double *r, int n) {
    register int i, j;
    for (i = 0; i < n; i++) {
        r[i] = 0;
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            r[j] += m[j][i]*v[i];
        }
    }
}

void matrix_min_matrix(double **m1, double **m2, double **r, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            r[i][j] = m1[i][j] - m2[i][j];
        }
    }
}

double norm_vector(double *v, int n) {
    double norm = fabs(v[0]);
    for (int i = 1; i < n; i++) {
        if (fabs(v[i]) > norm) {
            norm = fabs(v[i]);
        }
    }
    return norm;
}

double error_l(double *w, double *z, double l, int n) {
    double *r = vector(n);
    for (int i = 0; i < n; i++) {
        r[i] = w[i] - l * z[i];
    }
    return norm_vector(r, n);
}

double z_vector(double **a, double *z, int n) {
    double *w = vector(n);
    double l;
    double err = 2*eps;
    for (int i = 1; ;i++) {
        matrix_mul_vector(a, z, w, n);
        l = norm_vector(w, n);
        err = error_l(w, z, l, n);

        vector_div_numeric(w, l, z, n);
        if (err < eps || i > 10) {
            printf("\n-------Iteracja %d-------\n", i);
            printf("Wartosc wlasna: %f\n", l);
            printf("Wektor z:");
            print_vector(z, n);
            printf("Błąd:\t%f\n", err);
            break;
        }
    }
    return l;
}

void e_matrix(double **a, double *z, double **e, int n) {
    double *az = vector(n);
    matrix_mul_vector(a, z, az, n);
    double **azazt = matrix(n);
    vector_mul_vector(az, az, azazt, n);
    matrix_min_matrix(a, azazt, e, n);
}


void svd(double **a, int n) {
    double *z_prev = vector(n);
    double lambda = z_vector(a, z_prev, n);
    for (int i = 0; i < n; i++) {
        printf("\n\n=======SIGMA %d=======\n", i+1);
        printf("Vector z_prev:");
        print_vector(z_prev, n);
        double **e = matrix(n);
        double *z = vector(n);
        e_matrix(a, z_prev, e, n);
        printf("Macierz E:\n");
        print_matrix(e, n);
        lambda = z_vector(e, z, n);
        for (int i = 0; i < n; i++) {
            z_prev[i] = z[i];
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Require size of matrix!\n");
    }
    int n = atoi(argv[2]);
    double **a = read_matrix(argv[1], n);

    printf("Macierz odczytana z pliku %s:\n", argv[1]);
    print_matrix(a, n);
    printf("\n");
    svd(a, n);
    return 0;
}