#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define CLOCKS_PER_MILISEC CLOCKS_PER_SEC / 1000

double **reference;



// zwraca roznice pomiedzy dwiema macierzami
double check(double **a, double **b, int n) {
    double diff = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            diff += sqrt(a[i][j]*a[i][j] - b[i][j]*b[i][j]);
        }
    }
    return diff;
}

void print(double **a, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%'.2f\t", a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// zwraca tablice NxM
double **matrix(int n) {
    double **m = calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++) {
        m[i] = calloc(n, sizeof(double));
    }
    return m;
}

double drand(double min, double max) {
    return ((double) rand() * (max - min) ) / (double) RAND_MAX + min;
}

double **matrix_rand(int n) {
    double **m = matrix(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            m[i][j] = drand(-10, 10);
        }
    }
    return m;
}

void mul_0_ijk(double **a, double **b, double **c, int n) {
    register int i, j, k;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                c[i][j] += a[i][k]*b[k][j];
            }
        }
    }
}

void ijk(double **a, double **b, int n) {
    clock_t start, stop;
    double **ijk = matrix(n);

    start = clock();
    mul_0_ijk(a, b, ijk, n);
    stop = clock();

    printf("-------IJK-------\n");
    printf("Time: %'.2f ms\n", (double) (stop - start) / (CLOCKS_PER_SEC/1000));
    printf("Difference: %'.2f\n\n", check(reference, ijk, n));
}

void mul_1_ikj(double **a, double **b, double **c, int n) {
    register int i, j, k;
    for (i = 0; i < n; i++) {
        for (k = 0; k < n; k++) {
            for (j = 0; j < n; j++) {
                c[i][j] += a[i][k]*b[k][j];
            }
        }
    }
}

void ikj(double **a, double **b, int n) {
    clock_t start, stop;
    double **ikj = matrix(n);

    start = clock();
    mul_1_ikj(a, b, ikj, n);
    stop = clock();

    printf("-------IKJ-------\n");
    printf("Time: %'.2f ms\n", (double) (stop - start) / (CLOCKS_PER_SEC/1000));
    printf("Difference: %'.2f\n\n", check(reference, ikj, n));
}

void mul_2_jik(double **a, double **b, double **c, int n) {
    register int i, j, k;
    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) {
            for (k = 0; k < n; k++) {
                c[i][j] += a[i][k]*b[k][j];
            }
        }
    }
}

void jik(double **a, double **b, int n) {
    clock_t start, stop;
    double **jik = matrix(n);

    start = clock();
    mul_2_jik(a, b, jik, n);
    stop = clock();

    printf("-------JIK-------\n");
    printf("Time: %'.2f ms\n", (double) (stop - start) / (CLOCKS_PER_SEC/1000));
    printf("Difference: %'.2f\n\n", check(reference, jik, n));
}

void mul_3_jki(double **a, double **b, double **c, int n) {
    register int i, j, k;
    for (j = 0; j < n; j++) {
        for (k = 0; k < n; k++) {
            for (i = 0; i < n; i++) {
                c[i][j] += a[i][k]*b[k][j];
            }
        }
    }
}

void jki(double **a, double **b, int n) {
    clock_t start, stop;
    double **jki = matrix(n);

    start = clock();
    mul_3_jki(a, b, jki, n);
    stop = clock();

    printf("-------JKI-------\n");
    printf("Time: %'.2f ms\n", (double) (stop - start) / (CLOCKS_PER_SEC/1000));
    printf("Difference: %'.2f\n\n", check(reference, jki, n));
}

void mul_4_kij(double **a, double **b, double **c, int n) {
    register int i, j, k;
    for (k = 0; k < n; k++) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                c[i][j] += a[i][k]*b[k][j];
            }
        }
    }
}

void kij(double **a, double **b, int n) {
    clock_t start, stop;
    double **kij = matrix(n);

    start = clock();
    mul_4_kij(a, b, kij, n);
    stop = clock();

    printf("-------KIJ-------\n");
    printf("Time: %'.2f ms\n", (double) (stop - start) / (CLOCKS_PER_SEC/1000));
    printf("Difference: %'.2f\n\n", check(reference, kij, n));
}

void mul_5_kji(double **a, double **b, double **c, int n) {
    register int i, j, k;
    for (k = 0; k < n; k++) {
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {
                c[i][j] += a[i][k]*b[k][j];
            }
        }
    }
}

void kji(double **a, double **b, int n) {
    clock_t start, stop;
    double **kji = matrix(n);

    start = clock();
    mul_5_kji(a, b, kji, n);
    stop = clock();

    printf("-------KJI-------\n");
    printf("Time: %'.2f ms\n", (double) (stop - start) / (CLOCKS_PER_SEC/1000));
    printf("Difference: %'.2f\n\n", check(reference, kji, n));
}



int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Require size argument\n");
        return 1;
    }
    int n = atoi(argv[1]);
    srand(time(NULL));

    double **a = matrix_rand(n);
    // print(a, n);
    double **b = matrix_rand(n);
    // print(b, n);

    // Wynik wzorcowy
    reference = matrix(n); mul_0_ijk(a, b, reference, n);

    printf("===========%d===========\n", n);
    ijk(a, b, n);
    ikj(a, b, n);
    jik(a, b, n);
    jki(a, b, n);
    kij(a, b, n);
    kji(a, b, n);
    printf("========================\n\n\n");

    return 0;
}