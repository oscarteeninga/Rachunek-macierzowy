#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h> 

const double eps = 1e-15;

void print_matrix(double **a, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%'.2f\t", a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_vector(double *v, int n) {
    for (int i = 0; i < n; i++) {
        printf("%'.2f\t", v[i]);
    }
    printf("\n\n");
}

double *vector(int n) {
    double *v = calloc(n, sizeof(double));
    for (int i = 0; i < n; i++) {
        v[i] = 1.0;
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

double **copy_matrix(double **m, int n) {
    double **new = matrix(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            new[i][j] = m[i][j];
        }
    }
    return new;
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

void vector_mul_numeric(double *v, double a, int n) {
    for (int i = 0; i < n; i++) {
        v[i] *= a;
    }
}

 void matrix_mul_numeric(double **m, double a, int n) {
     for (int i = 0; i < n; i++) {
         vector_mul_numeric(m[i], a, n);
     }
 }

 void matrix_mul_matrix(double **a, double **b, double **c, int n) {
    register int i, j, k;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                c[i][j] += a[i][k]*b[k][j];
            }
        }
    }
}

double *matrix_mul_vector(double **m, double *v, int n) {
    double *r = calloc(n, sizeof(double));
    register int i, j, k;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            r[i] += m[i][j]*v[j];
        }
    }
    return r;
}

double norma_vector(double *v, int n, double p) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += pow(fabs(v[i]), p);
    }
    return pow(sum, 1.0/p);
}

double max_vector(double *v, int n) {
    double max = v[0];
    for (int i = 1; i < n; i++) {
        if (v[i] > max) {
            max = v[i];
        }
    }
    return max;
}

void invert_matrix(double **m, int n) {
    if (n == 2) {
        double det = m[0][0]*m[1][1]-m[0][1]*m[1][0];
        if (det == 0) {
            printf("Matrix det is 0\n");
            return;
        }
        double **tmp = matrix(2);
        tmp[0][0] = m[1][1];
        tmp[0][1] = -m[0][1];
        tmp[1][0] = -m[1][0];
        tmp[1][1] = m[0][0];
        matrix_mul_numeric(tmp, 1.0/det, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                m[i][j] = tmp[i][j];
            }
        }
    } else if (n == 3) {
        double det = m[0][0]*m[1][1]*m[2][2] + m[0][1]*m[1][2]*m[2][0] + m[0][2]*m[1][0]*m[2][1] - m[2][0]*m[1][1]*m[0][2] - m[2][1]*m[1][2]*m[0][0] - m[2][2]*m[1][0]*m[0][1];
        if (det == 0) {
            printf("Matrix det is 0\n");
            return;
        }
        double **tmp = matrix(3);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
			    tmp[i][j] = m[(j+1)%3][(i+1)%3] * m[(j+2)%3][(i+2)%3] - m[(j+1)%3][(i+2)%3] * m[(j+2)%3][(i+1)%3];
            }
        }
        matrix_mul_numeric(tmp, 1.0/det, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                m[i][j] = tmp[i][j];
            }
        }
    } else {
        printf("Implementation cannot handle matrix bigger than 3x3!\n");
        return;
    }
}

double norma_matrix(double **m, int n, double p) {
    double *x = vector(n);
    printf("Vector x:\n");
    print_vector(x, n);
    vector_mul_numeric(x, 1.0/norma_vector(x, n, p), n);
    printf("Vector x normalized:\n");
    print_vector(x, n);
    double *mx = matrix_mul_vector(m, x, n);
    printf("Vector Ax:\n");
    print_vector(mx, n);
    return max_vector(mx, n);
}

int main(int argc, char *argv[]) {
    if (argc == 3) {
        int n = atoi(argv[1]);
        double p = atof(argv[2]);

        srand(time(NULL));

        double **m = matrix_rand(n);
        double **inv_m = copy_matrix(m, n);
        printf("===Matrix M===\n");
        print_matrix(m, n);

        printf("Norm: %f\n\n", norma_matrix(m, n, p));

        if (n == 2 || n == 3) {
            invert_matrix(inv_m, n);
            printf("===Inverted M===\n");
            print_matrix(inv_m, n);

            printf("Norm: %f\n\n", norma_matrix(inv_m, n, p));

            printf("Check M*M^(-1):\n");
            double **check = matrix(n);
            matrix_mul_matrix(m, inv_m, check, n);
            print_matrix(check, n);
        }
        return 0;
    } else if (argc == 4) {
        int n = atoi(argv[1]);
        double p = atof(argv[2]);

        printf("!!! Macierz odczytana z pliku %s !!!\n\n", argv[3]);

        double **m = read_matrix(argv[3], n);
        double **inv_m = copy_matrix(m, n);
        printf("===Matrix M===\n");
        print_matrix(m, n);

        printf("Norm: %f\n\n", norma_matrix(m, n, p));

        if (n == 2 || n == 3) {
            invert_matrix(inv_m, n);
            printf("===Inverted M===\n");
            print_matrix(inv_m, n);

            printf("Norm: %f\n\n", norma_matrix(inv_m, n, p));

            printf("Check M*M^(-1):\n");
            double **check = matrix(n);
            matrix_mul_matrix(m, inv_m, check, n);
            print_matrix(check, n);
        }
    } else {
        printf("Require size and p arguments\n");
        return 1;
    }


    return 0;
}