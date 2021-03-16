#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h> 

const double eps = 1e-15;

void print(double **a, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%'.1f\t", a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
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

void gauss_ones(double **u, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = n-1; j >= i; j--) {
            u[i][j] /= u[i][i];
        }
        for (int j = i + 1; j < n && i != n - 1; j++) {
            if (fabs(u[i][i]) < eps) {
                printf("Wiersz %d wymaga pivotingu", i);
                return;
            }
            double m = -u[j][i];
            for (int k = i; k < n; k++) {
                u[j][k] += m * u[i][k];
            }
        }
    }
}

void gauss(double **u, int n) {
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (fabs(u[i][i]) < eps) {
                printf("Wiersz %d wymaga pivotingu", i);
                return;
            }
            double m = -u[j][i] / u[i][i];
            for (int k = i; k < n; k++) {
                u[j][k] += m * u[i][k];
            }
        }
    }
}

void gauss_pivot(double **u, int n) {
    for (int i = 0; i < n; i++) {
        int p = i;
        double max_p = 0.0;
        for (int j = i; j < n; j++) {
            if (max_p < fabs(u[j][i])) {
                max_p = fabs(u[j][i]);
                p = j;
            }
        }

        if (p != i) {
            for (int j = i; j < n; j++) {
                double tmp = u[p][j];
                u[p][j] = u[i][j];
                u[i][j] = tmp;
            }
        }

        for (int j = i +1; j < n; j++) {
            double m = -u[j][i] / u[i][i];
            for (int k = i; k < n; k++) {
                u[j][k] += m * u[i][k];
            }
        }
    }
}

void lu(double **a, double **l, double **u, int n) {
    for(int j = 0; j < n; j++) {
        for(int i = 0; i < n; i++) {
            if(i <= j) {
                u[i][j] = a[i][j];
                for (int k = 0; k < i-1; k++) {
                    u[i][j] -= l[i][k]*u[k][j];
                }
                if (i == j) l[i][j] = 1;
                else l[i][j] = 0;
            } else {
                l[i][j] = a[i][j];
                for(int k = 0; k <= j - 1; k++) {
                    l[i][j] -= l[i][k] * u[k][j];
                }
                l[i][j] /= u[j][j];
                u[i][j] = 0;
            }
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc == 2) {
        int n = atoi(argv[1]);

        srand(time(NULL));
        double **u = matrix_rand(n);

        print(u, n);
        printf("Algorytm eliminacji Gaussa bez pivotingu - jedynki na przekątnej\n");
        double **a = copy_matrix(u, n);
        gauss_ones(a, n);
        print(a, n);
        printf("Algorytm eliminacji Gaussa bez pivotingu\n");
        double **b = copy_matrix(u, n);
        gauss(b, n);
        print(b, n);
        printf("Algorytm eliminacji Gaussa z pivotigniem\n");
        double **c = copy_matrix(u, n);
        gauss_pivot(c, n);
        print(c, n);
        printf("Algorytm LU faktoryzacji bez pivotingu\n");
        double **d = copy_matrix(u, n);
        double **l = matrix(n);
        double **u1 = matrix(n);
        lu(d, l, u1, n);
        printf("L:\n");
        print(l, n);
        printf("U:\n");
        print(u1, n);
    } else if (argc == 3) {
        int n = atoi(argv[1]);

        printf("!!! Macierz odczytana z pliku %s !!!\n\n", argv[2]);
        double **u = read_matrix(argv[2], n);

        print(u, n);
        printf("Algorytm eliminacji Gaussa bez pivotingu - jedynki na przekątnej\n");
        double **a = copy_matrix(u, n);
        gauss_ones(a, n);
        print(a, n);
        printf("Algorytm eliminacji Gaussa bez pivotingu\n");
        double **b = copy_matrix(u, n);
        gauss(b, n);
        print(b, n);
        printf("Algorytm eliminacji Gaussa z pivotigniem\n");
        double **c = copy_matrix(u, n);
        gauss_pivot(c, n);
        print(c, n);
        printf("Algorytm LU faktoryzacji bez pivotingu\n");
        double **d = copy_matrix(u, n);
        double **l = matrix(n);
        double **u1 = matrix(n);
        lu(d, l, u1, n);
        printf("L:\n");
        print(l, n);
        printf("U:\n");
        print(u1, n);
    } else {
        printf("Require size argument\n");
        return 1;
    }


    return 0;
}