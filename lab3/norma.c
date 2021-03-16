#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h> 

#define STEPS 1000

const double eps = 1e-15;

void print_matrix(double **a, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%'.4f\t", a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_vector(double *v, int n) {
    for (int i = 0; i < n; i++) {
        printf("%'.4f\t", v[i]);
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

void copy_vector(double *v, double *cpy, int n) {
    for (int i = 0; i < n; i++) {
        cpy[i] = v[i];
    }
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

double norm_vector(double *v, int n, double p) {
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
        free(tmp);
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
        free(tmp);
    } else {
        printf("Najwieksza obslugiwana macierz to 3x3\n");
    }
}

double norm_matrix(double **m, int n, double p) {
    double *x = vector(n);
    double *mx = vector(n);
    double *best = vector(n);
    double max = 0.0;
    double norm = 0.0;
    if (n == 2) {
        // (a, b, c) będzie w kwadracie o wierzchołkach (-1/1, -1/1), poniewaz jezeli jakakolwiek wspolzedna bedzie wieksza niz 1 to norma bedzie wieksza od 1
        for (double a = -1.0; a <= 1.0; a += 2.0/STEPS/STEPS) {
            // (1 - a^p)^1/p = |b|
            x[0] = a;
            double b = pow(1.0-pow(fabs(a), p), 1.0/p);

            x[1] = b;
            matrix_mul_vector(m, x, mx, n);
            norm = norm_vector(mx, n, p);
            if (norm > max) {
                max = norm;
                copy_vector(x, best, n);
            }

            x[1] = -b;
            matrix_mul_vector(m, x, mx, n);
            norm = norm_vector(mx, n, p);
            if (norm > max) {
                max = norm;
                copy_vector(x, best, n);
            }
        }
    } else if (n == 3) {
        // (a, b, c) będzie w sześcianie o wierzchołkach (-1/1, -1/1, -1/1), poniewaz jezeli jakakolwiek wspolzedna bedzie wieksza niz 1 to norma bedzie wieksza od 1
        for (double a = -1.0; a <= 1.0; a += 2.0/STEPS) {
            for (double b = -1.0; b <= 1.0; b += 2.0/STEPS) {
                // (1 - a^p - b^p)^1/p = |c|
                x[1] = a;
                x[2] = b;
                double c = pow((1.0 - pow(fabs(a), p) - pow(fabs(b), p)), 1.0/p);

                x[0] = c;
                matrix_mul_vector(m, x, mx, n);
                norm = norm_vector(mx, n, p);
                if (norm > max) {
                    max = norm;
                    copy_vector(x, best, n);
                }

                x[0] = -c;
                matrix_mul_vector(m, x, mx, n);
                norm = norm_vector(mx, n, p);
                if (norm > max) {
                    max = norm;
                    copy_vector(x, best, n);
                }
            }
        }
    } else {
        // Nalezałoby stworzyc funkcje rekurencyjna dzialajaca analogicznie do tych wyzej tzn.
        // (1 - x1^p - x2^p - ... - xn^p)^(1-p) = |x|
        // W ten sposob tracilibysmy jeden stopien swobody, co dla n = 2 i n = 3 daje fajne efekty, ale wraz ze wzrostem rozmiaru macierzy zlozonosc rośnie wykładniczo.
        // Złonosc obliczeniowa bedzie STEPS^(n-1), więc wykładniczy wzrost względem rozmiaru macierzy.
        // Aby zwiekszyć szybkość i dokładność rozwiązania nalezaloby najpierw przemnozyc macierz przez wektor z n-1 niewiadomymi, a nastepnie znaleźć maksimum z normy z otrzymanego wektora.
        printf("Najwieksza obslugiwana macierz to 3x3\n");
    }
    printf("Najlepszy wektor:\n");
    print_vector(best, n);
    free(x);
    free(mx);
    free(best);
    return max;
}

int main(int argc, char *argv[]) {
    if (argc != 3 && argc != 4) {
        printf("Przyjmowane argumenty n, p oraz plik\n");
        return 1;
    }
    int n = atoi(argv[1]);
    double p = atof(argv[2]);
    if (n > 3) {
        printf("Najwieksza obslugiwana macierz to 3x3\n");
        return 1;
    }
    double **m;
    if (argc == 3) {
        srand(time(NULL));
        m = matrix_rand(n);
    } else if (argc == 4) {
        printf("\n!!! Macierz odczytana z pliku %s !!!\n\n", argv[3]);
        m = read_matrix(argv[3], n);
    }

    double **inv_m = copy_matrix(m, n);

    printf("Macierz M\n");
    print_matrix(m, n);

    printf("Norma: %f\n\n", norm_matrix(m, n, p));

    invert_matrix(inv_m, n);
    printf("Macierz odwrotna M\n");
    print_matrix(inv_m, n);

    printf("Norma: %f\n\n", norm_matrix(inv_m, n, p));

    printf("Sprawdzenie M*M^(-1):\n");
    double **check = matrix(n);
    matrix_mul_matrix(m, inv_m, check, n);
    print_matrix(check, n);

    free(m);
    free(inv_m);
    free(check);
    return 0;
}