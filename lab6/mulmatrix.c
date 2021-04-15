#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <clapacke.h>
#include <Accelerate/Accelerate.h>

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

void block(double **a, double **b, double **c, int n) {
    int s = 2;
    for (int jj = 0; jj < n; jj += s) {
        for (int kk = 0; kk < n; kk += s) {
            for (int i = 0; i < n; i++) {
                for (int j = jj; j < ((jj+s) > n?n: (jj+s)); j++) {
                    int temp = 0;
                    for (int k = kk; k < ((kk+s) > n?n:(kk+s)); k++) {
                        temp += a[i][k] * b[k][j];
                    }
                    c[i][j] += temp;
                }
            }
        }
    }
}


void SingularValueDecomposition(int m,     // number of rows in matrix
                                int n,     // number of columns in matrix
                                int lda,   // leading dimension of matrix
                                double *a) // pointer to top-left corner
{
    // Setup a buffer to hold the singular values:
    int numberOfSingularValues = m < n ? m : n;
    double *s = malloc(numberOfSingularValues * sizeof s[0]);

    // Setup buffers to hold the matrices U and Vt:
    double *u = malloc(m*m * sizeof u[0]);
    double *vt = malloc(n*n * sizeof vt[0]);

    // Workspace and status variables:
    double workSize;
    double *work = &workSize;
    int lwork = -1;
    int *iwork = malloc(8 * numberOfSingularValues * sizeof iwork[0]);
    int info = 0;

    // Call dgesdd_ with lwork = -1 to query optimal workspace size:
    dgesdd_("A", &m, &n, a, &lda, s, u, &m, vt, &n, work, &lwork, iwork, &info);
    if (info) // handle error conditions here

    // Optimal workspace size is returned in work[0].
    lwork = workSize;
    work = malloc(lwork * sizeof work[0]);

    // Call dgesdd_ to do the actual computation:
    dgesdd_("A", &m, &n, a, &lda, s, u, &m, vt, &n, work, &lwork, iwork, &info);
    if (info) // handle error conditions here

    // Cleanup workspace:
    free(work);
    free(iwork);

    // do something useful with U, S, Vt ...

    // and then clean them up too:
    free(s);
    free(u);
    free(vt);
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Require size of matrix!\n");
    }
    int n = atoi(argv[3]);
    double **m1 = read_matrix(argv[1], n);
    double **m2 = read_matrix(argv[2], n);
    double **m3 = matrix(n);

    block(m1, m2, m3, n);
    print_matrix(m3, n);

    SingularValueDecomposition(n, n, n, m3);
    print_matrix(m3, n);

    return 0;
}