#include <stdio.h>
#include <math.h>
#include "matrix_op.h"



//dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd

static void print_mat(const char *name, const double M[SIZE][SIZE]) {
    printf("%s =\n", name);
    int r,c;
    for (r = 0; r < SIZE; r++) {
        for (c = 0; c < SIZE; c++) {
            printf("%10.6f ", M[r][c]);
        }
        printf("\n");
    }
    printf("\n");
}

static int approx_identity(const double M[SIZE][SIZE], double eps) {
	int r,c;
    for (r = 0; r < SIZE; r++) {
        for (c = 0; c < SIZE; c++) {
            double target = (r == c) ? 1.0 : 0.0;
            if (fabs(M[r][c] - target) > eps) return 0;
        }
    }
    return 1;
}

int main(void) {
    double A[SIZE][SIZE] = {
        { 4,  7,  2},
        { 3,  6,  1},
        { 2,  5,  3}
    };

    double B[SIZE][SIZE] = {
        { 1,  0,  2},
        {-1,  3,  1},
        { 4,  2,  0}
    };

    double R[SIZE][SIZE];
    double T[SIZE][SIZE];
    double invA[SIZE][SIZE];
    double check[SIZE][SIZE];

    print_mat("A", A);
    print_mat("B", B);

    mat_add(A, B, R);
    print_mat("A + B", R);

    mat_sub(A, B, R);
    print_mat("A - B", R);

    mat_elem_mul(A, B, R);
    print_mat("A o B (element-wise)", R);

    mat_mul(A, B, R);
    print_mat("A * B", R);

    mat_transpose(A, T);
    print_mat("A^T", T);

    double detA = mat_det(A);
    printf("det(A) = %.6f\n\n", detA);

    if (!mat_inverse(A, invA)) {
        printf("A is not invertible (det(A) == 0)\n");
        return 0;
    }

    print_mat("inv(A)", invA);

    mat_mul(A, invA, check);
    print_mat("A * inv(A)", check);

    if (approx_identity(check, 1e-9)) {
        printf("Check PASSED: A * inv(A) is approximately Identity.\n");
    } else {
        printf("Check FAILED: A * inv(A) is NOT close enough to Identity.\n");
    }

    return 0;
}

