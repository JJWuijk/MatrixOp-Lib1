#include "matrix_op.h"

static double det2(double a, double b, double c, double d) {
    return a * d - b * c;
}

void mat_add(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double R[SIZE][SIZE]) {
	int r,c;
    for (r = 0; r < SIZE; r++)
        for (c = 0; c < SIZE; c++)
            R[r][c] = A[r][c] + B[r][c];
}

void mat_sub(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double R[SIZE][SIZE]) {
    int r,c;
	for (r = 0; r < SIZE; r++)
        for (c = 0; c < SIZE; c++)
            R[r][c] = A[r][c] - B[r][c];
}

void mat_elem_mul(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double R[SIZE][SIZE]) {
    int r,c;
	for (r = 0; r < SIZE; r++)
        for (c = 0; c < SIZE; c++)
            R[r][c] = A[r][c] * B[r][c];
}

void mat_mul(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double R[SIZE][SIZE]) {
    int r,c,k;
	for (r = 0; r < SIZE; r++) {
        for (c = 0; c < SIZE; c++) {
            double sum = 0.0;
            for (k = 0; k < SIZE; k++) sum += A[r][k] * B[k][c];
            R[r][c] = sum;
        }
    }
}

void mat_transpose(const double A[SIZE][SIZE], double R[SIZE][SIZE]) {
    int r,c;
	for (r = 0; r < SIZE; r++)
        for (c = 0; c < SIZE; c++)
            R[c][r] = A[r][c];
}

double mat_det(const double A[SIZE][SIZE]) {
    const double a = A[0][0], b = A[0][1], c = A[0][2];
    const double d = A[1][0], e = A[1][1], f = A[1][2];
    const double g = A[2][0], h = A[2][1], i = A[2][2];

    return a * (e * i - f * h)
         - b * (d * i - f * g)
         + c * (d * h - e * g);
}

void mat_adjoint(const double A[SIZE][SIZE], double adj[SIZE][SIZE]) {
    int r,c;
	double C[SIZE][SIZE];

    C[0][0] =  det2(A[1][1], A[1][2], A[2][1], A[2][2]);
    C[0][1] = -det2(A[1][0], A[1][2], A[2][0], A[2][2]);
    C[0][2] =  det2(A[1][0], A[1][1], A[2][0], A[2][1]);

    C[1][0] = -det2(A[0][1], A[0][2], A[2][1], A[2][2]);
    C[1][1] =  det2(A[0][0], A[0][2], A[2][0], A[2][2]);
    C[1][2] = -det2(A[0][0], A[0][1], A[2][0], A[2][1]);

    C[2][0] =  det2(A[0][1], A[0][2], A[1][1], A[1][2]);
    C[2][1] = -det2(A[0][0], A[0][2], A[1][0], A[1][2]);
    C[2][2] =  det2(A[0][0], A[0][1], A[1][0], A[1][1]);

    for (r = 0; r < SIZE; r++)
        for (c = 0; c < SIZE; c++)
            adj[r][c] = C[c][r];
}

int mat_inverse(const double A[SIZE][SIZE], double inv[SIZE][SIZE]) {
    int r,c;
	const double detA = mat_det(A);
    if (detA == 0.0) return 0;

    double adj[SIZE][SIZE];
    mat_adjoint(A, adj);

    for (r = 0; r < SIZE; r++)
        for (c = 0; c < SIZE; c++)
            inv[r][c] = adj[r][c] / detA;

    return 1;
}

