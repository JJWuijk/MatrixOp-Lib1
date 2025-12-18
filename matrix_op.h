#ifndef MATRIX_OP_H
#define MATRIX_OP_H

#define SIZE 3

void mat_add(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double R[SIZE][SIZE]);
void mat_sub(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double R[SIZE][SIZE]);
void mat_elem_mul(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double R[SIZE][SIZE]);

void mat_mul(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double R[SIZE][SIZE]);
void mat_transpose(const double A[SIZE][SIZE], double R[SIZE][SIZE]);

double mat_det(const double A[SIZE][SIZE]);
void mat_adjoint(const double A[SIZE][SIZE], double adj[SIZE][SIZE]);

int mat_inverse(const double A[SIZE][SIZE], double inv[SIZE][SIZE]);

#endif

