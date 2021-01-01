#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>

int allocate_matrix(matrix **mat, int rows, int cols) {

    if (rows < 1 || cols < 1) {
        // PyErr_SetString(PyExc_ValueError, "Invalid rows or cols");
        return -1;
    }
    matrix *new_mat = (matrix *) malloc(sizeof(matrix));
    if (new_mat == NULL) {
        // PyErr_SetString(PyExc_RuntimeError, "Out of memory");
        return -1;
    }

    new_mat->data_1d = (double *) calloc(rows * cols, sizeof(double));
    if (new_mat->data_1d == NULL) {
        free(new_mat);
        // PyErr_SetString(PyExc_RuntimeError, "Out of memory");
        return -1;
    }
    
    new_mat->data = (double **) malloc(sizeof(double *) * (rows));
    if (new_mat->data == NULL) {
        free(new_mat->data_1d);
        free(new_mat);
        // PyErr_SetString(PyExc_RuntimeError, "Out of memory");
        return -1;
    }

    for (int i = 0; i < rows; i += 4) {
        new_mat->data[i] = &(new_mat->data_1d[i * cols]);
    }

    new_mat->rows = rows;
    new_mat->cols = cols;

    if (rows == 1 || cols == 1) {
        new_mat->is_1d = 1;
    } else {
        new_mat->is_1d = 0;
    }

    new_mat->parent = NULL;
    new_mat->ref_cnt = 1;

    *mat = new_mat;
    return 0;
}

int allocate_matrix_ref(matrix **mat, matrix *from, int row_offset, int col_offset,
                        int rows, int cols) {

    matrix *new_mat = (matrix *) malloc(sizeof(matrix));
    if (new_mat == NULL) {
        return -1;
    }

    new_mat->data = (double **) malloc(sizeof(double *) * (rows));
    if (new_mat->data == NULL) {
        free(new_mat);
        return -1;
    }

    for (int i = 0; i < rows; i++) {
        new_mat->data[i] = &from->data[row_offset + i][col_offset];
    }

    new_mat->rows = rows;
    new_mat->cols = cols;

    if (rows == 1 || cols == 1) {
        new_mat->is_1d = 1;
    } else {
        new_mat->is_1d = 0;
    }



    while (from->parent != NULL) { // pointing directly to its parent
        from = from->parent;
    }
    new_mat->parent = from; 

    from->ref_cnt += 1;
    new_mat->ref_cnt = 1;

    // Don't set slice's parent to its direct parent,
    // end up with long chain, which is hard when deallocating it

    *mat = new_mat;
    return 0;
}

void deallocate_matrix(matrix *mat) {
    if (mat == NULL) {
        return;
    } else if (mat->parent == NULL) { // I am the parent/source matrix
        if (mat->ref_cnt == 1) { // I have no child
            free(mat->data);
            free(mat->data_1d);
            free(mat);
        } else { // I have child
            mat->ref_cnt--;
            free(mat->data);
            free(mat);
        }
        // else I have children, so free nothing
    } else if (mat->parent != NULL) { // I am a child/sliced matrix
        if (mat->ref_cnt == 1) { // I am one of my parent's child
            if (mat->parent->ref_cnt <= 1) { // I am my parent's only child and we already deallocated my parent
                free(mat->parent->data_1d);
            }
            mat->parent->ref_cnt--;
            free(mat->data);
            free(mat);
        }
    }
    // automatically called by python
    // ex. if c = a + b
    // and you set c = None
    // python will call deallocate on c for us
}

double get(matrix *mat, int row, int col) {
    return mat->data[row][col];
}

/*
 * Set the value at the given row and column to val. You may assume `row` and
 * `col` are valid
 */
void set(matrix *mat, int row, int col, double val) {
    mat->data[row][col] = val;
}

int add_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    int rows, cols;
    rows = mat1->rows;
    cols = mat1->cols;

    if ((mat1->cols != mat2->cols) || (mat1->rows != mat2->rows)) {
        // PyErr_SetString(PyExc_ValueError, "Cannot add matrix with invalid dimensions");
        return -1;
    }

    // #pragma omp parallel for
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result->data[i][j] = mat1->data[i][j] + mat2->data[i][j];
        }
    }
    return 0;
}

void add_test(void) {
    matrix *result = NULL;
    matrix *mat1 = NULL;
    matrix *mat2 = NULL;
    int rows = 8;
    int cols = 8;
    printf("Allocating result matrix should be 0; status is: %d\n", allocate_matrix(&result, rows,cols));
    printf("Allocating matrix 1 should be 0; status is: %d\n", allocate_matrix(&mat1, rows, cols));
    printf("Allocating matrix 2 should be 0; status is: %d\n", allocate_matrix(&mat2, rows, cols));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set(mat1, i, j, i * cols + j);
            set(mat2, i, j, i * cols + j);
        }
    }
    add_matrix(result, mat1, mat2);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("Value at i: %d, j: %d should be: %d, your value is: %f\n", i, j, 2 * (i * cols + j), get(result, i, j));
            // CU_ASSERT_EQUAL(get(result, i, j), 2 * (i * 2 + j));
        }
    }

    int rows1 = 4;
    int cols1 = 4;
    int start1 = 0;

    matrix *ref = NULL;
    printf("Allocating ref matrix from result [0:4, 0:4] should be 0; status is: %d\n", allocate_matrix_ref(&ref, result, start1, start1, rows1, cols1));

    for (int i = start1; i < start1 + rows1; i++) {
        for (int j = start1; j < start1 + cols1; j++) {
            printf("Value of ref at i: %d, j: %d should be: %d, your value is: %f\n", i, j, 2 * (i * cols + j), get(ref, i, j));
            // CU_ASSERT_EQUAL(get(result, i, j), 2 * (i * 2 + j));
        }
    }

    int rows2 = 4;
    int cols2 = 4;
    int start2 = 4;

    matrix *ref1 = NULL;
    printf("Allocating ref2 matrix from result [4:8, 4:8] should be 0; status is: %d\n", allocate_matrix_ref(&ref1, result, start2, start2, rows2, cols2));

    for (int i = start2; i < start2 + rows2; i++) {
        for (int j = start2; j < start2 + cols2; j++) {
            printf("Value of ref at i: %d, j: %d should be: %d, your value is: %f\n", i, j, 2 * (i * cols + j), get(ref1, i - start2, j - start2));
        }
    }

    int rows3 = 2;
    int cols3 = 2;
    int start3 = 0;

    matrix *nested_ref = NULL;
    printf("Allocating nested_ref matrix from ref [0:2, 0:2] should be 0; status is: %d\n", allocate_matrix_ref(&nested_ref, ref, start3, start3, rows3, cols3));

    for (int i = start3; i < start3 + rows3; i++) {
        for (int j = start3; j < start3 + cols3; j++) {
            printf("Value of ref at i: %d, j: %d should be: %d, your value is: %f\n", i, j, 2 * (i * cols1 + j), get(nested_ref, i, j));
            // CU_ASSERT_EQUAL(get(result, i, j), 2 * (i * 2 + j));
        }
    }
    int rows4 = 1;
    int cols4 = 2;
    int start4 = 1;

    matrix *nested_ref1 = NULL;
    printf("Allocating nested_ref2 matrix from nested_ref [0, 0:2] should be 0; status is: %d\n", allocate_matrix_ref(&nested_ref1, nested_ref, start4, start4, rows4, cols4));

    for (int i = start4; i < rows4 + start4; i++) {
        for (int j = start4; j < cols4 + start4; j++) {
            printf("Value of ref at i: %d, j: %d should be: %d, your value is: %f\n", i, j, 2 * (i * 4 + j), get(nested_ref1, i - start4, j - start4));
            // CU_ASSERT_EQUAL(get(result, i, j), 2 * (i * 2 + j));
        }
    }

    printf("Address of result: %p\n", result);
    printf("Address of ref->parent: %p\n", ref->parent);
    printf("Address of nested_ref->parent: %p\n", nested_ref->parent);


    deallocate_matrix(ref);
    deallocate_matrix(result);
    deallocate_matrix(nested_ref);
    deallocate_matrix(mat1);
    deallocate_matrix(ref1);
    deallocate_matrix(mat2);
    deallocate_matrix(nested_ref1);

}

int main(void) {
    add_test();
}









