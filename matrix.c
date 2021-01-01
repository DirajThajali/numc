#include "matrix.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Include SSE intrinsics
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#include <x86intrin.h>
#endif

#define min(a,b) (((a)<(b))?(a):(b))

/* Below are some intel intrinsics that might be useful
 * void _mm256_storeu_pd (double * mem_addr, __m256d a)
 * __m256d _mm256_set1_pd (double a)
 * __m256d _mm256_set_pd (double e3, double e2, double e1, double e0)
 * __m256d _mm256_loadu_pd (double const * mem_addr)
 * __m256d _mm256_add_pd (__m256d a, __m256d b)
 * __m256d _mm256_sub_pd (__m256d a, __m256d b)
 * __m256d _mm256_fmadd_pd (__m256d a, __m256d b, __m256d c)
 * __m256d _mm256_mul_pd (__m256d a, __m256d b)
 * __m256d _mm256_cmp_pd (__m256d a, __m256d b, const int imm8)
 * __m256d _mm256_and_pd (__m256d a, __m256d b)
 * __m256d _mm256_max_pd (__m256d a, __m256d b)
*/

/*
 * Generates a random double between `low` and `high`.
 */
double rand_double(double low, double high) {
    double range = (high - low);
    double div = RAND_MAX / range;
    return low + (rand() / div);
}

/*
 * Generates a random matrix with `seed`.
 */
void rand_matrix(matrix *result, unsigned int seed, double low, double high) {
    srand(seed);
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            set(result, i, j, rand_double(low, high));
        }
    }
}

/*
 * Allocate space for a matrix struct pointed to by the double pointer mat with
 * `rows` rows and `cols` columns. You should also allocate memory for the data array
 * and initialize all entries to be zeros. Remember to set all fieds of the matrix struct.
 * `parent` should be set to NULL to indicate that this matrix is not a slice.
 * You should return -1 if either `rows` or `cols` or both have invalid values, or if any
 * call to allocate memory in this function fails. If you don't set python error messages here upon
 * failure, then remember to set it in numc.c.
 * Return 0 upon success and non-zero upon failure.
 */
int allocate_matrix(matrix **mat, int rows, int cols) {

    if (rows < 1 || cols < 1) {
        PyErr_SetString(PyExc_ValueError, "Invalid rows or cols");
        return -1;
    }
    matrix *new_mat = (matrix *) malloc(sizeof(matrix));
    if (new_mat == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Out of memory");
        return -1;
    }

    new_mat->data_1d = (double *) calloc(rows * cols, sizeof(double));
    if (new_mat->data_1d == NULL) {
        free(new_mat);
        PyErr_SetString(PyExc_RuntimeError, "Out of memory");
        return -1;
    }

    new_mat->data = (double **) malloc(sizeof(double *) * (rows));
    if (new_mat->data == NULL) {
        free(new_mat->data_1d);
        free(new_mat);
        PyErr_SetString(PyExc_RuntimeError, "Out of memory");
        return -1;
    }

    for (int i = 0; i < rows; i++) {
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

/*
 * Allocate space for a matrix struct pointed to by `mat` with `rows` rows and `cols` columns.
 * This is equivalent to setting the new matrix to be
 * from[row_offset:row_offset + rows, col_offset:col_offset + cols]
 * If you don't set python error messages here upon failure, then remember to set it in numc.c.
 * Return 0 upon success and non-zero upon failure.
 */
int allocate_matrix_ref(matrix **mat, matrix *from, int row_offset, int col_offset,
                        int rows, int cols) {

    matrix *new_mat = (matrix *) malloc(sizeof(matrix));
    if (new_mat == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Out of memory");
        return -1;
    }

    new_mat->data = (double **) malloc(sizeof(double *) * (rows));
    if (new_mat->data == NULL) {
        free(new_mat);
        PyErr_SetString(PyExc_RuntimeError, "Out of memory");
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

    // from->ref_cnt += 1;
    // new_mat->ref_cnt = from->ref_cnt;

    // new_mat->parent = from;

    *mat = new_mat;
    return 0;
}

/*
 * This function will be called automatically by Python when a numc matrix loses all of its
 * reference pointers.
 * You need to make sure that you only free `mat->data` if no other existing matrices are also
 * referring this data array.
 * See the spec for more information.
 */
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

/*
 * Return the double value of the matrix at the given row and column.
 * You may assume `row` and `col` are valid.
 */
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

/*
 * Set all entries in mat to val
 */
void fill_matrix(matrix *mat, double val) {
    #pragma omp parallel for
    for (int i = 0; i < mat->rows * mat->cols; i++) {
        mat->data_1d[i] = val;
    }
}

/*
 * Store the result of adding mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int add_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    int rows, cols, total_elems;
    rows = mat1->rows;
    cols = mat1->cols;
    total_elems = rows * cols;

    if ((mat1->cols != mat2->cols) || (mat1->rows != mat2->rows)) {
        PyErr_SetString(PyExc_ValueError, "Cannot add matrix with invalid dimensions");
        return -1;
    }
    if (total_elems < 200) {
        for (int i = 0; i < total_elems; i++) {
            result->data_1d[i] = mat1->data_1d[i] + mat2->data_1d[i];
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < total_elems / 16 * 16; i += 16) {
        __m256d a = _mm256_loadu_pd(&(mat1->data_1d[i]));
        __m256d b = _mm256_loadu_pd(&(mat2->data_1d[i]));
        __m256d add_v1 = _mm256_add_pd(a, b);
        _mm256_storeu_pd(&(result->data_1d[i]), add_v1);

        __m256d c = _mm256_loadu_pd(&(mat1->data_1d[i + 4]));
        __m256d d = _mm256_loadu_pd(&(mat2->data_1d[i + 4]));
        __m256d add_v2 = _mm256_add_pd(c, d);
        _mm256_storeu_pd(&(result->data_1d[i + 4]), add_v2);

        __m256d e = _mm256_loadu_pd(&(mat1->data_1d[i + 8]));
        __m256d f = _mm256_loadu_pd(&(mat2->data_1d[i + 8]));
        __m256d add_v3 = _mm256_add_pd(e, f);
        _mm256_storeu_pd(&(result->data_1d[i + 8]), add_v3);

        __m256d g = _mm256_loadu_pd(&(mat1->data_1d[i + 12]));
        __m256d h = _mm256_loadu_pd(&(mat2->data_1d[i + 12]));
        __m256d add_v4 = _mm256_add_pd(g, h);
        _mm256_storeu_pd(&(result->data_1d[i + 12]), add_v4);
    }

    #pragma omp parallel for
    for (int i = total_elems / 16 * 16; i < total_elems; i++) { // tail case
        result->data_1d[i] = mat1->data_1d[i] + mat2->data_1d[i];
    }

    return 0;
}

/*
 * Store the result of subtracting mat2 from mat1 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int sub_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    int rows, cols, total_elems;
    rows = mat1->rows;
    cols = mat1->cols;
    total_elems = rows * cols;

    if ((mat1->cols != mat2->cols) || (mat1->rows != mat2->rows)) {
        PyErr_SetString(PyExc_ValueError, "Cannot subtract matrix with invalid dimensions");
        return -1;
    }
    #pragma omp parallel for
    for (int i = 0; i < total_elems / 16 * 16; i += 16) {
        __m256d a = _mm256_loadu_pd(&(mat1->data_1d[i]));
        __m256d b = _mm256_loadu_pd(&(mat2->data_1d[i]));
        __m256d sub_v1 = _mm256_sub_pd(a, b);
        _mm256_storeu_pd(&(result->data_1d[i]), sub_v1);

        __m256d c = _mm256_loadu_pd(&(mat1->data_1d[i + 4]));
        __m256d d = _mm256_loadu_pd(&(mat2->data_1d[i + 4]));
        __m256d sub_v2 = _mm256_sub_pd(c, d);
        _mm256_storeu_pd(&(result->data_1d[i + 4]), sub_v2);

        __m256d e = _mm256_loadu_pd(&(mat1->data_1d[i + 8]));
        __m256d f = _mm256_loadu_pd(&(mat2->data_1d[i + 8]));
        __m256d sub_v3 = _mm256_sub_pd(e, f);
        _mm256_storeu_pd(&(result->data_1d[i + 8]), sub_v3);

        __m256d g = _mm256_loadu_pd(&(mat1->data_1d[i + 12]));
        __m256d h = _mm256_loadu_pd(&(mat2->data_1d[i + 12]));
        __m256d sub_v4 = _mm256_sub_pd(g, h);
        _mm256_storeu_pd(&(result->data_1d[i + 12]), sub_v4);
    }
    #pragma omp parallel for
    for (int i = total_elems / 16 * 16; i < total_elems; i++) { // tail case
        result->data_1d[i] = mat1->data_1d[i] - mat2->data_1d[i];
    }
    return 0;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 */
int mul_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    int I = mat1->rows;
    int J = mat2->cols;
    int K = mat2->rows;
    int L = mat1->cols;

    if (L != K) {
        PyErr_SetString(PyExc_ValueError, "Cannot multiply matrix with invalid dimensions");
        return -1;
    }

    double *mat1_data = mat1->data_1d;
    double *mat2_data = mat2->data_1d;
    double *result_data = result->data_1d;
//    omp_set_num_threads(16);
    
    #pragma omp parallel for
    for(int i = 0; i < I; i++){
        for(int k = 0; k < K; k++){
            __m256d vA = _mm256_set1_pd(mat1_data[i * K + k]);
            for(int j = 0; j < J / 16 * 16; j += 16){
                __m256d sum = _mm256_loadu_pd(result_data + i * J + j);
                __m256d vB = _mm256_loadu_pd(mat2_data + k * J + j);
                sum = _mm256_fmadd_pd(vA, vB, sum);
                _mm256_storeu_pd(result_data + i * J + j, sum);

                __m256d sum1 = _mm256_loadu_pd(result_data + i * J + j + 4);
                __m256d vB1 = _mm256_loadu_pd(mat2_data + k * J + j + 4);
                sum1 = _mm256_fmadd_pd(vA, vB1, sum1);
                _mm256_storeu_pd(result_data + i * J + j + 4, sum1);

                __m256d sum2 = _mm256_loadu_pd(result_data + i * J + j + 8);
                __m256d vB2 = _mm256_loadu_pd(mat2_data + k * J + j + 8);
                sum2 = _mm256_fmadd_pd(vA, vB2, sum2);
                _mm256_storeu_pd(result_data + i * J + j + 8, sum2);

                __m256d sum3 = _mm256_loadu_pd(result_data + i * J + j + 12);
                __m256d vB3 = _mm256_loadu_pd(mat2_data + k * J + j + 12);
                sum3 = _mm256_fmadd_pd(vA, vB3, sum3);
                _mm256_storeu_pd(result_data + i * J + j + 12, sum3);
             }
             double independent = mat1_data[i * K + k];
             for(int j = J / 16 * 16; j < J; j++){ // tail case
                 result_data[i * J + j] += independent * mat2_data[k * J + j];
             }
         }
     }
    return 0;
}

/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that pow is defined with matrix multiplication, not element-wise multiplication.
 */
int pow_matrix(matrix *result, matrix *mat, int pow) {
    // iterative solution is faster because you don't have to do all the setting up register args for function calls
    int rows, cols, total_elems;
    rows = mat->rows;
    cols = mat->cols;
    total_elems = rows * cols;
    if (mat->rows != mat->cols) {
        PyErr_SetString(PyExc_ValueError, "Pow matrix must be a square matrix");
        return -1;
    }
    if (pow < 0) {
        PyErr_SetString(PyExc_ValueError, "Pow must be a non-negative integer");
        return -1;
    }
    fill_matrix(result, 0);
    for (int i = 0; i < rows; i++) {
        set(result, i, i, 1);
    }

    matrix *temp_mat;
    int allocation_failed = allocate_matrix(&temp_mat, rows, cols);
    if (allocation_failed) {
        return allocation_failed;
    }
    for (int i = 0; i < total_elems; i++) {
        temp_mat->data_1d[i] = mat->data_1d[i];
    }

    matrix *temp;
    int allocation_failed1 = allocate_matrix(&temp, rows, cols);
    if (allocation_failed1) {
        free(temp_mat);
        return allocation_failed1;
    }

    while (pow > 0) {
        if (pow & 1) { // pow is odd
            // result = result * temp_mat
            for (int i = 0; i < total_elems; i++) {
                temp->data_1d[i] = result->data_1d[i];
            }
            fill_matrix(result, 0);
            int operation_failed = mul_matrix(result, temp, temp_mat);
            if (operation_failed) {
                free(temp_mat);
                free(temp);
                return operation_failed;
            }
            if (pow == 1) {
                free(temp_mat);
                free(temp);
                return 0;
            }
        }
        pow = pow / 2;
        // temp_mat = temp_mat ^ 2
        fill_matrix(temp, 0); // power is even
        int operation_failed = mul_matrix(temp, temp_mat, temp_mat); // temp = temp_mat * temp_mat
        if (operation_failed) {
            free(temp_mat);
            free(temp);
            return operation_failed;
        }
        for (int i = 0; i < total_elems; i++) {
            temp_mat->data_1d[i] = temp->data_1d[i];
        }
    }
    free(temp);
    free(temp_mat);
    return 0;
}


/*
 * Store the result of element-wise negating mat's entries to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int neg_matrix(matrix *result, matrix *mat) {
    int total_elems = mat->rows * mat->cols;

    double *mat_data = mat->data_1d;
    double *result_data = result->data_1d;

    __m256d zero = _mm256_set1_pd(0);
    #pragma omp parallel for
    for (int i = 0; i < total_elems / 4 * 4; i += 4) {
        __m256d elem_v = _mm256_loadu_pd(&(mat_data[i]));
        __m256d neg_v = _mm256_sub_pd(zero, elem_v);
        _mm256_storeu_pd(&(result_data[i]), neg_v);

        // __m256d elem_v1 = _mm256_loadu_pd(&(mat_data[i + 4]));
        // __m256d neg_v1 = _mm256_sub_pd(zero, elem_v1);
        // _mm256_storeu_pd(&(result_data[i + 4]), neg_v1);

        // __m256d elem_v2 = _mm256_loadu_pd(&(mat_data[i + 8]));
        // __m256d neg_v2 = _mm256_sub_pd(zero, elem_v2);
        // _mm256_storeu_pd(&(result_data[i + 8]), neg_v2);

        // __m256d elem_v3 = _mm256_loadu_pd(&(mat_data[i + 12]));
        // __m256d neg_v3 = _mm256_sub_pd(zero, elem_v3);
        // _mm256_storeu_pd(&(result_data[i + 12]), neg_v3);
    }
    
    for (int i = total_elems / 4 * 4; i < total_elems; i++) { // tail case
        result_data[i] = 0 - mat_data[i];
    }
    return 0;
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int abs_matrix(matrix *result, matrix *mat) {
    
    int total_elems = mat->rows * mat->cols;

    double *mat_data = mat->data_1d;
    double *result_data = result->data_1d;

    __m256d zero = _mm256_set1_pd(0);
    #pragma omp parallel for
    for (int i = 0; i < total_elems / 4 * 4; i += 4) {
        __m256d value = _mm256_loadu_pd(&(mat_data[i]));
        __m256d sub = _mm256_sub_pd(zero, value);
        __m256d res_v = _mm256_and_pd(value, sub);
        _mm256_storeu_pd(&(result_data[i]), res_v);

        // __m256d value1 = _mm256_loadu_pd(&(mat_data[i + 4]));
        // __m256d sub1 = _mm256_sub_pd(zero, value1);
        // __m256d res_v1 = _mm256_and_pd(value1, sub1);
        // _mm256_storeu_pd(&(result_data[i + 4]), res_v1);

        // __m256d value2 = _mm256_loadu_pd(&(mat_data[i + 8]));
        // __m256d sub2 = _mm256_sub_pd(zero, value2);
        // __m256d res_v2 = _mm256_and_pd(value2, sub2);
        // _mm256_storeu_pd(&(result_data[i + 8]), res_v2);

        // __m256d value3 = _mm256_loadu_pd(&(mat_data[i + 12]));
        // __m256d sub3 = _mm256_sub_pd(zero, value3);
        // __m256d res_v3 = _mm256_and_pd(value3, sub3);
        // _mm256_storeu_pd(&(result_data[i + 12]), res_v3);
    }
    
    for (int i = total_elems / 4 * 4; i < total_elems; i++) { // tail case
        double value = mat_data[i];
        if (value < 0) {
            result_data[i] = value * (-1);
        } else {
            result_data[i] = value;
        }
    }
    return 0;
}



/*##########################################################################*/
// Possible Other techniques to speed up matrix Multiplication

// Memory operations are the most expensive operations, so make use of caches.

// transpose second matrix and multiply row by row


    // int i,j,k;
    // double tmp[K][J];
    // for (i = 0; i < K; ++i)
    //     for (j = 0; j < J; ++j) 
    //         tmp[i][j] = mat2->data[j][i];


    // for (i = 0; i < I; ++i) 
    //     for (j = 0; j < J; ++j)
    //         for (k = 0; k < K; ++k)
    //             result->data[i][j] += mat1->data[i][k] * tmp[j][k];

// cache blocking + simd 
// reduce the number of stores <= expensive
    // int bs;

    // if (I > 288 && J > 288) {
    //     bs = 192;
    // } else if (J > 63 && I > 63) {
    //     bs = 32;
    // } else {
    //     bs = 0;
    // }

    // int i,j,k,ii,jj,kk;
    
    // int bs = 64 / sizeof(double);
    // int bs = 128;

    // #pragma omp parallel for
    // for(jj = 0; jj < J; jj += bs) {
    //     for(ii = 0; ii < I; ii += bs) {
	// 	    for(kk = 0; kk < K; kk += bs) {
    //             for(j = jj; j < min(J, jj+bs) / 16 * 16; j += 16) {
    //                 for(i = ii; i < min(I, ii+bs); i++) {
    //                     __m256d sum = _mm256_loadu_pd(result->data_1d + i * J + j);
    //                     __m256d sum1 = _mm256_loadu_pd(result->data_1d + i * J + j + 4);
    //                     __m256d sum2 = _mm256_loadu_pd(result->data_1d + i * J + j + 8);
    //                     __m256d sum3 = _mm256_loadu_pd(result->data_1d + i * J + j + 12);
    //                     for(k = kk; k < min(K, kk+bs); k++) {
    //                         __m256d vA = _mm256_set1_pd(mat1->data_1d[i * K + k]);
    //                         __m256d vB = _mm256_loadu_pd(mat2->data_1d + k * J + j);
    //                         sum = _mm256_fmadd_pd(vA, vB, sum);

    //                         __m256d vA1 = _mm256_set1_pd(mat1->data_1d[i * K + k + 4]);
    //                         __m256d vB1 = _mm256_loadu_pd(mat2->data_1d + k * J + j + 4);
    //                         sum1 = _mm256_fmadd_pd(vA1, vB1, sum1);

    //                         __m256d vA2 = _mm256_set1_pd(mat1->data_1d[i * K + k + 8]);
    //                         __m256d vB2 = _mm256_loadu_pd(mat2->data_1d + k * J + j + 8);
    //                         sum2 = _mm256_fmadd_pd(vA2, vB2, sum2);

    //                         __m256d vA3 = _mm256_set1_pd(mat1->data_1d[i * K + k + 12]);
    //                         __m256d vB3 = _mm256_loadu_pd(mat2->data_1d + k * J + j + 12);
    //                         sum3 = _mm256_fmadd_pd(vA3, vB3, sum3);
    //                     }
    //                     _mm256_storeu_pd(result->data_1d + i * J + j, sum);
    //                     _mm256_storeu_pd(result->data_1d + i * J + j + 4, sum1);
    //                     _mm256_storeu_pd(result->data_1d + i * J + j + 8, sum2);
    //                     _mm256_storeu_pd(result->data_1d + i * J + j + 12, sum3);
    //                 }
    //             }
    //             for (j = min(J, jj+bs) / 16 * 16; j < min(J, jj+bs); j++) {
    //                 for (i = ii; i < min(I, ii+bs); i++) {
    //                     double val = 0;
    //                     for (k = kk; k < min(K, kk+bs); k++) {
    //                         val += mat1->data_1d[i * K + k] * mat2->data_1d[k * J + j];
    //                     }
    //                     result->data_1d[i * J + j] = val;
    //                 }
    //             }
    //         }
    //     }
    // }

    // // ikj or kij is optimal
    // #pragma omp parallel for
    // for(ii = 0; ii < I; ii += bs) {
	// 	for(kk = 0; kk < K; kk += bs) {
	// 		for(jj = 0; jj < J; jj += bs) {
    //             for(i = ii; i < min(I, ii+bs); i++) {
	// 				for(k = kk; k < min(K, kk+bs); k++) {
    //                     __m256d vA = _mm256_set1_pd(mat1->data_1d[i * K + k]);
	// 					for(j = jj; j < min(J, jj+bs) / 16 * 16; j += 16) {
    //                         __m256d sum = _mm256_loadu_pd(result->data_1d + i * J + j);
    //                         __m256d vB = _mm256_loadu_pd(mat2->data_1d + k * J + j);
    //                         sum = _mm256_fmadd_pd(vA, vB, sum);
    //                         _mm256_storeu_pd(result->data_1d + i * J + j, sum);

    //                         __m256d sum1 = _mm256_loadu_pd(result->data_1d + i * J + j + 4);
    //                         __m256d vB1 = _mm256_loadu_pd(mat2->data_1d + k * J + j + 4);
    //                         sum1 = _mm256_fmadd_pd(vA, vB1, sum1);
    //                         _mm256_storeu_pd(result->data_1d + i * J + j + 4, sum1);

    //                         __m256d sum2 = _mm256_loadu_pd(result->data_1d + i * J + j + 8);
    //                         __m256d vB2 = _mm256_loadu_pd(mat2->data_1d + k * J + j + 8);
    //                         sum2 = _mm256_fmadd_pd(vA, vB2, sum2);
    //                         _mm256_storeu_pd(result->data_1d + i * J + j + 8, sum2);

    //                         __m256d sum3 = _mm256_loadu_pd(result->data_1d + i * J + j + 12);
    //                         __m256d vB3 = _mm256_loadu_pd(mat2->data_1d + k * J + j + 12);
    //                         sum3 = _mm256_fmadd_pd(vA, vB3, sum3);
    //                         _mm256_storeu_pd(result->data_1d + i * J + j + 12, sum3);
    //                     }
    //                     double independent = mat1 -> data_1d[i * K + k];
    //                     for(int x = min(J, jj+bs) / 16 * 16; x < min(J, jj+bs); x++){
    //                         result->data_1d[i * J + x] += independent * mat2 -> data_1d[k * J + x];
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }


/*##########################################################################*/


// Recursive: Repeated squrare algorithm
//    if ((pow % 2) == 0) { // even power
//        // mat^(pow/2) * mat^(pow/2)
//        while (pow > 1) {
//
//        }
//
//    } else { // odd power
//        // mat^(pow/2) * mat^(pow/2) * mat
//    }

/*##########################################################################*/
