#include "numc.h"
#include <structmember.h>

PyTypeObject Matrix61cType;

/* Helper functions for initalization of matrices and vectors */

/*
 * Return a tuple given rows and cols
 */
PyObject *get_shape(int rows, int cols) {
  if (rows == 1 || cols == 1) {
    return PyTuple_Pack(1, PyLong_FromLong(rows * cols));
  } else {
    return PyTuple_Pack(2, PyLong_FromLong(rows), PyLong_FromLong(cols));
  }
}
/*
 * Matrix(rows, cols, low, high). Fill a matrix random double values
 */
int init_rand(PyObject *self, int rows, int cols, unsigned int seed, double low,
              double high) {
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    rand_matrix(new_mat, seed, low, high);
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * Matrix(rows, cols, val). Fill a matrix of dimension rows * cols with val
 */
int init_fill(PyObject *self, int rows, int cols, double val) {
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed)
        return alloc_failed;
    else {
        fill_matrix(new_mat, val);
        ((Matrix61c *)self)->mat = new_mat;
        ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    }
    return 0;
}

/*
 * Matrix(rows, cols, 1d_list). Fill a matrix with dimension rows * cols with 1d_list values
 */
int init_1d(PyObject *self, int rows, int cols, PyObject *lst) {
    if (rows * cols != PyList_Size(lst)) {
        PyErr_SetString(PyExc_ValueError, "Incorrect number of elements in list");
        return -1;
    }
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    int count = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set(new_mat, i, j, PyFloat_AsDouble(PyList_GetItem(lst, count)));
            count++;
        }
    }
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * Matrix(2d_list). Fill a matrix with dimension len(2d_list) * len(2d_list[0])
 */
int init_2d(PyObject *self, PyObject *lst) {
    int rows = PyList_Size(lst);
    if (rows == 0) {
        PyErr_SetString(PyExc_ValueError,
                        "Cannot initialize numc.Matrix with an empty list");
        return -1;
    }
    int cols;
    if (!PyList_Check(PyList_GetItem(lst, 0))) {
        PyErr_SetString(PyExc_ValueError, "List values not valid");
        return -1;
    } else {
        cols = PyList_Size(PyList_GetItem(lst, 0));
    }
    for (int i = 0; i < rows; i++) {
        if (!PyList_Check(PyList_GetItem(lst, i)) ||
                PyList_Size(PyList_GetItem(lst, i)) != cols) {
            PyErr_SetString(PyExc_ValueError, "List values not valid");
            return -1;
        }
    }
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set(new_mat, i, j,
                PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(lst, i), j)));
        }
    }
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * This deallocation function is called when reference count is 0
 */
void Matrix61c_dealloc(Matrix61c *self) {
    deallocate_matrix(self->mat);
    Py_TYPE(self)->tp_free(self);
}

/* For immutable types all initializations should take place in tp_new */
PyObject *Matrix61c_new(PyTypeObject *type, PyObject *args,
                        PyObject *kwds) {
    /* size of allocated memory is tp_basicsize + nitems*tp_itemsize*/
    Matrix61c *self = (Matrix61c *)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

/*
 * This matrix61c type is mutable, so needs init function. Return 0 on success otherwise -1
 */
int Matrix61c_init(PyObject *self, PyObject *args, PyObject *kwds) {
    /* Generate random matrices */
    if (kwds != NULL) {
        PyObject *rand = PyDict_GetItemString(kwds, "rand");
        if (!rand) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
        if (!PyBool_Check(rand)) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
        if (rand != Py_True) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }

        PyObject *low = PyDict_GetItemString(kwds, "low");
        PyObject *high = PyDict_GetItemString(kwds, "high");
        PyObject *seed = PyDict_GetItemString(kwds, "seed");
        double double_low = 0;
        double double_high = 1;
        unsigned int unsigned_seed = 0;

        if (low) {
            if (PyFloat_Check(low)) {
                double_low = PyFloat_AsDouble(low);
            } else if (PyLong_Check(low)) {
                double_low = PyLong_AsLong(low);
            }
        }

        if (high) {
            if (PyFloat_Check(high)) {
                double_high = PyFloat_AsDouble(high);
            } else if (PyLong_Check(high)) {
                double_high = PyLong_AsLong(high);
            }
        }

        if (double_low >= double_high) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }

        // Set seed if argument exists
        if (seed) {
            if (PyLong_Check(seed)) {
                unsigned_seed = PyLong_AsUnsignedLong(seed);
            }
        }

        PyObject *rows = NULL;
        PyObject *cols = NULL;
        if (PyArg_UnpackTuple(args, "args", 2, 2, &rows, &cols)) {
            if (rows && cols && PyLong_Check(rows) && PyLong_Check(cols)) {
                return init_rand(self, PyLong_AsLong(rows), PyLong_AsLong(cols), unsigned_seed, double_low,
                                 double_high);
            }
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
    }
    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    PyObject *arg3 = NULL;
    if (PyArg_UnpackTuple(args, "args", 1, 3, &arg1, &arg2, &arg3)) {
        /* arguments are (rows, cols, val) */
        if (arg1 && arg2 && arg3 && PyLong_Check(arg1) && PyLong_Check(arg2) && (PyLong_Check(arg3)
                || PyFloat_Check(arg3))) {
            if (PyLong_Check(arg3)) {
                return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), PyLong_AsLong(arg3));
            } else
                return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), PyFloat_AsDouble(arg3));
        } else if (arg1 && arg2 && arg3 && PyLong_Check(arg1) && PyLong_Check(arg2) && PyList_Check(arg3)) {
            /* Matrix(rows, cols, 1D list) */
            return init_1d(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), arg3);
        } else if (arg1 && PyList_Check(arg1) && arg2 == NULL && arg3 == NULL) {
            /* Matrix(rows, cols, 1D list) */
            return init_2d(self, arg1);
        } else if (arg1 && arg2 && PyLong_Check(arg1) && PyLong_Check(arg2) && arg3 == NULL) {
            /* Matrix(rows, cols, 1D list) */
            return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), 0);
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return -1;
    }
}

/*
 * List of lists representations for matrices
 */
PyObject *Matrix61c_to_list(Matrix61c *self) {
    int rows = self->mat->rows;
    int cols = self->mat->cols;
    PyObject *py_lst = NULL;
    if (self->mat->is_1d) {  // If 1D matrix, print as a single list
        py_lst = PyList_New(rows * cols);
        int count = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                PyList_SetItem(py_lst, count, PyFloat_FromDouble(get(self->mat, i, j)));
                count++;
            }
        }
    } else {  // if 2D, print as nested list
        py_lst = PyList_New(rows);
        for (int i = 0; i < rows; i++) {
            PyList_SetItem(py_lst, i, PyList_New(cols));
            PyObject *curr_row = PyList_GetItem(py_lst, i);
            for (int j = 0; j < cols; j++) {
                PyList_SetItem(curr_row, j, PyFloat_FromDouble(get(self->mat, i, j)));
            }
        }
    }
    return py_lst;
}

PyObject *Matrix61c_class_to_list(Matrix61c *self, PyObject *args) {
    PyObject *mat = NULL;
    if (PyArg_UnpackTuple(args, "args", 1, 1, &mat)) {
        if (!PyObject_TypeCheck(mat, &Matrix61cType)) {
            PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
            return NULL;
        }
        Matrix61c* mat61c = (Matrix61c*)mat;
        return Matrix61c_to_list(mat61c);
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return NULL;
    }
}

/*
 * Add class methods
 */
PyMethodDef Matrix61c_class_methods[] = {
    {"to_list", (PyCFunction)Matrix61c_class_to_list, METH_VARARGS, "Returns a list representation of numc.Matrix"},
    {NULL, NULL, 0, NULL}
};

/*
 * Matrix61c string representation. For printing purposes.
 */
PyObject *Matrix61c_repr(PyObject *self) {
    PyObject *py_lst = Matrix61c_to_list((Matrix61c *)self);
    return PyObject_Repr(py_lst);
}

/* NUMBER METHODS */

/*
 * Add the second numc.Matrix (Matrix61c) object to the first one. The first operand is
 * self, and the second operand can be obtained by casting `args`.
 */
PyObject *Matrix61c_add(Matrix61c* self, PyObject* args) {

    if (!PyObject_TypeCheck(args, &Matrix61cType)) {
        PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
        return NULL;
    }
    Matrix61c *other = (Matrix61c *) args;

    Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);

    matrix *result;

    int alloc_failed = allocate_matrix(&result, self->mat->rows, self->mat->cols);
    if (alloc_failed) {
        Matrix61c_dealloc(rv);
        return NULL;
    }

    int operation_failed = add_matrix(result, self->mat, other->mat);
    if (operation_failed) {
        deallocate_matrix(result);
        Matrix61c_dealloc(rv);
        return NULL;
    }

    rv->mat = result;
    rv->shape = get_shape(result->rows, result->cols);

    return (PyObject *) rv;
}

/*
 * Subtract the second numc.Matrix (Matrix61c) object from the first one. The first operand is
 * self, and the second operand can be obtained by casting `args`.
 */
PyObject *Matrix61c_sub(Matrix61c* self, PyObject* args) {
    if (!PyObject_TypeCheck(args, &Matrix61cType)) {
        PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
        return NULL;
    }
    Matrix61c *other = (Matrix61c *) args;

    Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);

    matrix *result;

    int alloc_failed = allocate_matrix(&result, self->mat->rows, other->mat->cols);
    if (alloc_failed) {
        Matrix61c_dealloc(rv);
        return NULL;
    }

    int operation_failed = sub_matrix(result, self->mat, other->mat);
    if (operation_failed) {
        Matrix61c_dealloc(rv);
        return NULL;
    }

    rv->mat = result;
    rv->shape = get_shape(result->rows, result->cols);

    return (PyObject *) rv;
}

/*
 * NOT element-wise multiplication. The first operand is self, and the second operand
 * can be obtained by casting `args`.
 */
PyObject *Matrix61c_multiply(Matrix61c* self, PyObject *args) {
    if (!PyObject_TypeCheck(args, &Matrix61cType)) {
        PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
        return NULL;
    }
    Matrix61c *other = (Matrix61c *) args;

    Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);

    matrix *result;

    int alloc_failed = allocate_matrix(&result, self->mat->rows, other->mat->cols);
    if (alloc_failed) {
        Matrix61c_dealloc(rv);
        return NULL;
    }

    int operation_failed = mul_matrix(result, self->mat, other->mat);
    if (operation_failed) {
        Matrix61c_dealloc(rv);
        return NULL;
    }

    rv->mat = result;
    rv->shape = get_shape(result->rows, result->cols);

    return (PyObject *) rv;
}

/*
 * Negates the given numc.Matrix.
 */
PyObject *Matrix61c_neg(Matrix61c* self) {

    Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);

    matrix *result;

    int alloc_failed = allocate_matrix(&result, self->mat->rows, self->mat->cols);
    if (alloc_failed) {
        Matrix61c_dealloc(rv);
        return NULL;
    }

    int operation_failed = neg_matrix(result, self->mat);
    if (operation_failed) {
        Matrix61c_dealloc(rv);
        return NULL;
    }

    rv->mat = result;
    rv->shape = get_shape(result->rows, result->cols);

    return (PyObject *) rv;
}

/*
 * Take the element-wise absolute value of this numc.Matrix.
 */
PyObject *Matrix61c_abs(Matrix61c *self) {

    Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);

    matrix *result;

    int alloc_failed = allocate_matrix(&result, self->mat->rows, self->mat->cols);
    if (alloc_failed) {
        Matrix61c_dealloc(rv);
        return NULL;
    }

    int operation_failed = abs_matrix(result, self->mat);
    if (operation_failed) {
        Matrix61c_dealloc(rv);
        return NULL;
    }

    rv->mat = result;
    rv->shape = get_shape(result->rows, result->cols);

    return (PyObject *) rv;
}

/*
 * Raise numc.Matrix (Matrix61c) to the `pow`th power. You can ignore the argument `optional`.
 */
PyObject *Matrix61c_pow(Matrix61c *self, PyObject *pow, PyObject *optional) {

    if (!PyLong_Check(pow)) {
        PyErr_SetString(PyExc_TypeError, "Pow must be an integer");
        return NULL;
    }

    Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);

    matrix *result;

    int alloc_failed = allocate_matrix(&result, self->mat->rows, self->mat->cols);
    if (alloc_failed) {
        Matrix61c_dealloc(rv);
        return NULL;
    }

    int operation_failed = pow_matrix(result, self->mat, PyLong_AsLong(pow));
    if (operation_failed) {
        Matrix61c_dealloc(rv);
        return NULL;
    }

    rv->mat = result;
    rv->shape = get_shape(result->rows, result->cols);

    return (PyObject *) rv;
}

/*
 * Create a PyNumberMethods struct for overloading operators with all the number methods you have
 * define. You might find this link helpful: https://docs.python.org/3.6/c-api/typeobj.html
 */
PyNumberMethods Matrix61c_as_number = {
        .nb_add = (binaryfunc) Matrix61c_add,
        .nb_subtract = (binaryfunc) Matrix61c_sub,
        .nb_multiply = (binaryfunc) Matrix61c_multiply,
        .nb_negative = (unaryfunc) Matrix61c_neg,
        .nb_absolute = (unaryfunc) Matrix61c_abs,
        .nb_power = (ternaryfunc) Matrix61c_pow
};


/* INSTANCE METHODS */

/*
 * Given a numc.Matrix self, parse `args` to (int) row, (int) col, and (double/int) val.
 * Return None in Python (this is different from returning null).
 */
PyObject *Matrix61c_set_value(Matrix61c *self, PyObject* args) {
    int row = 0;
    int col = 0;
    double val = 0;

    if (PyTuple_Size(args) != 3) {
        PyErr_SetString(PyExc_TypeError, "Invalid number of arguments");
        return NULL;
    }

    PyObject *first = PyTuple_GetItem(args, 0);
    if (first != NULL) {
        if (PyLong_Check(first)) {
            row = PyLong_AsLong(first);
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments. 1");
            return NULL;
        }
    }
    PyObject *second = PyTuple_GetItem(args, 1);
    if (second != NULL) {
        if (PyLong_Check(second)) {
            col = PyLong_AsLong(second);
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments. 2");
            return NULL;
        }
    }

    PyObject *third = PyTuple_GetItem(args, 2);
    if (third != NULL) {
        if (PyLong_Check(third)) {
            val = PyLong_AsDouble(third);
        } else if (PyFloat_Check(third)) {
            val = PyFloat_AsDouble(third);
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments 3");
            return NULL;
        }
    }

    if (row >= self->mat->rows || col >= self->mat->cols) {
        PyErr_SetString(PyExc_IndexError, "My error: Index out of range");
        return NULL;
    }

    set(self->mat, row, col, val);

    return Py_None;
}

/*
 * Given a numc.Matrix `self`, parse `args` to (int) row and (int) col.
 * Return the value at the `row`th row and `col`th column, which is a Python
 * float/int.
 */
PyObject *Matrix61c_get_value(Matrix61c *self, PyObject* args) {
    long row = 0;
    long col = 0;

    if (PyTuple_Size(args) != 2) {
        PyErr_SetString(PyExc_TypeError, "Invalid number of arguments");
        return NULL;
    }

    PyObject *first = PyTuple_GetItem(args, 0);
    if (first != NULL) {
        if (PyLong_Check(first)) {
            row = PyLong_AsLong(first);
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments 4");
            return NULL;
        }
    }
    PyObject *second = PyTuple_GetItem(args, 1);
    if (second != NULL) {
        if (PyLong_Check(second)) {
            col = PyLong_AsLong(second);
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments 5");
            return NULL;
        }
    }

    if (row >= self->mat->rows || col >= self->mat->cols) {
        PyErr_SetString(PyExc_IndexError, "My error: Index out of range");
        return NULL;
    }

    return PyFloat_FromDouble(get(self->mat, row, col));
    // we are not expected to support for 1d matrices
}

/*
 * Create an array of PyMethodDef structs to hold the instance methods.
 * Name the python function corresponding to Matrix61c_get_value as "get" and Matrix61c_set_value
 * as "set"
 * You might find this link helpful: https://docs.python.org/3.6/c-api/structures.html
 */
PyMethodDef Matrix61c_methods[] = {
    {"set", (PyCFunction)&Matrix61c_set_value, METH_VARARGS, "Sets self's key equal to value"},
    {"get", (PyCFunction)&Matrix61c_get_value, METH_VARARGS, "Gets the self's value at key index"},
    {NULL, NULL, 0, NULL}
};

/* INDEXING */

/*
 * Given a numc.Matrix `self`, index into it with `key`. Return the indexed result.
 */
PyObject *Matrix61c_subscript(Matrix61c* self, PyObject* key) {

    long row_index, col_index,
        row_start_index, row_end_index, row_step_index,
        col_start_index, col_end_index, col_step_index,
        row_length, col_length, row_slice_length, col_slice_length;

    matrix *result;

    Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);

    int allocation_failed;

    if (self->mat->rows == 1 || self->mat->cols == 1) { // if the matrix is 1d
        if (PyLong_Check(key)) { // int ->test passed
            col_length = self->mat->rows * self->mat->cols;
            col_index = PyLong_AsLong(key);
            if (col_index >= col_length) {
                PyErr_SetString(PyExc_IndexError, "My error: Key is out of range");
                return NULL;
            }
            // this is right. 1 by 1
            if (self->mat->rows == 1) {
                return PyFloat_FromDouble(self->mat->data[0][col_index]);
            } else {
                return PyFloat_FromDouble(self->mat->data[col_index][0]);
            }

        } else if (PySlice_Check(key)) { // slice -> test passed
            col_length = self->mat->rows * self->mat->cols;
            PySlice_GetIndicesEx(key, col_length, &col_start_index, &col_end_index, &col_step_index, &col_slice_length);
            if (col_step_index != 1) {
                PyErr_SetString(PyExc_ValueError, "My error: Numc only supports step size == 1");
                return NULL;
            }
            if (col_slice_length < 1) {
                PyErr_SetString(PyExc_ValueError, "My error: Slice must be of length >= 1");
                return NULL;
            }
            if (col_slice_length == 1) { // result is 1 by 1
                // this is right
                if (self->mat->rows == 1) {
                    return PyFloat_FromDouble(self->mat->data[0][col_start_index]);
                } else {
                    return PyFloat_FromDouble(self->mat->data[col_start_index][0]);
                }
            } else { // result is not 1 by 1
                if (self->mat->rows == 1) {
                    allocation_failed = allocate_matrix_ref(&result, self->mat, 0, col_start_index, 1, col_slice_length);
                } else {
                    allocation_failed = allocate_matrix_ref(&result, self->mat, col_start_index, 0, col_slice_length, 1);
                }
                rv->shape = get_shape(1, col_slice_length);
            }

        } else if (PyTuple_Check(key)) { // tuple -> test passed
            PyErr_SetString(PyExc_TypeError, "My error: 1D matrices only support single slice");
            return NULL;
        } else {
            PyErr_SetString(PyExc_TypeError, "My error: Invalid arguments");
            return NULL;
        }
    } else { // if the matrix is 2d
        if (PyLong_Check(key)) { // int -> test passed
            row_index = PyLong_AsLong(key);
            if (row_index >= self->mat->rows) {
                PyErr_SetString(PyExc_IndexError, "My error: Key is out of range");
                return NULL;
            }
            // this is right
            allocation_failed = allocate_matrix_ref(&result, self->mat, row_index, 0, 1, self->mat->cols);
            rv->shape = get_shape(1, self->mat->cols);
        } else if (PySlice_Check(key)) { // slice -> test passed
            row_length = self->mat->rows;
            PySlice_GetIndicesEx(key, row_length, &row_start_index, &row_end_index, &row_step_index, &row_slice_length);
            if (row_slice_length < 1) {
                PyErr_SetString(PyExc_ValueError, "My error: Slice must be of length >= 1");
                return NULL;
            }
            if (row_step_index != 1) {
                PyErr_SetString(PyExc_ValueError, "My error: Numc only supports step size == 1");
                return NULL;
            }
            // this is right
            allocation_failed = allocate_matrix_ref(&result, self->mat, row_start_index, 0, row_slice_length, self->mat->cols);
            rv->shape = get_shape(row_slice_length, self->mat->cols);
        } else if (PyTuple_Check(key)) { // () <- tuple
            PyObject *first = NULL;
            PyObject *second = NULL;
            PyArg_UnpackTuple(key, "args", 2, 2, &first, &second);
            if (PyLong_Check(first) && PyLong_Check(second)) { // (int, int)
                row_index = PyLong_AsLong(first);
                col_index = PyLong_AsLong(second);
                if ((row_index >= self->mat->rows) || (col_index >= self->mat->cols)) {
                    PyErr_SetString(PyExc_IndexError, "My error: Key is out of range");
                    return NULL;
                }
                // this is right. 1 by 1
                return Matrix61c_get_value(self, PyTuple_Pack(2, PyLong_FromLong(row_index), PyLong_FromLong(col_index)));
            } else if (PyLong_Check(first) && PySlice_Check(second)) { // (int, slice)
                row_index = PyLong_AsLong(first);
                if (row_index >= self->mat->rows) {
                    PyErr_SetString(PyExc_IndexError, "My error: Key is out of range");
                    return NULL;
                }
                col_length = self->mat->cols;
                PySlice_GetIndicesEx(second, col_length, &col_start_index, &col_end_index, &col_step_index, &col_slice_length);
                if (col_step_index != 1) {
                    PyErr_SetString(PyExc_ValueError, "My error: Numc only supports step size == 1");
                    return NULL;
                }
                if (col_slice_length < 1) {
                    PyErr_SetString(PyExc_ValueError, "My error: Slice must be of length >= 1");
                    return NULL;
                }
                if (col_slice_length == 1) { // result is 1 by 1
                    // this is right
                    return Matrix61c_get_value(self, PyTuple_Pack(2, PyLong_FromLong(row_index), PyLong_FromLong(col_start_index)));
                } else { // result is not 1 by 1
                    // this is right
                    allocation_failed = allocate_matrix_ref(&result, self->mat, row_index, col_start_index, 1, col_slice_length);
                    rv->shape = get_shape(1, col_slice_length);
                }

            } else if (PySlice_Check(first) && PyLong_Check(second)) { // (slice, int)
                row_length = self->mat->rows;
                PySlice_GetIndicesEx(first, row_length, &row_start_index, &row_end_index, &row_step_index, &row_slice_length);
                if (row_step_index != 1) {
                    PyErr_SetString(PyExc_ValueError, "My error: Numc only supports step size == 1");
                    return NULL;
                }
                col_index = PyLong_AsLong(second);
                if (col_index >= self->mat->cols) {
                    PyErr_SetString(PyExc_IndexError, "My error: Key is out of range");
                    return NULL;
                }
                if (row_slice_length < 1) {
                    PyErr_SetString(PyExc_ValueError, "My error: Slice must be of length >= 1");
                    return NULL;
                }
                if (row_slice_length == 1) { // result is 1 by 1
                    // this is right
                    return Matrix61c_get_value(self, PyTuple_Pack(2, PyLong_FromLong(row_start_index), PyLong_FromLong(col_index)));
                } else { // result is not 1 by 1
                    // this is right
                    allocation_failed = allocate_matrix_ref(&result, self->mat, row_start_index, col_index, row_slice_length, 1);
                    rv->shape = get_shape(row_slice_length, 1);
                }
            } else if (PySlice_Check(first) && PySlice_Check(second)) { // (slice, slice)
                row_length = self->mat->rows;
                col_length = self->mat->cols;
                PySlice_GetIndicesEx(first, row_length, &row_start_index, &row_end_index, &row_step_index, &row_slice_length);
                PySlice_GetIndicesEx(second, col_length, &col_start_index, &col_end_index, &col_step_index, &col_slice_length);
                if ((row_step_index != 1) || (col_step_index != 1)) {
                    PyErr_SetString(PyExc_ValueError, "My error: Numc only supports step size == 1");
                    return NULL;
                }
                if ((row_slice_length < 1) || (col_slice_length < 1)) {
                    PyErr_SetString(PyExc_ValueError, "My error: Slice must be of length >= 1");
                    return NULL;
                }
                if (row_slice_length == 1 && col_slice_length == 1) { // result is 1 by 1
                    // this is right
                    return Matrix61c_get_value(self, PyTuple_Pack(2, PyLong_FromLong(row_start_index), PyLong_FromLong(col_start_index)));
                } else { // result is not 1 by 1
                    // this is right
                    allocation_failed = allocate_matrix_ref(&result, self->mat, row_start_index, col_start_index, row_slice_length, col_slice_length);
                    rv->shape = get_shape(row_slice_length, col_slice_length);
                }
            } else {
                PyErr_SetString(PyExc_TypeError, "My error: Invalid arguments");
                return NULL;
            }
        } else {
            PyErr_SetString(PyExc_TypeError, "My error: Invalid arguments");
            return NULL;
        }
    }

    if (allocation_failed) {
        Matrix61c_dealloc(rv);
        return NULL;
    }

    rv->mat = result;

    return (PyObject *) rv;
}

static int set_1d_slice_to_v(Matrix61c* self, int row_start_index, int col_start_index, int col_slice_length, int int_slice, PyObject *v) {
    int col_count = 0;
    if (PyList_Size(v) == col_slice_length) {
        for (int i = col_start_index; i < (col_start_index + col_slice_length); i++) {
            PyObject *elem = PyList_GetItem(v, col_count); // we don't want the item of v at i!!!! Initialize some counter or something
            col_count++;
            if (PyFloat_Check(elem) || PyLong_Check(elem)) {
                // row = row_start_index, col = i, val = elem
                PyObject *set_failed;
                if (self->mat->rows == 1) {
                    if (PyFloat_Check(elem)) {
                        set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(row_start_index), PyLong_FromLong(i), PyFloat_FromDouble(PyFloat_AsDouble(elem))));
                    } else {
                        set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(row_start_index), PyLong_FromLong(i), PyLong_FromDouble(PyLong_AsDouble(elem))));
                    }
                } else if (self->mat->cols == 1) {
                    if (PyFloat_Check(elem)) {
                        set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(i), PyLong_FromLong(row_start_index), PyFloat_FromDouble(PyFloat_AsDouble(elem))));
                    } else {
                        set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(i), PyLong_FromLong(row_start_index), PyLong_FromDouble(PyLong_AsDouble(elem))));
                    }
                } else {
                    if (int_slice == 1) { // int_slice
                        if (PyFloat_Check(elem)) {
                            set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(row_start_index), PyLong_FromLong(i), PyFloat_FromDouble(PyFloat_AsDouble(elem))));
                        } else {
                            set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(row_start_index), PyLong_FromLong(i), PyLong_FromDouble(PyLong_AsDouble(elem))));
                        }
                    } else if (int_slice == 0) { // slice_int
                        if (PyFloat_Check(elem)) {
                            set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(i), PyLong_FromLong(row_start_index), PyFloat_FromDouble(PyFloat_AsDouble(elem))));
                        } else {
                            set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(i), PyLong_FromLong(row_start_index), PyLong_FromDouble(PyLong_AsDouble(elem))));
                        }
                    } else { // 2d -> int
                        if (PyFloat_Check(elem)) {
                            set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(row_start_index), PyLong_FromLong(i), PyFloat_FromDouble(PyFloat_AsDouble(elem))));
                        } else {
                            set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(row_start_index), PyLong_FromLong(i), PyLong_FromDouble(PyLong_AsDouble(elem))));
                        }
                    }
                }
                if (set_failed == NULL) return -1;
            } else {
                PyErr_SetString(PyExc_ValueError, "My error: Resulting slice is 1D, but element of v is not a float or int");
                return -1;
            }
        }
        return 0;
    } else {
        PyErr_SetString(PyExc_ValueError, "My error: Resulting slice is 1D, but v has the wrong length");
        return -1;
    }
}


static int set_2d_slice_to_v(Matrix61c* self, int row_start_index, int col_start_index, int row_slice_length, int col_slice_length, PyObject *v) {
    int row_count = 0;
    if (PyList_Size(v) == row_slice_length) {
        for (int i = row_start_index; i < (row_slice_length + row_start_index); i++) {
            PyObject *elem = PyList_GetItem(v, row_count);
            row_count++;
            int col_count = 0;
            if (PyList_Check(elem)) {
                if (PyList_Size(elem) == col_slice_length) {
                    for (int j = col_start_index; j < ( col_start_index + col_slice_length); j++) {
                        PyObject *sub_elem = PyList_GetItem(elem, col_count);
                        col_count++;
                        if (PyFloat_Check(sub_elem) || PyLong_Check(sub_elem)) {
                            // row = i, col = j, val = sub_elem
                            PyObject *set_failed;
                            if (PyFloat_Check(sub_elem)) {
                                set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(i), PyLong_FromLong(j), PyFloat_FromDouble(PyFloat_AsDouble(sub_elem))));
                            } else {
                                set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(i), PyLong_FromLong(j), PyLong_FromDouble(PyLong_AsDouble(sub_elem))));
                            }
                            if (set_failed == NULL) return -1;
                        } else {
                            PyErr_SetString(PyExc_ValueError, "My error: Resulting slice is 2D, but element of an element of v is not a float or int");
                            return -1;
                        }
                    }
                } else {
                    PyErr_SetString(PyExc_ValueError, "My error: Resulting slice is 2D, but an element of v has the wrong length");
                    return -1;
                }
            } else {
                PyErr_SetString(PyExc_TypeError, "My error: Resulting slice is 2D, but element of v is not list. Wasn't mentioned in spec. I created this. Check back.");
                return -1;
            }
        }
        return 0;
    } else {
        PyErr_SetString(PyExc_ValueError, "My error: Resulting slice is 2D, but v has the wrong length"); // not necessarily 2d
        return -1;
    }
}

/*
 * Given a numc.Matrix `self`, index into it with `key`, and set the indexed result to `v`.
 */
int Matrix61c_set_subscript(Matrix61c* self, PyObject *key, PyObject *v) {

    long row_index, col_index,
            row_start_index, row_end_index, row_step_index,
            col_start_index, col_end_index, col_step_index,
            row_length, col_length, row_slice_length, col_slice_length;

    if (self->mat->rows == 1 || self->mat->cols == 1) { // if the matrix is 1d
        if (PyLong_Check(key)) { // int
            col_length = self->mat->rows * self->mat->cols;
            col_index = PyLong_AsLong(key);
            if (col_index >= col_length) {
                PyErr_SetString(PyExc_IndexError, "My error: Key is out of range");
                return -1;
            }
            // Resulting slice is 1 by 1, but v is not a float or int
            if (PyFloat_Check(v) || PyLong_Check(v)) {
                if (PyFloat_Check(v)) {
                    if (self->mat->rows == 1) {
                        self->mat->data[0][col_index] = PyFloat_AsDouble(v);
                    } else {
                        self->mat->data[col_index][0] = PyFloat_AsDouble(v);
                    }
                } else {
                    if (self->mat->rows == 1) {
                        self->mat->data[0][col_index] = PyLong_AsDouble(v);
                    } else {
                        self->mat->data[col_index][0] = PyLong_AsDouble(v);
                    }
                }
                return 0;
            } else {
                PyErr_SetString(PyExc_TypeError, "My error: Resulting slice is 1 by 1, but v is not a float or int");
                return -1;
            }
        } else if (PySlice_Check(key)) { // slice
            col_length = self->mat->rows * self->mat->cols;
            PySlice_GetIndicesEx(key, col_length, &col_start_index, &col_end_index, &col_step_index, &col_slice_length);
            if (col_step_index != 1) {
                PyErr_SetString(PyExc_ValueError, "My error: Numc only supports step size == 1");
                return -1;
            }
            if (col_slice_length < 1) {
                PyErr_SetString(PyExc_ValueError, "My error: Slice must be of length >= 1");
                return -1;
            }
            if (col_slice_length == 1) { // result is 1 by 1
                if (PyFloat_Check(v) || PyLong_Check(v)) {
                    if (PyFloat_Check(v)) {
                        if (self->mat->rows == 1) {
                            self->mat->data[0][col_start_index] = PyFloat_AsDouble(v);
                        } else {
                            self->mat->data[col_start_index][0] = PyFloat_AsDouble(v);
                        }
                    } else {
                        if (self->mat->rows == 1) {
                            self->mat->data[0][col_start_index] = PyLong_AsDouble(v);
                        } else {
                            self->mat->data[col_start_index][0] = PyLong_AsDouble(v);
                        }
                    }
                    return 0;
                } else {
                    PyErr_SetString(PyExc_TypeError, "My error: Resulting slice is 1 by 1, but v is not a float or int");
                    return -1;
                }
            } else { // result is not 1 by 1
                // Resulting slice is not 1 by 1, but v is not a list.
                if (PyList_Check(v)) {
                    int setting_v_failed;
                    setting_v_failed = set_1d_slice_to_v(self, 0, col_start_index, col_slice_length, 1000, v);
                     if (setting_v_failed) return -1;
                } else {
                    PyErr_SetString(PyExc_TypeError, "My error: Resulting slice is not 1 by 1, but v is not a list.");
                    return -1;
                }
            }
        } else if (PyTuple_Check(key)) { // tuple
            PyErr_SetString(PyExc_TypeError, "My error: 1D matrices only support single slice");
            return -1;
        } else {
            PyErr_SetString(PyExc_TypeError, "My error: Invalid arguments");
            return -1;
        }
    } else { // if the matrix is 2d
        if (PyLong_Check(key)) { // int
            row_index = PyLong_AsLong(key);
            if (row_index >= self->mat->rows) {
                PyErr_SetString(PyExc_IndexError, "My error: Key is out of range");
                return -1;
            }
            // Resulting slice is not 1 by 1, but v is not a list.
            if (PyList_Check(v)) {
                int setting_v_failed = set_1d_slice_to_v(self, row_index, 0, self->mat->cols, 1000, v);
                if (setting_v_failed) return -1;
            } else {
                PyErr_SetString(PyExc_TypeError, "My error: Resulting slice is not 1 by 1, but v is not a list.");
                return -1;
            }
        } else if (PySlice_Check(key)) { // slice
            row_length = self->mat->rows;
            PySlice_GetIndicesEx(key, row_length, &row_start_index, &row_end_index, &row_step_index, &row_slice_length);
            if (row_slice_length < 1) {
                PyErr_SetString(PyExc_ValueError, "My error: Slice must be of length >= 1");
                return -1;
            }
            if (row_step_index != 1) {
                PyErr_SetString(PyExc_ValueError, "My error: Numc only supports step size == 1");
                return -1;
            }
            // might have to add checks for where row_slice_length == 1. => if so, copy from 1d slice
            // Resulting slice is not 1 by 1, but v is not a list.
            if (PyList_Check(v)) {
                int setting_v_failed = set_2d_slice_to_v(self, row_start_index, 0, row_slice_length, self->mat->cols, v);
                if (setting_v_failed) return -1;
            } else {
                PyErr_SetString(PyExc_TypeError, "My error: Resulting slice is not 1 by 1, but v is not a list.");
                return -1;
            }
        } else if (PyTuple_Check(key)) { // () <- tuple
            PyObject *first = NULL;
            PyObject *second = NULL;
            PyArg_UnpackTuple(key, "args", 2, 2, &first, &second);
            if (PyLong_Check(first) && PyLong_Check(second)) { // (int, int)
                row_index = PyLong_AsLong(first);
                col_index = PyLong_AsLong(second);
                if ((row_index >= self->mat->rows) || (col_index >= self->mat->cols)) {
                    PyErr_SetString(PyExc_IndexError, "My error: Key is out of range");
                    return -1;
                }
                // Resulting slice is 1 by 1, but v is not a float or int
                if (PyFloat_Check(v) || PyLong_Check(v)) {
                    PyObject *set_failed;
                    if (PyFloat_Check(v)) {
                        set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(row_index), PyLong_FromLong(col_index), PyFloat_FromDouble(PyFloat_AsDouble(v))));
                    } else {
                        set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(row_index), PyLong_FromLong(col_index), PyLong_FromDouble(PyLong_AsDouble(v))));
                    }
                    if (set_failed == NULL) return -1;
                    return 0;
                } else {
                    PyErr_SetString(PyExc_TypeError, "My error: Resulting slice is 1 by 1, but v is not a float or int");
                    return -1;
                }
            } else if (PyLong_Check(first) && PySlice_Check(second)) { // (int, slice)
                row_index = PyLong_AsLong(first);
                if (row_index >= self->mat->rows) {
                    PyErr_SetString(PyExc_IndexError, "My error: Key is out of range");
                    return -1;
                }
                col_length = self->mat->cols;
                PySlice_GetIndicesEx(second, col_length, &col_start_index, &col_end_index, &col_step_index, &col_slice_length);
                if (col_step_index != 1) {
                    PyErr_SetString(PyExc_ValueError, "My error: Numc only supports step size == 1");
                    return -1;
                }
                if (col_slice_length < 1) {
                    PyErr_SetString(PyExc_ValueError, "My error: Slice must be of length >= 1");
                    return -1;
                }
                if (col_slice_length == 1) { // result is 1 by 1
                    // Resulting slice is 1 by 1, but v is not a float or int
                    if (PyFloat_Check(v) || PyLong_Check(v)) {
                        // row = row_index, col = col_start_index, val = v
                        PyObject *set_failed;
                        if (PyFloat_Check(v)) {
                            set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(row_index), PyLong_FromLong(col_start_index), PyFloat_FromDouble(PyFloat_AsDouble(v))));
                        } else {
                            set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(row_index), PyLong_FromLong(col_start_index), PyLong_FromDouble(PyLong_AsDouble(v))));
                        }
                        if (set_failed == NULL) return -1;
                        return 0;
                    } else {
                        PyErr_SetString(PyExc_TypeError, "My error: Resulting slice is 1 by 1, but v is not a float or int");
                        return -1;
                    }
                } else { // result is not 1 by 1
                    // Resulting slice is not 1 by 1, but v is not a list.
                    if (PyList_Check(v)) {
                        int setting_v_failed = set_1d_slice_to_v(self, row_index, col_start_index, col_slice_length, 1, v);
                        if (setting_v_failed) return -1;
                    } else {
                        PyErr_SetString(PyExc_TypeError, "My error: Resulting slice is not 1 by 1, but v is not a list.");
                        return -1;
                    }
                }

            } else if (PySlice_Check(first) && PyLong_Check(second)) { // (slice, int)
                row_length = self->mat->rows;
                PySlice_GetIndicesEx(first, row_length, &row_start_index, &row_end_index, &row_step_index, &row_slice_length);
                if (row_step_index != 1) {
                    PyErr_SetString(PyExc_ValueError, "My error: Numc only supports step size == 1");
                    return -1;
                }
                col_index = PyLong_AsLong(second);
                if (col_index >= self->mat->cols) {
                    PyErr_SetString(PyExc_IndexError, "My error: Key is out of range");
                    return -1;
                }
                if (row_slice_length < 1) {
                    PyErr_SetString(PyExc_ValueError, "My error: Slice must be of length >= 1");
                    return -1;
                }
                if (row_slice_length == 1) { // result is 1 by 1
                    // Resulting slice is 1 by 1, but v is not a float or int
                    if (PyFloat_Check(v) || PyLong_Check(v)) {
                        // row = row_start_index, col = col_index, val = v
                        PyObject *set_failed;
                        if (PyFloat_Check(v)) {
                            set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(row_start_index), PyLong_FromLong(col_index), PyFloat_FromDouble(PyFloat_AsDouble(v))));
                        } else {
                            set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(row_start_index), PyLong_FromLong(col_index), PyLong_FromDouble(PyLong_AsDouble(v))));
                        }
                        if (set_failed == NULL) return -1;
                        return 0;
                    } else {
                        PyErr_SetString(PyExc_TypeError, "My error: Resulting slice is 1 by 1, but v is not a float or int");
                        return -1;
                    }
                } else { // result is not 1 by 1
                    // Resulting slice is not 1 by 1, but v is not a list.
                    if (PyList_Check(v)) {
                        int setting_v_failed = set_1d_slice_to_v(self, col_index, row_start_index, row_slice_length, 0, v);
                        if (setting_v_failed) return -1;
                    } else {
                        PyErr_SetString(PyExc_TypeError, "My error: Resulting slice is not 1 by 1, but v is not a list.");
                        return -1;
                    }

                }
            } else if (PySlice_Check(first) && PySlice_Check(second)) { // (slice, slice)
                row_length = self->mat->rows;
                col_length = self->mat->cols;
                PySlice_GetIndicesEx(first, row_length, &row_start_index, &row_end_index, &row_step_index, &row_slice_length);
                PySlice_GetIndicesEx(second, col_length, &col_start_index, &col_end_index, &col_step_index, &col_slice_length);
                if ((row_step_index != 1) || (col_step_index != 1)) {
                    PyErr_SetString(PyExc_ValueError, "My error: Numc only supports step size == 1");
                    return -1;
                }
                if ((row_slice_length < 1) || (col_slice_length < 1)) {
                    PyErr_SetString(PyExc_ValueError, "My error: Slice must be of length >= 1");
                    return -1;
                }
                if (row_slice_length == 1 && col_slice_length == 1) { // result is 1 by 1
                    // Resulting slice is 1 by 1, but v is not a float or int
                    if (PyFloat_Check(v) || PyLong_Check(v)) {
                        // row = row_start_index, col = col_start_index, val = v
                        PyObject *set_failed;
                        if (PyFloat_Check(v)) {
                            set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(row_start_index), PyLong_FromLong(col_start_index), PyFloat_FromDouble(PyFloat_AsDouble(v))));
                        } else {
                            set_failed = Matrix61c_set_value(self, PyTuple_Pack(3, PyLong_FromLong(row_start_index), PyLong_FromLong(col_start_index), PyLong_FromDouble(PyLong_AsDouble(v))));
                        }
                        if (set_failed == NULL) return -1;
                        return 0;
                    } else {
                        PyErr_SetString(PyExc_TypeError, "My error: Resulting slice is 1 by 1, but v is not a float or int");
                        return -1;
                    }
                } else { // result is not 1 by 1
                    // Resulting slice is not 1 by 1, but v is not a list.
                    if (PyList_Check(v)) {
                        int setting_v_failed = set_2d_slice_to_v(self, row_start_index, col_start_index, row_slice_length, col_slice_length, v);
                        if (setting_v_failed) return -1;
                    } else {
                        PyErr_SetString(PyExc_TypeError, "My error: Resulting slice is not 1 by 1, but v is not a list.");
                        return -1;
                    }
                }
            } else {
                PyErr_SetString(PyExc_TypeError, "My error: Invalid arguments");
                return -1;
            }
        } else {
            PyErr_SetString(PyExc_TypeError, "My error: Invalid arguments");
            return -1;
        }
    }
    return 0;
}

PyMappingMethods Matrix61c_mapping = {
    NULL,
    (binaryfunc) Matrix61c_subscript,
    (objobjargproc) Matrix61c_set_subscript,
};

/* INSTANCE ATTRIBUTES*/
PyMemberDef Matrix61c_members[] = {
    {
        "shape", T_OBJECT_EX, offsetof(Matrix61c, shape), 0,
        "(rows, cols)"
    },
    {NULL}  /* Sentinel */
};

PyTypeObject Matrix61cType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "numc.Matrix",
    .tp_basicsize = sizeof(Matrix61c),
    .tp_dealloc = (destructor)Matrix61c_dealloc,
    .tp_repr = (reprfunc)Matrix61c_repr,
    .tp_as_number = &Matrix61c_as_number,
    .tp_flags = Py_TPFLAGS_DEFAULT |
    Py_TPFLAGS_BASETYPE,
    .tp_doc = "numc.Matrix objects",
    .tp_methods = Matrix61c_methods,
    .tp_members = Matrix61c_members,
    .tp_as_mapping = &Matrix61c_mapping,
    .tp_init = (initproc)Matrix61c_init,
    .tp_new = Matrix61c_new
};


struct PyModuleDef numcmodule = {
    PyModuleDef_HEAD_INIT,
    "numc",
    "Numc matrix operations",
    -1,
    Matrix61c_class_methods
};

/* Initialize the numc module */
PyMODINIT_FUNC PyInit_numc(void) {
    PyObject* m;

    if (PyType_Ready(&Matrix61cType) < 0)
        return NULL;

    m = PyModule_Create(&numcmodule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&Matrix61cType);
    PyModule_AddObject(m, "Matrix", (PyObject *)&Matrix61cType);
    printf("numc imported!\n");
    fflush(stdout);
    return m;
}