# numc - a version of numpy
A C and a performance project to build a version of numpy. Numc is not as fast as numpy, but much faster than the naive implementation of matrix operations.

## Contents
* [Prerequisite](#prerequisite)
* [numc as a Python Library](#numc-as-a-python-library)
    * [installing numc](#installing-numc)
    * [Importing numc.Matrix](#importing-numcmatrix)
    * [Initializing numc.Matrix](#initializing-numcmatrix)
    * [Instance Attributes](#instance-attributes)
    * [Python/C API Reference](#pythonc-api-reference)
* [Libraries and Technologies used](#libraries-and-technologies-used)
* [Python-C Interface](#python-c-interface)
    * [Number Methods](#number-methods)
    * [Instance Methods](#instance-methods)
    * [Indexing](#indexing)
* [Matrix Operations Speedup](#matrix-operations-speedup)
    * [Add, Sub, Neg and Abs](#add-sub-neg-and-abs)
    * [Mul](#mul)
    * [Pow](#pow)
    

## Prerequisite
In order for numc to work, your local computer must have support for OpenMP library and Intel AVX intrinsics. I also used python version 3.6 throughout the project, so I would recommend you to use the same version just in case to avoid any unexpected bugs.     

## numc as a Python Library
### Installing `numc`
* You should be able to install `numc` by running:
```bash
$ make      # this will install the numc library. You might need to install make beforehand if you haven't already
$ python3.6     # this will open up python interpreter on your terminal and you are now ready to import numc just like any other python library!!
Python 3.6.2
[GCC 4.2.1 (Apple Inc. build 5666) (dot 3)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>>   
```
* The command above works for me on mac. If you are on a different OS you might need to experment between using `python3` and `python`. If `python3` doesn't works for you then you might need to change line 11 on `Makefile` -- change from `python3` to `python`. And similarly to initialize a python interpreter on terminal use `python` instead of `python3`
* Upon successful installation, `numc.Matrix` will be initialized and ready to import.
* You can uninstall `numc` module by running:
```bash
$ make unisntall
```
### Importing `numc.Matrix`
* Here are several ways of importing `numc.Matrix`
```bash
from numc import Matrix

import numc
numc.Matrix 

import numc as nc
nc.Matrix
```
### Initializing `numc.Matrix`
* Here are all different ways of creating `numc.Matrix` objects
```bash
>>> import numc as nc
numc imported!
>>> nc.Matrix(3, 3) # This creates a 3 * 3 matrix with entries all zeros
[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
>>> nc.Matrix(3, 3, 1) # This creates a 3 * 3 matrix with entries all ones
[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]
>>> nc.Matrix([[1, 2, 3], [4, 5, 6]]) # This creates a 2 * 3 matrix with first row 1, 2, 3, second row 4, 5, 6
[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
>>> nc.Matrix(1, 2, [4, 5]) # This creates a 1 * 2 matrix with entries 4, 5
[4.0, 5.0]
```

### Instance attributes
* The matrices and vectors have an attribute shape. For a 2D matrix, the shape is a tuple `(rows, cols)`. For a 1D matrix, it is a one-element tuple of (number of elements, ). Example is given below.
```bash
>>> import numc as nc
numc imported!
>>> mat = nc.Matrix(3, 3)
>>> mat.shape
(3, 3)
>>> mat = nc.Matrix(3, 1)
>>> mat.shape
(3,)
```

### Python/C API Reference
If you ever feel lost on how `numc.c`/Python-C interface is implemented, this [full reference manual](https://docs.python.org/3.6/c-api/index.html) is a great resource. 

## Libraries and Technologies used
* C
* Python 3.6
* [Intel AVX Intrinsics](https://software.intel.com/sites/landingpage/IntrinsicsGuide/)
* [OpenMP](https://www.openmp.org//wp-content/uploads/OpenMP3.0-SummarySpec.pdf)

## Python-C Interface

All these methods below are called through Python-C interface. In other words, these functions are called when you do matrix operations with `numc.Matrix`.

### Number Methods
* All these number methods returns a `numc.Matrix` object.
* After these number methods are implemented, PyNumberMethods struct `Matrix61c_as_number` is filled out which is used to define the object type `numc.Matrix`.
* Refer to the [official documentation](https://docs.python.org/3/c-api/typeobj.html#c.PyNumberMethods) of PyNumberMethods struct if this was confusing. 

Operator | Function | Description
-------- | -------- | ---------
a + b | `Matrix61c_add` | Element-wise sum of a and b.
a - b | `Matrix61c_sub` | Element-wise subtraction of a and b.
a * b | `Matrix61c_multiply` | Matrix multiplication of a and b. Remember that this is a matrix multiplication, not an element-wise multiplication.
- a | `Matrix61c_neg` | Element-wise negation of a.
abs(a) | `Matrix61c_abs` | Element-wise absolute value of a.
a** pow | `Matrix61c_pow` | Raise a to the powth power. a to the 0th power is the identity matrix (1 on the top left to bottom right diagonal and 0 everywhere else). This operator is defined in terms of matrix multiplication, not element-wise multiplication.

### Instance Methods
* After implementing these instance methods below, PyMethodDef[] structs `Matrix61C_methods` is filled out, which is used to define the object `numc.Matrix`.
* Refer to the [PyMethodDef documentation](https://docs.python.org/3/c-api/structures.html) if you need more clarification on how the instance method are linked on Python/C interface.
* For this project, we only have two instance methods: `set` and `get`

Python method | C Function | Description
-------- | --------- | ------
`set` | `Matrix61c_set_value` | Set self’s entry at the `i`th row and `j`th column to `val`.
`get` | `Matrix61c_get_value` | Returns the entry at the `i`th row and `j`th column.

### Indexing
#### Matrix61c_subscript
* Takes in a `numc` matrix and a key to index into the matrix. Below are examples of different cases.
```bash
>>> import numc as nc
  numc imported!
  >>> a = nc.Matrix(3, 3) 
  >>> a[0] # Key is a single number
  [0.0, 0.0, 0.0]
  >>> a[0:2] # key is a single slice
  [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  >>> a[0:2, 0:2] # key is a tuple of two slices
  [[0.0, 0.0], [0.0, 0.0]]
  >>> a[0:2, 0] # key is a tuple of (slice, int)
  [0.0, 0.0]
  >>> a[0, 0:2] # key is a tuple of (int, slice)
  [0.0, 0.0]
  >>> a[0, 0] # key is a tuple of (int, int)
  0.0
  >>> b = nc.Matrix(1, 3) # b is a 1D matrix
  >>> b[0]
  0.0
  >>> b[0:2] # Number of rows/cols does not matter now. You are slicing it as if it were a list
  [0.0, 0.0]
  >>> b[0:1, 0:1] # This is invalid!
  Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
  TypeError: 1D matrices only support single slice!
```
* returns a single number if the resulting slice is 1 by 1, otherwise it returns a new matrix that shares its data with its parent matrix. Here are more examples.
```bash
  >>> import numc as nc
  numc imported!
  >>> a = nc.Matrix(3, 3)
  >>> a[0][1]
  0.0
  >>> a[0:1, 0:1]
  0.0
```
* More examples. start:stop
```bash
>>> import numc as nc
numc imported!
>>> a = nc.Matrix(3, 1, [1, 2, 3])
>>> a[0]
1.0
>>> b = nc.Matrix(1, 3, [1, 2, 3])
>>> b[0]
1.0
>>> a[1:3]
[2.0, 3.0]
>>> b[1:3]
[2.0, 3.0]
```
* `numc` doesn't support the following
```bash
  >>> import numc as nc
  numc imported!
  >>> a = nc.Matrix(4, 4)
  >>> a[0:4:2] # Step size != 1
  Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
  ValueError: Slice info not valid!
  >>> a[0:0] # Slice has length < 1
  Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
  ValueError: Slice info not valid!
```

#### Matrix61c_set_subscript
* takes in a `numc` matrix, a `key` to index into the matrix and a value `v` to which to set the new slice. Below are examples of different cases:
```bash
>>> import numc as nc
numc imported!
>>> a = nc.Matrix(3, 3)
>>> a[0:1, 0:1] = 0.0 # Resulting slice is 1 by 1
>>> a[:, 0] = [1, 1, 1] # Resulting slice is 1D
>>> a
[[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
>>> a[0, :] = [2, 2, 2] # Resulting slice is 1D
>>> a
[[2.0, 2.0, 2.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
>>> a[0:2, 0:2] = [[1, 2], [3, 4]] # Resulting slice is 2D
>>> a
[[1.0, 2.0, 2.0], [3.0, 4.0, 0.0], [1.0, 0.0, 0.0]]
```
* slices share data with the original matrices which means that changing the values of slices will also change the value of original matrices. Here are some examples:
```bash
 >>> import numc as nc
numc imported!
>>> a = nc.Matrix(2, 2)
>>> a[0:1, 0:1] = 1.0
>>> a
[[1.0, 0.0], [0.0, 0.0]]
>>> a[1] = [2, 2]
>>> a
[[1.0, 0.0], [2.0, 2.0]]
>>> b = a[1]
>>> b[1] = 3
>>> a
[[1.0, 0.0], [2.0, 3.0]]
 ```
 * it is also possible to have nested slices, and changing nested slice's data should also change the values of original matrices. For example:
```bash
 >>> import numc as nc
numc imported!
>>> a = nc.Matrix(4, 4)
>>> b = a[0:3, 0:3]
>>> c = b[1:3, 1:3]
>>> c[0] = [2, 2] # Changing c should change both a and b
>>> c
[[2.0, 2.0], [0.0, 0.0]]
>>> b
[[0.0, 0.0, 0.0], [0.0, 2.0, 2.0], [0.0, 0.0, 0.0]]
>>> a
[[0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 2.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]]
```

## Matrix Operations Speedup
### Add, Sub, Neg and Abs
* avoid if/else statements and function calls when possible
    * if/else statements and function calls takes time setting up arguments and jumping to different part of the code 
* unroll loop when possible 
    * loop unrolling mostly/always helpful because less iterations implies less branch instructions which implies less stalling which then implies less time waiting

### Mul
* take advantage of cache
    * different loop orderings access items differnely. [Here](https://introcs.cs.princeton.edu/java/95linear/MatrixMultiplication.java.html) is list of different loop orderings and retrieve time comparision between them.
* a lot of iterations 
    * unrolling is always helpful
* same operations done simultaneously on many elements one by one,
    * why not do do them in vectors (four elements at a time)
  
### Pow
* Repeated square algorithmn
    * ~1700X speedup compared to naive implementation on a machine which has 4 cores with two hyperthreads each.
* Pseudocode:
 ```python
 def pow(base, exp):
    result = 1;
    while (exp > 0):
        if (exp & 1):
            result *= base
        exp /= 2
        base = base * base
    return result
 ```
