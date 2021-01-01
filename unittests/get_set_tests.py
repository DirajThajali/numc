from utils import *
from unittest import TestCase



"""
For each operation, you should write tests to test  on matrices of different sizes.
Hint: use dp_mc_matrix to generate dumbpy and numc matrices with the same data and use
      cmp_dp_nc_matrix to compare the results
"""
class TestGet(TestCase):
    # 1d or 2d _ key _ result
    # int
    def test_1d_int_num(self): # passed
        dp_mat, nc_mat = rand_dp_nc_matrix(1, 3, seed=0)
        rand_col = np.random.randint(dp_mat.shape[0])
        try:
            self.assertEqual(round(nc_mat[rand_col], decimal_places), round(dp_mat[rand_col], decimal_places))
        except IndexError as e:
            # print(e)
            pass

        dp_mat1, nc_mat1 = rand_dp_nc_matrix(3, 1, seed=0)
        rand_col1 = np.random.randint(dp_mat.shape[0])
        try:
            self.assertEqual(round(nc_mat1[rand_col1], decimal_places), round(dp_mat1[rand_col1], decimal_places))
        except IndexError as e:
            # print(e)
            pass

    # slice
    def test_1d_slice_list(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(1, 3, seed=0)
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[0:2], nc_mat[0:2]))


    def test_1d_slice_num(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(1, 3, seed=0)
        self.assertEqual(dp_mat[0:1], nc_mat[0:1])

    # tuple
    def test_1d_tuple_err(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(1, 3, seed=0)
        try:
            self.assertEqual(dp_mat[0:1, 0:1], nc_mat[0:1, 0:1])
        except TypeError as e:
            pass

    # int
    def test_2d_int_list(self):
        dp_mat, nc_mat = dp_nc_matrix([[1,2,3],[4,5,6],[7,8,9]])
        rand_row = np.random.randint(dp_mat.shape[0])
        self.assertEqual(dp_mat[rand_row], nc_mat[rand_row])

        out_of_range = dp_mat.shape[0] + 1
        try:
            self.assertEqual(dp_mat[out_of_range], nc_mat[out_of_range])
        except IndexError as e:
            pass
        
    # slice
    def test_2d_slice_list_of_list(self):
        dp_mat, nc_mat = dp_nc_matrix([[1,2,3],[4,5,6],[7,8,9]])
        rand_row_start = np.random.randint(dp_mat.shape[0])
        rand_row_end = np.random.randint(dp_mat.shape[0])
        try:
            self.assertEqual(dp_mat[rand_row_start:rand_row_end], nc_mat[rand_row_start:rand_row_end])
        except ValueError as e:
            pass

        dp_mat1, nc_mat1 = rand_dp_nc_matrix(4, 4, seed=0)
        try:
            self.assertEqual(dp_mat1[0:4:2], nc_mat1[0:4:2])
        except ValueError as e:
            pass

        try:
            self.assertEqual(dp_mat1[0:0], nc_mat1[0:0])
        except ValueError as e:
            pass

    # tuple
    def test_2d_int_int_num(self):
        dp_mat, nc_mat = dp_nc_matrix([[1,2,3],[4,5,6],[7,8,9]])
        # print(dp_mat.get(0,0))
        # print(nc_mat.get(0,0))

        rand_row = np.random.randint(dp_mat.shape[0])
        rand_col = np.random.randint(dp_mat.shape[1])
        self.assertEqual(dp_mat[rand_row][rand_col], nc_mat[rand_row][rand_col])

    def test_2d_int_slice_num(self):
        dp_mat, nc_mat = dp_nc_matrix([[1,2,3],[4,5,6],[7,8,9]])
        rand_row = np.random.randint(dp_mat.shape[0])
        rand_col_start = np.random.randint(dp_mat.shape[1])
        rand_col_end = np.random.randint(dp_mat.shape[1])
        try:
            self.assertEqual(dp_mat[rand_row][rand_col_start:rand_col_end], nc_mat[rand_row][rand_col_start:rand_col_end])
        except ValueError as e:
            pass

    def test_2d_slice_int_num(self):
        dp_mat, nc_mat = dp_nc_matrix([[1,2,3],[4,5,6],[7,8,9]])
        rand_row_start = np.random.randint(dp_mat.shape[0])
        rand_row_end = np.random.randint(dp_mat.shape[0])
        rand_col = np.random.randint(dp_mat.shape[1])
        try:
            self.assertEqual(dp_mat[rand_row_start:rand_row_end][rand_col], nc_mat[rand_row_start:rand_row_end][rand_col])
        except ValueError as e:
            pass

    def test_2d_slice_slice_num(self):
        dp_mat, nc_mat = dp_nc_matrix([[1,2,3],[4,5,6],[7,8,9]])
        rand_row_start = np.random.randint(dp_mat.shape[0])
        rand_row_end = np.random.randint(dp_mat.shape[0])
        rand_col_start = np.random.randint(dp_mat.shape[1])
        rand_col_end = np.random.randint(dp_mat.shape[1])
        try:
            self.assertEqual(dp_mat[rand_row_start:rand_row_end][rand_col_start:rand_col_end], nc_mat[rand_row_start:rand_row_end][rand_col_start:rand_col_end])
        except ValueError as e:
            pass

    def test_get_rand(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(2, 2, seed=0)
        rand_row = np.random.randint(dp_mat.shape[0])
        rand_col = np.random.randint(dp_mat.shape[1])
        self.assertEqual(round(dp_mat[rand_row][rand_col], decimal_places),
            round(nc_mat[rand_row][rand_col], decimal_places))

    def test_integration(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(2, 2, seed=0)
        rand_row = np.random.randint(dp_mat.shape[0])
        rand_col = np.random.randint(dp_mat.shape[1])
        self.assertEqual(round(dp_mat[rand_row][rand_col], decimal_places),
            round(nc_mat[rand_row][rand_col], decimal_places))

        dp_mat, nc_mat = rand_dp_nc_matrix(1000, 1, seed=0)
        rand_row = np.random.randint(dp_mat.shape[0])
        self.assertEqual(round(dp_mat[rand_row], decimal_places),
            round(nc_mat[rand_row], decimal_places))

        dp_mat, nc_mat = rand_dp_nc_matrix(14, 31, seed=0)
        rand_row = np.random.randint(dp_mat.shape[0])
        rand_col = np.random.randint(dp_mat.shape[1])
        self.assertEqual(round(dp_mat[rand_row][rand_col], decimal_places),
            round(nc_mat[rand_row][rand_col], decimal_places))
        self.assertEqual(round(dp_mat[0][0], decimal_places),
            round(nc_mat[0][0], decimal_places))
        self.assertEqual(round(dp_mat[13][30], decimal_places),
            round(nc_mat[13][30], decimal_places))

        a = nc.Matrix(3, 3, 1.5)
        
    def test_error(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(14, 31, seed=0)
        try:
            nc_mat.get(14, 2)
        except IndexError:
            pass

        try:
            nc_mat.get(12, 31)
        except IndexError:
            pass
        try:
            nc_mat.get(14, 2, 2)
        except TypeError:
            pass
        try:
            nc_mat.get(3.2, 2)
        except TypeError:
            pass
        
class TestSet(TestCase):
    def test_set(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(10, 10, seed=0)
        rand_row = np.random.randint(10)
        rand_col = np.random.randint(10)
        num = np.random.randint(100)
        dp_mat.set(rand_row, rand_col, num)
        nc_mat.set(rand_row, rand_col, num)
        self.assertEqual(round(dp_mat[rand_row][rand_col], decimal_places),
            round(nc_mat[rand_row][rand_col], decimal_places))

class TestShape(TestCase):
    def test_shape(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(2, 2, seed=0)
        self.assertTrue(dp_mat.shape == nc_mat.shape)
