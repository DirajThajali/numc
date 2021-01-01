from utils import *
from unittest import TestCase

"""
For each operation, you should write tests to test  on matrices of different sizes.
Hint: use dp_mc_matrix to generate dumbpy and numc matrices with the same data and use
      cmp_dp_nc_matrix to compare the results
"""
class TestAdd(TestCase):
    def test_small_add(self):
        dp_mat1, nc_mat1 = dp_nc_matrix([[1,2],[3,4]])
        dp_mat2, nc_mat2 = dp_nc_matrix([[1,2],[3,4]])
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
    def test_small_add(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(200, 200, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(200, 200, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_add(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(3000, 3000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(3000, 3000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_mid_large_add(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(8000, 8000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(8000, 8000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_add(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(10000, 10000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(10000, 10000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

class TestSub(TestCase):
    def test_small_sub(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(200, 200, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(200, 200, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "sub")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_sub(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(4000, 4000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(4000, 4000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "sub")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_mid_large_sub(self):
        # TODO: YOUR CODE HERE
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(8000, 8000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(8000, 8000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "sub")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_sub(self):
        # TODO: YOUR CODE HERE
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(10000, 10000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(10000, 10000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "sub")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

class TestAbs(TestCase):
    def test_simple_abs(self):
        # TODO: YOUR CODE HERE
        dp_mat = nc.Matrix(3,2, -8)
        nc_mat = nc.Matrix(3,2, -8)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
        dp_mat = nc.Matrix(100,200, 8)
        nc_mat = nc.Matrix(200,100, 8)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_rand_med_abs(self):
        # TODO: YOUR CODE HERE
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(4000, 4000, seed=0)
        is_correct, speed_up = compute([dp_mat1], [nc_mat1], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(8000, 8000, seed=0)
        is_correct, speed_up = compute([dp_mat1], [nc_mat1], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_abs(self):
        # TODO: YOUR CODE HERE
        dp_mat = nc.Matrix(10000,10000, -8)
        nc_mat = nc.Matrix(10000,10000, -8)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)


class TestNeg(TestCase):
    def test_small_neg(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(200, 200, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
        dp_mat, nc_mat = rand_dp_nc_matrix(500, 500, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_neg(self):
        # TODO: YOUR CODE HERE
        dp_mat = nc.Matrix(4000,4000, 1)
        nc_mat = nc.Matrix(4000,4000, 1)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_mid_large_neg(self):
        # TODO: YOUR CODE HERE
        dp_mat = nc.Matrix(8000,8000, -1)
        nc_mat = nc.Matrix(8000,8000, -1)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_neg(self):
        # TODO: YOUR CODE HERE
        dp_mat = nc.Matrix(10000,10000, -1)
        nc_mat = nc.Matrix(10000,10000, -1)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

class TestMul(TestCase):
    def test_small_mul(self):
        # TODO: YOUR CODE HERE
        dp_mat1, nc_mat1 = dp_nc_matrix([
            [1,2,3,4],
            [5,6,7,8]])
        dp_mat2, nc_mat2 = dp_nc_matrix([
            [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],
            [17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32],
            [33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48],
            [49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64]])
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
        #print(dp_mat1 * dp_mat2)
        #print(nc_mat1 * nc_mat2)
        self.assertTrue(is_correct)
        print_speedup(speed_up)

        dp_mat1, nc_mat1 = rand_dp_nc_matrix(320, 420, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(420, 500, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_mul(self):
        # TODO: YOUR CODE HERE
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(1000, 1000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(1000, 1000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_med_large_mul(self):
        # TODO: YOUR CODE HERE
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(2000, 2000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(2000, 2000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_mul(self):
        # TODO: YOUR CODE HERE
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(5000, 5000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(5000, 5000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

class TestPow(TestCase):
    def test_known_pow(self):
        dp_mat, nc_mat = dp_nc_matrix([[1,2],[3,4]])
        is_correct, speed_up = compute([dp_mat, 4], [nc_mat, 4], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_pow_edge(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(10, 10, seed=10)
        is_correct, speed_up = compute([dp_mat, 0], [nc_mat, 0], "pow")
        print_speedup(speed_up)
        self.assertTrue(is_correct)
        is_correct, speed_up = compute([dp_mat, 1], [nc_mat, 1], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
        dp_mat, nc_mat = rand_dp_nc_matrix(1, 1, seed=10)
        is_correct, speed_up = compute([dp_mat, 1], [nc_mat, 1], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_small_pow(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(200, 200, seed=0)
        is_correct, speed_up = compute([dp_mat, 100], [nc_mat, 100], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
    
    def test_medium_pow(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(500, 500, seed=10)
        is_correct, speed_up = compute([dp_mat, 1000], [nc_mat, 1000], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_pow(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(500, 500, seed=10)
        is_correct, speed_up = compute([dp_mat, 2000], [nc_mat, 2000], "pow")
        self.assertTrue(is_correct)