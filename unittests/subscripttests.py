from utils import *
from unittest import TestCase

"""
For each operation, you should write tests to test  on matrices of different sizes.
Hint: use dp_mc_matrix to generate dumbpy and numc matrices with the same data and use
      cmp_dp_nc_matrix to compare the results
"""
class TestSubscript(TestCase):
      def test_row_1D(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 1, seed=0)
        nc_mat[-1:]

      def test_subscript_1D(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(1, 3, seed=0)
        #self.assertTrue(dp_mat[0] == nc_mat[0])
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[0:2], nc_mat[0:2]))
        
        self.assertTrue(isinstance(dp_mat[1:2], float))
        self.assertTrue(isinstance(nc_mat[1:2], float))
        self.assertTrue(dp_mat[1:2] == nc_mat[1:2])
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 1, seed=0)
        self.assertTrue(dp_mat[0] == nc_mat[0])
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[0:2],nc_mat[0:2]))
        self.assertTrue(isinstance(dp_mat[1:2], float))
        self.assertTrue(isinstance(nc_mat[1:2], float))
        self.assertTrue(dp_mat[1:2] == nc_mat[1:2])

      def test_subscript_1D_slice_all(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 1, seed=0)
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[:],nc_mat[:]))
        dp_mat, nc_mat = rand_dp_nc_matrix(1, 3, seed=0)
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[:],nc_mat[:]))

      def test_subscript_1D_error(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(1, 3, seed=0)
        try:
            nc_mat[3]
            self.assertTrue(False)
        except IndexError:
            pass
        try:
            nc_mat[3.3]
            self.assertTrue(False)
        except TypeError:
            pass
        try:
            nc_mat[1:1]
            self.assertTrue(False)
        except ValueError:
            pass
        try:
            nc_mat[0:2:2]
            self.assertTrue(False)
        except ValueError:
            pass
        
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 1, seed=0)
        try:
            nc_mat[3]
            self.assertTrue(False)
        except IndexError:
            pass
        try:
            nc_mat[3.3]
            self.assertTrue(False)
        except TypeError:
            pass
        try:
            nc_mat[1:1]
            self.assertTrue(False)
        except ValueError:
            pass
        try:
            nc_mat[0:2:2]
            self.assertTrue(False)
        except ValueError:
            pass
        
      def test_subscript_1D_neg(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(4, 1, seed=0)
        
        self.assertTrue(isinstance(dp_mat[-1:], float))
        self.assertTrue(isinstance(nc_mat[-1:], float))
        self.assertTrue(dp_mat[-1:] == nc_mat[-1:])
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[:-1],nc_mat[:-1]))
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[1:-1],nc_mat[1:-1]))

        dp_mat, nc_mat = rand_dp_nc_matrix(1, 4, seed=0)
        self.assertTrue(isinstance(dp_mat[-1:], float))
        self.assertTrue(isinstance(nc_mat[-1:], float))
        self.assertTrue(dp_mat[-1:] == nc_mat[-1:])
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[:-1],nc_mat[:-1]))
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[1:-1],nc_mat[1:-1]))
    
      def test_subscript_2D_error(self):
        a = nc.Matrix(3, 2, 5)
        try:
            a[1.2]
            self.assertTrue(False)
        except TypeError:
            pass
        try:
            a[0:2:2]
            self.assertTrue(False)
        except ValueError:
            pass
        try:
            a[0,8]
            self.assertTrue(False)
        except IndexError:
            pass
        try:
            a[6]
            self.assertTrue(False)
        except IndexError:
            pass
        a = nc.Matrix(2, 3, 5)
        try:
            a[1.2]
            self.assertTrue(False)
        except TypeError:
            pass
        try:
            a[0:2:2]
            self.assertTrue(False)
        except ValueError:
            pass
        try:
            a[6]
            self.assertTrue(False)
        except IndexError:
            pass
        try:
            a[0,8]
            self.assertTrue(False)
        except IndexError:
            pass

      def test_subscript_2D(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 3, seed=0)
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[0],nc_mat[0]))
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[0:2],nc_mat[0:2]))
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[0:2, 0:2],nc_mat[0:2, 0:2]))
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[0:2, 0],nc_mat[0:2, 0]))
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[0, 0:2],nc_mat[0, 0:2]))
        self.assertTrue(isinstance(dp_mat[0,0], float))
        self.assertTrue(isinstance(nc_mat[0,0], float))
        self.assertTrue(dp_mat[0,0] == nc_mat[0,0])
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 5, seed=0)
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[0],nc_mat[0]))
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[0:2],nc_mat[0:2]))
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[0:2, 0:2],nc_mat[0:2, 0:2]))
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[0:2, 0],nc_mat[0:2, 0]))
        self.assertTrue(cmp_dp_nc_matrix(dp_mat[0, 0:2],nc_mat[0, 0:2]))
        self.assertTrue(isinstance(dp_mat[0,0], float))
        self.assertTrue(isinstance(nc_mat[0,0], float))
        self.assertTrue(dp_mat[0,0] == nc_mat[0,0])

      def test_subscript_2D_error(self):
        a = nc.Matrix(3, 3, 5)
        try:
            a[1.2]
            self.assertTrue(False)
        except TypeError:
            pass

        try:
            a[0:2:2]
            self.assertTrue(False)
        except ValueError:
            pass
        try:
            a[0:2:1, 0:2:2]
            self.assertTrue(False)
        except ValueError:
            pass

        try:
            a[0:0]
            self.assertTrue(False)
        except ValueError:
            pass

        try:
            a[0:2, 1:1]
            self.assertTrue(False)
        except ValueError:
            pass
        try:
            a[6]
            self.assertTrue(False)
        except IndexError:
            pass
        a = nc.Matrix(3, 4, 5)
        try:
            a[1.2]
            self.assertTrue(False)
        except TypeError:
            pass

        try:
            a[0:2:2]
            self.assertTrue(False)
        except ValueError:
            pass
        try:
            a[0:2:1, 0:2:2]
            self.assertTrue(False)
        except ValueError:
            pass

        try:
            a[0:0]
            self.assertTrue(False)
        except ValueError:
            pass

        try:
            a[0:2, 1:1]
            self.assertTrue(False)
        except ValueError:
            pass
        try:
            a[6]
            self.assertTrue(False)
        except IndexError:
            pass


class TestSetSubscript(TestCase):
      def test_set_subscript_1D(self):
            dp_mat, nc_mat = rand_dp_nc_matrix(1, 5, seed=0)
            dp_mat[1:3] = [1, 2]
            nc_mat[1:3] = [1, 2]
            self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
            dp_mat, nc_mat = rand_dp_nc_matrix(1, 5, seed=0)
            dp_mat[1:3] = [1.5, 2.5]
            nc_mat[1:3] = [1.5, 2.5]
            self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
            dp_mat[1] = 2
            nc_mat[1] = 2
            self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
            dp_mat, nc_mat = rand_dp_nc_matrix(5, 1, seed=0)
            dp_mat[1:3] = [1, 2]
            nc_mat[1:3] = [1, 2]
            self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
            dp_mat[1] = 2
            nc_mat[1] = 2
            self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))

      def test_cmp(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(5, 1, seed=0)
        dp_mat[1:3] = [1, 2]
        nc_mat[1:3] = [1, 2]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))

      def test_set_subscript_1D_neg(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(1, 5, seed=0)
        dp_mat[2:-1] = [1, 2]
        nc_mat[2:-1] = [1, 2]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
        dp_mat, nc_mat = rand_dp_nc_matrix(5, 1, seed=0)
        dp_mat[2:-1] = [1, 2]
        nc_mat[2:-1] = [1, 2]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
        
    
      def test_set_subscript_1D_error(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(1, 5, seed=0)
        try:
            nc_mat[3.3] = 2
            self.assertTrue(False)
        except TypeError:
            pass
        
        try:
            nc_mat[1:3:2] = 2
            self.assertTrue(False)
        except ValueError:
            pass

        try:
            nc_mat[1:1] = 2
            self.assertTrue(False)
        except ValueError:
            pass

        try:
            nc_mat[1:7] = [2, 3, 4, 5, 6, 7]
            self.assertTrue(False)
        except ValueError:
            pass

        try:
            nc_mat[8] = 1
            self.assertTrue(False)
        except IndexError:
            pass
        
        try:
            nc_mat[2] = [1,2]
        except TypeError:
            pass
    
        try:
            nc_mat[1:3] = 3.3
        except TypeError:
            pass

        try:
            nc_mat[1:3] = [1, 2, 3]
        except ValueError:
            pass
        
        try:
            nc_mat[1:3] = ['a', 2]
        except ValueError:
            pass
        
        dp_mat, nc_mat = rand_dp_nc_matrix(5, 1, seed=0)
        try:
            nc_mat[3.3] = 2
            self.assertTrue(False)
        except TypeError:
            pass
        
        try:
            nc_mat[1:3:2] = 2
            self.assertTrue(False)
        except ValueError:
            pass

        try:
            nc_mat[1:1] = 2
            self.assertTrue(False)
        except ValueError:
            pass

        try:
            nc_mat[1:7] = [2, 3, 4, 5, 6, 7]
            self.assertTrue(False)
        except ValueError:
            pass

        try:
            nc_mat[8] = 1
            self.assertTrue(False)
        except IndexError:
            pass
        
        try:
            nc_mat[2] = [1,2]
        except TypeError:
            pass
    
        try:
            nc_mat[1:3] = 3.3
        except TypeError:
            pass

        try:
            nc_mat[1:3] = [1, 2, 3]
        except ValueError:
            pass
        
        try:
            nc_mat[1:3] = ['a', 2]
        except ValueError:
            pass

      def test_set_subscript_2D_1_by_1(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 3, seed=0)
        dp_mat[0:1, 0:1] = 0.0
        nc_mat[0:1, 0:1] = 0.0
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 4, seed=0)
        dp_mat[0:1, 0:1] = 0.0
        nc_mat[0:1, 0:1] = 0.0
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))

      def test_set_subscript_2D_all_rows(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 3, seed=0)
        dp_mat[:, 0] = [1, 1, 1]
        nc_mat[:, 0] = [1, 1, 1]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
        dp_mat[:, 2] = [1, 1, 1]
        nc_mat[:, 2] = [1, 1, 1]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))

      def test_set_subscript_2D_all_cols(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 3, seed=0)
        dp_mat[0, :] = [1, 1, 1]
        nc_mat[0, :] = [1, 1, 1]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
        dp_mat[2, :] = [1, 1, 1]
        nc_mat[2, :] = [1, 1, 1]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))

      def test_set_subscript_2D_rows_and_cols(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 3, seed=0)
        dp_mat[0:2, 0:2] = [[1, 2], [3, 4]]
        nc_mat[0:2, 0:2] = [[1, 2], [3, 4]]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 5, seed=0)
        dp_mat[0:2, 0:3] = [[1, 2, 3], [3, 4, 5]]
        nc_mat[0:2, 0:3] = [[1, 2, 3], [3, 4, 5]]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
    
      def test_set_subscript_2D_one_slice_all_cols(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 3, seed=0)
        dp_mat[1] = [2, 2, 2]
        nc_mat[1] = [2, 2, 2]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 5, seed=0)
        dp_mat[2] = [1, 2, 3, 4, 5]
        nc_mat[2] = [1, 2, 3, 4, 5]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
        dp_mat[0:2] = [[1, 2, 3, 4, 5], [1, 2, 3, 4, 5]]
        nc_mat[0:2] = [[1, 2, 3, 4, 5], [1, 2, 3, 4, 5]]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))

      def test_is_2D(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 5, seed=0)
        dp_mat[0:2] = [[1, 2, 3, 4, 5], [1, 2, 3, 4, 5]]
        nc_mat[0:2] = [[1, 2, 3, 4, 5], [1, 2, 3, 4, 5]]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))

      def test_set_subscript_2D_nested(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(2, 2, seed=0)
        dp_mat[0:1, 0:1] = 0.0
        nc_mat[0:1, 0:1] = 0.0
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
        dp_mat[1] = [2, 2]
        nc_mat[1] = [2, 2]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
        a = dp_mat[1]
        a[1] = 3
        b = nc_mat[1]
        b[1] = 3
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
        dp_mat, nc_mat = rand_dp_nc_matrix(4, 4, seed=0)
        a = dp_mat[0:3, 0:3]
        b = a[1:3, 1:3]
        b[0] = [2, 2]
        c = nc_mat[0:3, 0:3]
        d = c[1:3, 1:3]
        d[0] = [2, 2]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))

      def test_set_subscript_2D_neg(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 3, seed=0)

        dp_mat[0:-1] = [[1, 2, 3], [1, 2, 3]]
        nc_mat[0:-1] = [[1, 2, 3], [1, 2, 3]]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))

        dp_mat[1:-1, 1:-1] = 0.0
        nc_mat[1:-1, 1:-1] = 0.0
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))

        dp_mat[:,0:-1] = [[1, 1],[1, 1],[1, 1]]
        nc_mat[:,0:-1] = [[1, 1],[1, 1],[1, 1]]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
        dp_mat[2, :] = [1, 1, 1]
        nc_mat[2, :] = [1, 1, 1]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))

        dp_mat, nc_mat = rand_dp_nc_matrix(3, 3, seed=0)
        dp_mat[0:-1, 0:-1] = [[1, 2], [3, 4]]
        nc_mat[0:-1, 0:-1] = [[1, 2], [3, 4]]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))

        dp_mat, nc_mat = rand_dp_nc_matrix(3, 5, seed=0)
        dp_mat[0:-1, 0:-2] = [[1, 2, 3], [3, 4, 5]]
        nc_mat[0:-1, 0:-2] = [[1, 2, 3], [3, 4, 5]]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))

      def test_2d_neg_index(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 5, seed=0)
        dp_mat[0:-1, 0:-2] = [[1, 2, 3], [3, 4, 5]]
        nc_mat[0:-1, 0:-2] = [[1, 2, 3], [3, 4, 5]]
        self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))

      def test_set_subscript_error(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 3, seed=0)
        try:
            nc_mat["a"] = 3
        except TypeError:
            pass
        
        try:
            nc_mat[1:3:2, 1:3] = [2, 2]
        except ValueError:
            pass
        
        try:
            nc_mat[1:3, 1:3:2] = [2, 2]
        except ValueError:
            pass
        
        try:
            nc_mat[5] = [2,2,2]
        except IndexError:
            pass
        
        try:
            nc_mat[7, 1:3] = [2, 2]
        except IndexError:
            pass
        
        try:
            nc_mat[1:3, 7] = [2, 2]
        except IndexError:
            pass

        try:
            nc_mat[1:7, 1:2] = [2, 2, 3, 4, 5, 6]
        except ValueError:
            pass
        try:
            nc_mat[1:2, 1:7] = [2, 2, 3, 4, 5, 6]
        except ValueError:
            pass
        try:
            nc_mat[1:2, 1:2] = [1, 2]
        except TypeError:
            pass
        try:
            nc_mat[1:3, 1:7] = 5
        except TypeError:
            pass

        try:
            nc_mat[0:2, 0:2] = [[0, 2], [1, 2], [1, 2]]
        except ValueError:
            pass
        try:
            nc_mat[0:2, 0:2] = [[0, 2, 3], [1, 2, 3]]
        except ValueError:
            pass
        try:
            nc_mat[0:2, 0:2] = [[0, 2, 3], [1, 2]]
        except ValueError:
            pass
        try:
            nc_mat[0:2, 0:2] = [["a", 2], [1, 2]]
        except ValueError:
            pass
    