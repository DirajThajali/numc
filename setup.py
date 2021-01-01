from distutils.core import setup, Extension
import sysconfig

def main():
	CFLAGS = ['-g', '-Wall', '-std=c99', '-fopenmp', '-mavx', '-mfma', '-pthread', '-O3'] # change -O0 back to -O3 before submitting and after debugging
	LDFLAGS = ['-fopenmp']
	# Use the setup function we imported and set up the modules.
	# You may find this reference helpful: https://docs.python.org/3.6/extending/building.html
	module1 = Extension('numc',
						sources = ['numc.c', 'matrix.c'],
						extra_compile_args = CFLAGS,
						extra_link_args = LDFLAGS)

	setup (name = 'numc',
		   version = '1.0',
		   description = 'Numc - slower version of numpy',
		   ext_modules = [module1])

if __name__ == "__main__":
	main()
