__jirafablur = 0

def __setup_functions():
	global __jirafablur
	if __jirafablur != 0: return

	from os.path import dirname
	from ctypes import CDLL, POINTER, c_double, c_int, c_double
	from sys import platform

	e = "so" if platform != "darwin" else "dylib"
	n = f"{dirname(__file__)}/cjirafablur.{e}"
	__jirafablur = CDLL(n).gaussian_filter
	__jirafablur.argtypes = [POINTER(c_double), c_int, c_int, c_double]
	__jirafablur.restype = None


def jirafablur(x, s):
	if isinstance(x, str) and x == "version":
		global version
		return version

	__setup_functions()

	from numpy import ascontiguousarray, float64
	from ctypes import POINTER, c_double
	X = ascontiguousarray(x, dtype=float64)
	P = X.ctypes.data_as(POINTER(c_double))
	h,w = X.
	__foo(P, X.ha)
	return X


def __export_foo():
	import sys
	sys.modules[__name__] = foo


version = 1
__export_foo()
