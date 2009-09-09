from distutils.core import setup, Extension

psippmodule = Extension ( "_psipy",
		sources = [
# 			"bootstrap.h",
# 			"core.h",
# 			"data.h",
# 			"errors.h",
# 			"mclist.h",
# 			"mcmc.h",
# 			"optimizer.h",
# 			"prior.h",
# 			"psipp.h",
# 			"psychometric.h",
# 			"rng.h",
# 			"sigmoid.h",
# 			"special.h",
			"bootstrap.cc",
			"core.cc",
			"data.cc",
			"mclist.cc",
			"mcmc.cc",
			"optimizer.cc",
			"psipy.cc",
			"psychometric.cc",
			"rng.cc",
			"sigmoid.cc",
			"special.cc"]
		)

setup ( name = "psipy",
		version = "0.1",
		description = "psipy psychometric functions in python",
		ext_modules = [psippmodule] )

