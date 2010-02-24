import swignifit, numpy, pylab, unittest

class TestSigmoid(unittest.TestCase):
    """ test that all sigmoids have been wrapped and can be executed """

    def all_methods(self, sigmoid):
        sigmoid.f(0.0)
        sigmoid.df(0.0)
        sigmoid.ddf(0.0)
        sigmoid.inv(0.1)

    def test_cauchy(self):
        self.all_methods(swignifit.PsiCauchy())

    def test_exponential(self):
        self.all_methods(swignifit.PsiExponential())

    def test_gauss(self):
        self.all_methods(swignifit.PsiGauss())

    def test_gumbell(self):
        self.all_methods(swignifit.PsiGumbelL())

    def test_gumbelr(self):
        self.all_methods(swignifit.PsiGumbelR())

    def test_logistic(self):
        self.all_methods(swignifit.PsiLogistic())

    def test_exponential_exception(self):
        s = swignifit.PsiExponential()
        self.assertRaises(ValueError, s.inv, 0)
        self.assertRaises(ValueError, s.inv, 1)

class TestData(unittest.TestCase):

    def test_data(self):
        x = swignifit.vector_double([0.,2.,4.,6., 8., 10.])
        k = swignifit.vector_int([24, 32, 40,48, 50,48])
        n = swignifit.vector_int(6*[50])
        data = swignifit.PsiData(x, n, k, 2)

class TestCore(unittest.TestCase):
    x = swignifit.vector_double([0.,2.,4.,6., 8., 10.])
    k = swignifit.vector_int([24, 32, 40,48, 50,48])
    n = swignifit.vector_int(6*[50])
    data = swignifit.PsiData(x, n, k, 2)

    def all_methods(self, core):
        params = swignifit.vector_double([1.0,1.0])
        core.g(0.0, params)
        core.dg(0.0,params,0)
        core.dg(0.0,params,1)
        core.ddg(0.0,params,0,0)
        core.ddg(0.0,params,0,1)
        core.ddg(0.0,params,1,0)
        core.ddg(0.0,params,1,1)
        core.inv(0.0,params)
        core.dinv(0.0,params,0)
        core.dinv(0.0,params,1)
        core.transform(2,1.0,1.0)

    def test_ab_core(self):
        self.all_methods(swignifit.abCore())

    def test_linear_core(self):
        self.all_methods(swignifit.linearCore())

    def test_log_core(self):
        self.all_methods(swignifit.logCore(TestCore.data))

    def test_mw_core(self):
        # mwCore constructor is a bit different than the rest
        self.all_methods(swignifit.mwCore(1))

    def test_poly_core(self):
        self.all_methods(swignifit.polyCore(TestCore.data))

    def test_weibull_core(self):
        self.all_methods(swignifit.weibullCore(TestCore.data))

#core = swignifit.abCore()
#sigmoid = swignifit.PsiLogistic()
#psi = swignifit.PsiPsychometric(2,core,sigmoid)
#params = swignifit.vector_double([0.5,0.5,0.01])
#x = numpy.arange(0,10,0.1)
#y = numpy.zeros(len(x))
#for i,val in enumerate(x):
#    y[i] = psi.evaluate(val,params)
#pylab.plot(x,y)
#pylab.show()

if __name__ == "__main__":
    unittest.main()

