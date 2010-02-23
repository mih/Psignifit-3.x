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


#core = swignifit.abCore()
#sigmoid = swignifit.PsiLogistic()
#psi = swignifit.PsiPsychometric(2,core,sigmoid)
#params = swignifit.vectord([0.5,0.5,0.01])
#x = numpy.arange(0,10,0.1)
#y = numpy.zeros(len(x))
#for i,val in enumerate(x):
#    y[i] = psi.evaluate(val,params)
#pylab.plot(x,y)
#pylab.show()

if __name__ == "__main__":
    unittest.main()

