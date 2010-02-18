import swignifit, numpy, pylab
core = swignifit.abCore()
sigmoid = swignifit.PsiLogistic()
psi = swignifit.PsiPsychometric(2,core,sigmoid)
params = swignifit.vectord([0.5,0.5,0.01])
x = numpy.arange(0,10,0.1)
y = numpy.zeros(len(x))
for i,val in enumerate(x):
    y[i] = psi.evaluate(val,params)
pylab.plot(x,y)
pylab.show()




