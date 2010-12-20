cython -o google_perftools.c google_perftools.pyx
gcc -o google_perftools.so -shared google_perftools.c -fPIC -lpython2.6 -lprofiler -I/usr/include/python2.6
