Metalin
=========================

Introducing Metalin, a python script... for generating c++ code... for doing linear algebra... for graphics.

In the repo, you will find `lingen.py` which is a python script which generates `lin.h`, a header file.  You can see the code it outputs by typing:

	python lingen.py

You will also see `test.cpp`, which includes `lin.h` and uses it.  In the `makefile`, the test depends on `lin.h` which in turn, depends on the python script.  If you build, by typing:

	make test

Make will run the python script to generate `lin.h`, and then build the test using c++:
	
	$ make test
	touch lin.h; rm lin.h; python lingen.py >> lin.h
	c++ test.cpp -o test

You can then run the test...

	./test

...and what you will see is the result of a bunch of matrix/vector arithmetic operations.  See the source file `test.cpp`.

You can then include `lin.h` into your own c++ project and use the classes to do whatever you need.

	#include "lin.h"
	
	int main()
	{
		Vec3 X(2, 3, 4);
		Vec3 Y(-2, 0, 1);

		(X*Y).display();

		...etc...
	}

Enjoy.
