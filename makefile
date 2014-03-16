
current: test

lin.h: lingen.py
	touch lin.h; rm lin.h; python lingen.py >> lin.h

test: test.cpp lin.h
	c++ test.cpp -o test

quaterniontest: quaternion.h
	c++ quaterniontest.cpp -o quaterniontest

