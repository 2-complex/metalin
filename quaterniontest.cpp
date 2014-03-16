
#include "quaternion.h"

/*	Creates a 4-by-4 matrix which rotates around the given axis by the given.*/
Mat3 axisRotate(const Vec3& axis, double theta)
{
	double x = axis.x, y = axis.y, z = axis.z;
	double n = axis.mag();
	x /= n;
	y /= n;
	z /= n;
	double xx = x * x, yy = y * y, zz = z * z;
	double c = cos(theta);
	double s = sin(theta);
	double omc = 1 - c; // one minus cosine
	
	return Mat3(xx + (1 - xx) * c,
		x * y * omc + z * s,
		x * z * omc - y * s,
		x * y * omc - z * s,
		yy + (1.0 - yy) * c,
		y * z * omc + x * s,
		x * z * omc + y * s,
		y * z * omc - x * s,
		zz + (1.0 - zz) * c);
};

int main(int argc, char** args)
{
	Quaternion z(1,2,3,4), w(2,-3,1,-2);
	Quaternion one(0,0,0,1), i(1,0,0,0), j(0,1,0,0), k(0,0,1,0);
	
	z.display();
	
	(z+w).display();
	(z-w).display();
	(z*w).display();
	
	(z*8).display();
	(8*z).display();
	(z/8).display();
	
	printf( "i*j=" ); (i*j).display();
	printf( "j*i=" ); (j*i).display();
	
	printf( "j*k=" ); (j*k).display();
	printf( "k*j=" ); (k*j).display();
	
	printf( "i*k=" ); (i*k).display();
	printf( "k*i=" ); (k*i).display();
	
	printf( "k*one=" ); (k*one).display();
	printf( "one*j=" ); (one*j).display();
	
	printf( "z*z.inverse() = " ); (z*z.inverse()).display();
	printf( "w*w.inverse() = " ); (w*w.inverse()).display();
	
	printf( "\n" );
	
	Vec3 V(1,2,3);
	Mat3 M = axisRotate(Vec3(1, 1, 1), 2.0 * M_PI / 3.0);
	
	printf("M=\n");
	M.display();
	
	printf("V = "); V.display();
	printf("M*V = "); (M*V).display();
	
	Quaternion qm(M);
	
	printf("qm = "); qm.display();
	
	printf("qm.matrix() = \n");
	(qm.matrix()).display();
	
	printf("qm.matrix()*qm.matrix().transpose() = \n");
	(qm.matrix()*qm.matrix().transpose()).display();
	
	printf("qm.rotate(V) = "); (qm.rotate(V)).display();

	
	return 0;
}

