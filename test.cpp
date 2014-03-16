

#include "lin.h"


int main(int argc, char** args)
{
	{
		printf( "default constructor\n" );
		Vec2 X, Y;
		X.display();
		Y.display();
		printf( "\n" );
	}
	
	{
		printf( "arithmetic\n" );
		Vec2 X(2,3), Y(4,8);
		X.display();
		Y.display();
		(X+Y).display();
		(X-Y).display();
		(X*Y).display();
		(X/Y).display();
		(3*X).display();
		(X*3).display();
		(X/3).display();
		(3/X).display();
		printf( "\n" );
	}
	
	{
		printf( "in-place arithmetic\n" );
		Vec2 X(2,3), Y(4,8);
		X /= Y;
		X.display();
		X *= 2;
		X.display();
		X /= 2;
		X.display();
		X += Y;
		X.display();
		X -= Y;
		X.display();
		printf( "\n" );
	}
	
	{
		printf( "matrix mul and inverses\n");
		Mat2 M1(1,2, 3,4);
		Mat2 M2(Vec2(1,2), Vec2(3,4));
		M1.display();
		M2.display();
		
		Mat2 A, B(1,2, 3,4);
		A.display();
		
		(A*B).display();
		
		A = B;
		A.inverse().display();
		Mat2 Ainv(A.inverse());
		(Ainv * A).print(); printf( " should be identity.\n" );
		(A * Ainv).display();
		
		M1.transpose().display();
		M1.transpose().inverse().print(); printf( " should be the inverse transpose.\n" );
		printf( "\n" );
	}
	
	{
		printf( "three dimensional vector arithmetic\n" );
		Vec3 X, Y;
		X.display();
		Y.display();
		X = Vec3(1, 2, 3);
		Y = Vec3(-2, 3, 1);
		(X + Y).display();
		(X * Y).display();
		(X / Y).display();
		
		// normalize(Vec3(1,1,1)).print(); printf( " should be unit length.\n" );
		// printf( "%f", normalize(Vec3(1,1,1)).mag()); printf( " should be 1.\n" );
		
		Vec3(1, 0, 0).cross(Vec3(0, 1, 0)).print(); printf( " should be (0, 0, 1)\n" );
		printf( "\n" );
	}

	{
		printf( "four dimensional vector arithmetic\n" );
		Vec4 X, Y;
		X.display();
		Y.display();
		X = Vec4(1, 2, 3, 0.1);
		Y = Vec4(-2, 3, 1, 3);
		(X + Y).display();
		(X * Y).display();
		(X / Y).display();
		
		// normalize(Vec4(1,1,1,2)).print(); printf( " should be unit length.\n" );
		// printf( "%f", normalize(Vec4(1, 1, 1, 4)).mag()); printf( " should be 1.\n" );
		printf( "\n" );
	}
	
	{
		printf( "arbitrary matrix, inverse, transpose, determinant and inverse-determinant are reciprocals\n" );
		Mat4 M(-1, 2, -1, 3,    1, 1, 1 ,-1,   -2, -1, 1, 3,   1, 2, 3, 1);
		printf( "M   = " ); M.display();
		printf( "Mt  = " ); M.transpose().display();
		printf( "Mi  = " ); M.inverse().display();
		printf( "Mit = " ); (M.inverse() * M).display();
		printf( "det(M ) = %f\n" , M.det() );
		printf( "det(Mi) = %f\n" , M.inverse().det() );
		printf( "\n" );
	}
}
