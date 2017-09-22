#include <Math/Matrix/Matrix.h>


 void Matrix::fillMatrixWithVector(dMatrix & m, const dVector &a, int j, int n){
	int i;
	for (i = 0; i < n; i++){
		m[i][j] = a[i];
	}
}
 
double Matrix::dotProduct(dVector u, dVector v)
{
	int n = u.size();
	double sum = 0.0;
	for (int i = 0; i < n; i++){
		sum += u[i] * v[i];
	}

	return sum;
}
