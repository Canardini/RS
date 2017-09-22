#ifndef MATRIX_H
#define MATRIX_H

#include<vector>

#pragma once


    typedef std::vector<double> dVector;
    typedef std::vector<dVector> dMatrix;

	namespace Matrix{
	  double dotProduct(dVector u, dVector v);
	  void fillMatrixWithVector(dMatrix & m, const dVector & a, int i, int n);
	}
 

#endif

