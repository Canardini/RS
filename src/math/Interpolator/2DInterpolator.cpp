#include<Math\Interpolator\2DInterpolator.h>


void TwoDInterpolator::findNeighboursPoints(std::vector<int> & outNeighbours, int &Num, double inX,double inY){
	Num = 4;
	
	int lIndexX = locate(mX, inX);
	int lIndexY = locate(mY, inY);

	if ((lIndexX == -1) | (lIndexX == mX.size() - 1)){
		Num/=2;
	}
	if ((lIndexY == -1) | (lIndexY == mY.size() - 1)){
		Num /= 2;
	}
	outNeighbours.resize(2 * Num);


	if (Num == 4){
		outNeighbours[0] = lIndexX;
		outNeighbours[1] = lIndexY;
		outNeighbours[2] = lIndexX;
		outNeighbours[3] = lIndexY+1;
		outNeighbours[4] = lIndexX + 1;
		outNeighbours[5] = lIndexY;
		outNeighbours[6] = lIndexX + 1;
		outNeighbours[7] = lIndexY+1;
	}
	else if(Num == 2){
		if (lIndexX == -1){
			outNeighbours[0] = 0;
			outNeighbours[1] = lIndexY;
			outNeighbours[2] = 0;
			outNeighbours[3] = lIndexY+1;
		}
		else if (lIndexY == -1){
			outNeighbours[0] = lIndexX;
			outNeighbours[1] = 0;
			outNeighbours[2] = lIndexX+1;
			outNeighbours[3] =0;
		}
		else if (lIndexX == mX.size()-1){
			outNeighbours[0] = lIndexX;
			outNeighbours[1] = lIndexY;
			outNeighbours[2] = lIndexX;
			outNeighbours[3] = lIndexY+1;
		}

		else if (lIndexY == mY.size() - 1){
			outNeighbours[0] = lIndexX;
			outNeighbours[1] = lIndexY;
			outNeighbours[2] = lIndexX+1;
			outNeighbours[3] = lIndexY;
		}
	}
	else{
		if (lIndexX == -1){
			if (lIndexY == -1){
				outNeighbours[0] = 0;
				outNeighbours[1] = 0;
			}
			else{
				outNeighbours[0] = 0;
				outNeighbours[1] = lIndexY;
			}
		}
		else{
			if (lIndexY== -1){
				outNeighbours[0] = lIndexX;
				outNeighbours[1] = 0;
			}
			else{
				outNeighbours[0] = lIndexX;
				outNeighbours[1] = lIndexY;
			}
		}
	}

}

int TwoDInterpolator::locate(const std::vector<double> & xx, double inX0){
	int ju, jm, jl;
	int mm = 1;
	int n = xx.size();
	if (inX0 > xx[n - 1]){
		return n - 1;
	}
	if (inX0 < xx[0]){
		return -1;
	}

	bool ascnd = (xx[n - 1] >= xx[0]);
	jl = 0;
	ju = n - 1;
	while (ju - jl > 1) {
		jm = (ju + jl) >> 1;
		if (inX0 >= xx[jm] == ascnd)
			jl = jm;
		else
			ju = jm;
	}

	return jl;
}


void TwoDInterpolator::shiftXYZ(int i, int j, double shift){
	mZ[i][j] += shift;
}