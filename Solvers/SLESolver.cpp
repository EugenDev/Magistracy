#include "SLESolver.hpp"

void SolveSLE(const Matrix &A, Matrix &Z, const Matrix &b, float initialZ)
{
	int aCols = A.GetColsCount();
	int aRows = A.GetRowsCount();
	int zCols = Z.GetColsCount();
	int zRows = Z.GetRowsCount();
	int bCols = b.GetColsCount();
	int bRows = b.GetRowsCount();

	if (aCols != aRows || zCols != 1 || bCols != 1 || zRows != bRows || zRows != aCols)
		throw("SLE solver: This is not a SLE!\n");
	
	printf(" Solving SLE ... \n");
	clock_t time = clock();

	int m = aCols;
	float normb = b.Norm();

	Matrix rk(m, 1);
	Matrix pk(m, 1);
	Matrix Ark(m, 1);
	Matrix Apk(m, 1);
	Matrix zp(m, 1);
	Matrix zpp(m, 1);
	Matrix tmpV(m, 1);

	zpp.Fill(initialZ);
	pk.Fill(0.0f);
	rk.Fill(0.0f);
	Ark.Fill(0.0f);
	Apk.Fill(0.0f);
	zp.Fill(0.0f);
	tmpV.Fill(0.0f);

	rk = A * zpp - b;
	float normr = rk.Norm();

	Ark = A * rk;
	float d = Ark.Dot(rk);

	zp = zpp - rk * (normr * normr / d);

	
	int flag = 1;
	int iterations = 1;
	
	while (flag == 1)
	{
		rk = A * zp - b;
		normr = rk.Norm();

		pk = zp - zpp;

		Ark = A * rk;
		Apk = A * pk;

		float dot1 = Ark.Dot(pk);
		float dot2 = rk.Dot(pk);
		float dot3 = Ark.Dot(rk);
		float dot4 = Apk.Dot(pk);
		d = dot3 * dot4 - dot1 * dot1;

		float gamma = ((normr * normr) * dot4 - dot2 * dot1) / d;
		float beta = ((normr * normr) * dot1 - dot2 * dot3) / d;
			
		zpp = zp;
		zp -= rk * gamma - pk * beta; 

		tmpV = A * zp - b;
		double norm = tmpV.Norm();
			
		double error = norm / normb;
		
		printf("   Iteration:%d\terror:%f\n", iterations, error);
			
		if (error < 0.001)
			flag = 0;

		iterations++;
	}
	
	printf(" SLE solved with %d iterations for %f\n", iterations, (double)(clock() - time) / (CLOCKS_PER_SEC * 60));

	Z = zp;
	
	return;
}