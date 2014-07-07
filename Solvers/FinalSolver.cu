#include "FinalSolver.hpp"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "CudaHelpers.hpp"
#include "CudaOperatorReplacer.hpp"
#include "SLESolver.hpp"
#define GRAVY_CONST 6.673

GridData LightLinearisedMinimalError(Task task, float alpha)
{
	GridParameters gp = task.grid.GetGridParameters();
	int M = gp.NX * gp.NY;

	Matrix A(M, 1), F(M, 1), Z(M, 1);
	Z.Fill(task.initialZ);
	task.grid.FillMatrix(F);

	float error;
	int iteration = 1;

	CudaDirectSolver dslvr(gp, false);
	CudaOperatorReplacer oper(task);
	float *devF;
	float *devA;
	float *devZ;
	float *devTmp1;
					
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devF));
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devA));
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devZ));
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devTmp1));
		
	printf("Solving with Linearised Minimal Error Method\n");

	clock_t t0 = clock();

	while (true)
	{
		dslvr.SolveDirectTask(Z, A, task);

		error = (F - A).Norm() / F.Norm();

		printf("Iteration #%d.\tError = %f.\n", iteration, error);

		if (error < task.precision)
			break;

		CheckCublas(cublasSetVector(M, sizeof(float), A.elements, 1, devA, 1));
		CheckCublas(cublasSetVector(M, sizeof(float), F.elements, 1, devF, 1));
		CheckCublas(cublasSetVector(M, sizeof(float), Z.elements, 1, devZ, 1));
		
		//A = A - F
		cublasSaxpy(M, -1.0f, devF, 1, devA, 1);
		CheckCublas(cublasGetError());

		// ||(A - F)||
		float n1 = cublasSnrm2(M, devA, 1);

		//S = AnkT(A - F)
		oper.CalcATX(devZ, devA, devTmp1);

		// ||S||
		float n2 = cublasSnrm2(M, devTmp1, 1);

		float gamma = (alpha * n1 * n1) / (n2 * n2);

		cublasSaxpy(M, gamma, devTmp1, 1, devZ, 1);

		cublasGetVector(M, sizeof(float), devZ, 1, Z.elements, 1);

		iteration++;
	}

	printf("Solving finished for %f sec\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
	Z *= -1.0f; 
	GridData out(gp);
	memcpy(out.data, Z.elements, M * sizeof(float));

	CheckCublas(cublasFree(devTmp1));
	CheckCublas(cublasFree(devZ));
	CheckCublas(cublasFree(devA));
	CheckCublas(cublasFree(devF));

	return out;
}

GridData LightLinearisedSpeedDescent(Task task, float alpha)
{
	GridParameters gp = task.grid.GetGridParameters();
	int M = gp.NX * gp.NY;

	Matrix A(M, 1), F(M, 1), Z(M, 1);
	Z.Fill(task.initialZ);
	task.grid.FillMatrix(F);

	float error;
	int iteration = 1;

	CudaDirectSolver dslvr(gp, false);
	CudaOperatorReplacer oper(task);

	float *devF;
	float *devA;
	float *devZ;
	float *devS;
	float *devC;
					
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devF));
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devA));
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devZ));
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devS));
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devC));
		
	CheckCublas(cublasSetVector(M, sizeof(float), F.elements, 1, devF, 1));
	CheckCublas(cublasSetVector(M, sizeof(float), Z.elements, 1, devZ, 1));

	printf("Solving with Linearised Speed Descent Method\n");

	clock_t t0 = clock();

	while (true)
	{
		dslvr.SolveDirectTask(Z, A, task);

		CheckCublas(cublasSetVector(M, sizeof(float), A.elements, 1, devA, 1));
				
		//A = A - F
		cublasSaxpy(M, -1.0f, devF, 1, devA, 1);
		CheckCublas(cublasGetError());

		//---------------------Ошибка---------------------
		float err1 = cublasSnrm2(M, devA, 1);
		CheckCublas(cublasGetError());
		float err2 = cublasSnrm2(M, devF, 1);
		CheckCublas(cublasGetError());
		error = err1 / err2;

		printf("Iteration #%d.\tError = %f.\n", iteration, error);

		if (error < task.precision)
			break;
		//---------------------Ошибка---------------------
		
		//S = AnkT(A - F)
		oper.CalcATX(devZ, devA, devS);

		// C = Ank * S
		oper.CalcAX(devZ, devS, devC);
		
		// ||S||
		float n1 = cublasSnrm2(M, devS, 1);

		// ||C||
		float n2 = cublasSnrm2(M, devC, 1);

		float gamma = (alpha * n1 * n1) / (n2 * n2);

		cublasSaxpy(M, gamma, devS, 1, devZ, 1);

		cublasGetVector(M, sizeof(float), devZ, 1, Z.elements, 1);

		iteration++;
	}

	printf("Solving finished for %f sec\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
	Z *= -1.0f; 
	GridData out(gp);
	memcpy(out.data, Z.elements, M * sizeof(float));

	CheckCublas(cublasFree(devC));
	CheckCublas(cublasFree(devS));
	CheckCublas(cublasFree(devZ));
	CheckCublas(cublasFree(devA));
	CheckCublas(cublasFree(devF));

	return out;
}

//GridData LightLevenbergMarkvardt(Task t, float alpha)
//{
//	GridParameters gp = t.grid.GetGridParameters();
//	int M = gp.NX * gp.NY;
//
//	Matrix A(M, 1), F(M, 1), Z(M, 1), Ank(M, M), Tmp1(M, 1), Tmp2(M, 1);
//	Z.Fill(t.initialZ);
//	t.grid.FillMatrix(F);
//
//	float error;
//	int iteration = 1;
//
//	CudaDirectSolver dslvr(gp, true);
//	CudaOperatorReplacer oper(t);
//
//	float *devF;
//	float *devA;
//	float *devZ;
//	float *devTmp1;
//	float *devTmp2;
//	float *devTmp3;
//					
//	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devF));
//	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devA));
//	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devZ));
//	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devTmp1));
//	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devTmp2));
//	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devTmp3));
//	
//	CheckCublas(cublasSetVector(M, sizeof(float), F.elements, 1, devF, 1));
//	CheckCublas(cublasSetVector(M, sizeof(float), Z.elements, 1, devZ, 1));	
//
//	float beta = 0.01f;
//	float gamma = 0.5f;
//
//	printf("Solving with Levenberg-Markvardt Method\n");
//
//	clock_t t0 = clock();
//
//	while (true)
//	{
//		dslvr.SolveDirectTask(Z, A, Ank, t);
//		
//		CheckCublas(cublasSetVector(M, sizeof(float), A.elements, 1, devA, 1));
//		
//		//A = A - F;
//		//A = A - F
//		cublasSaxpy(M, -1.0f, devF, 1, devA, 1);
//		CheckCublas(cublasGetError());		
//		
//		//---------------------Ошибка---------------------
//		float err1 = cublasSnrm2(M, devA, 1);
//		CheckCublas(cublasGetError());
//		float err2 = cublasSnrm2(M, devF, 1);
//		CheckCublas(cublasGetError());
//		error = err1 / err2;
//
//		printf("Iteration #%d.\tError = %f.\n", iteration, error);
//
//		if (error < t.precision)
//			break;
//		//---------------------Ошибка---------------------
//
//		//Tmp1 = Ank * A;
//		//Tmp1 = AnkT(A - F)
//		oper.CalcATX(devZ, devA, devTmp1);
//
//		//Tmp *= -gamma;
//
//		cublasSscal(M, -gamma, devTmp1, 1);
//		CheckCublas(cublasGetError());
//			
//		//Z *= beta;
//		//Z += Tmp1;
//		cublasSaxpy(M, beta, devZ, 1, devTmp1, 1);
//		CheckCublas(cublasGetError());
//
//		//Tmp2 = Ank * Z;
//		//Z = Ank.Transpose() * Tmp2;
//		oper.CalcAX(devZ, devZ, devTmp2);
//		oper.CalcATX(devZ, devTmp2, devTmp3);
//
//		cublasSaxpy(M, 1.0f, devTmp3, 1, devTmp1, 1);
//		//Отсюда в devTmp1 вся правая часть
//
//		//******************************Решение системы*************************************************
//		// devTmp2 будет за Z, devTmp1 за правую часть 
//
//
//		//**********************************************************************************************
//
//		cublasGetVector(M, sizeof(float), devZ, 1, Z.elements, 1);
//
//		iteration++;
//	}
//
//	printf("Solving finished for %f sec\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
//	Z *= -1.0f; 
//	GridData out(gp);
//	memcpy(out.data, Z.elements, M * sizeof(float));
//
//	CheckCublas(cublasFree(devTmp3));
//	CheckCublas(cublasFree(devTmp2));
//	CheckCublas(cublasFree(devTmp1));
//	CheckCublas(cublasFree(devZ));
//	CheckCublas(cublasFree(devA));
//	CheckCublas(cublasFree(devF));
//
//	return out;
//}


//void CublasSolveSLE(float *devA, float *devz, float *devb, int M)
//{
//	float *devzpp;
//	float *devrk;
//	float *devArk;
//	float *devzp;
//	float *devpk;
//	float *devApk;
//	float *devTmp;
//	cublasStatus state;
//	float normb = cublasSnrm2(M, devb, 1);
//	
//	state = cublasAlloc(M, sizeof(float), (void**)&devzpp);
//	cublasScopy(M, devz, 1, devzpp, 1);
//	state = cublasGetError();
//	
//	state = cublasAlloc(M, sizeof(float), (void**)&devrk);
//	cublasScopy(M, devb, 1, devrk, 1);
//	state = cublasGetError();
//
//	cublasSgemv('T', M, M, 1.0f, devA, M, devzpp, 1, -1.0f, devrk, 1);
//	state = cublasGetError();
//
//	float normr = cublasSnrm2(M, devrk, 1);
//
//	state = cublasAlloc(M, sizeof(float), (void**)&devArk);
//	cublasSgemv('T', M, M, 1.0f, devA, M, devrk, 1, 0.0f, devArk, 1);
//	state = cublasGetError();
//
//	float d = cublasSdot(M, devArk, 1, devrk, 1);
//	state = cublasGetError();
//
//	state = cublasAlloc(M, sizeof(float), (void**)&devzp);
//	cublasScopy(M, devzpp, 1, devzp, 1);
//	state = cublasGetError();
//
//	cublasSaxpy(M, - (normr * normr / d), devrk, 1, devzp, 1);
//	state = cublasGetError();
//
//	state = cublasAlloc(M, sizeof(float), (void**)&devpk);
//	state = cublasAlloc(M, sizeof(float), (void**)&devApk);
//	state = cublasAlloc(M, sizeof(float), (void**)&devTmp);
//
//	int flag = 1;
//	int iterations = 1;
//
//	while (flag == 1)
//	{
//		cublasScopy(M, devb, 1, devrk, 1);
//		state = cublasGetError();
//		cublasSgemv('T', M, M, 1.0f, devA, M, devzp, 1, -1.0f, devrk, 1);
//		state = cublasGetError();
//
//		normr = cublasSnrm2(M, devrk, 1);
//		state = cublasGetError();
//
//		cublasScopy(M, devzp, 1, devpk, 1);
//		state = cublasGetError();
//
//		cublasSaxpy(M, -1.0f, devzpp, 1, devpk, 1);
//		state = cublasGetError();
//
//		cublasSgemv('T', M, M, 1.0f, devA, M, devrk, 1, 0.0f, devArk, 1);
//		state = cublasGetError();
//		cublasSgemv('T', M, M, 1.0f, devA, M, devpk, 1, 0.0f, devApk, 1);
//		state = cublasGetError();
//
//		float dot1 = cublasSdot(M, devArk, 1, devpk, 1);
//		state = cublasGetError();
//		float dot2 = cublasSdot(M, devrk, 1, devpk, 1);
//		state = cublasGetError();
//		float dot3 = cublasSdot(M, devArk, 1, devrk, 1);
//		state = cublasGetError();
//		float dot4 = cublasSdot(M, devApk, 1, devpk, 1);
//		state = cublasGetError();
//
//		d = dot3 * dot4 - dot1 * dot1;
//
//		float gamma = ((normr * normr) * dot4 - dot2 * dot1) / d;
//		float beta = ((normr * normr) * dot1 - dot2 * dot3) / d;
//
//		cublasScopy(M, devzp, 1, devzpp, 1);
//		state = cublasGetError();
//
//		cublasSaxpy(M, -gamma, devrk, 1, devzp, 1);
//		state = cublasGetError();
//		cublasSaxpy(M, beta, devpk, 1, devzp, 1);
//		state = cublasGetError();
//
//		cublasScopy(M, devb, 1, devTmp, 1);
//		state = cublasGetError();
//
//		cublasSgemv('T', M, M, 1.0f, devA, M, devzp, 1, -1.0f, devTmp, 1);
//		state = cublasGetError();
//
//		double norm = cublasSnrm2(M, devTmp, 1);
//		state = cublasGetError();
//
//		double error = norm / normb;
//
//		printf("   Iteration:%d\terror:%f\n", iterations, error);
//
//		if (error < 0.001)
//			flag = 0;
//
//		iterations++;
//	}
//	
//	cublasFree(devzp);
//	cublasFree(devzpp);
//	cublasFree(devArk);
//	cublasFree(devApk);
//	cublasFree(devrk);
//	cublasFree(devpk);
//	cublasFree(devTmp);
//
//	return;
//}

GridData LightLevenbergMarkvardt(Task t, float alpha)
{
	GridParameters gp = t.grid.GetGridParameters();
	DirectSolver dslvr(gp);
	int M = gp.NX * gp.NY;
	int i, j, l, k;
	Matrix F(M, 1), A(M, 1), Z(M, 1), TMP1(M, 1), TMP2(M, 1);
	Matrix Ank(M, M), Bnk(M, M);
	Matrix exactZ(M, 1);
	exactZ *= -1.0f;

	t.grid.FillMatrix(F);
	F *= -1.0f;
	t.exactSolution.FillMatrix(exactZ);
	Z.Fill(t.initialZ);
		
	cublasStatus state;
	float *devA;
	float *devAnk;
	float *devF;
	float *devZ;
	float *devTmpV1;
	float *devB;
					
	state = cublasAlloc(M, sizeof(float), (void**)&devA);
	state = cublasAlloc(M * M, sizeof(float), (void**)&devAnk);
	state = cublasAlloc(M, sizeof(float), (void**)&devF);
	state = cublasAlloc(M, sizeof(float), (void**)&devZ);
	state = cublasAlloc(M, sizeof(float), (void**)&devTmpV1);
	state = cublasAlloc(M * M, sizeof(float), (void**)&devB);

	float beta = 0.01f;
	float gamma = 0.5f;
		
	int iteration = 1;
		

	while (true)
	{
		dslvr.SolveDirectTask(Z, A, Ank, t);
		//A *= -1.0f;
		state = cublasSetVector(M, sizeof(float), A.elements, 1, devA, 1);
		state = cublasSetMatrix(M, M, sizeof(float), Ank.elements, M, devAnk, M);
		state = cublasSetVector(M, sizeof(float), Z.elements, 1, devZ, 1); 
		state = cublasSetVector(M, sizeof(float), F.elements, 1, devF, 1);
		float err = (Z + exactZ).Norm() / exactZ.Norm();
		float u = (A + F).Norm() / F.Norm();
		cout << "Levenberg-Markvardt error = " << u << " 2 = " << err << endl;
		if (u < t.precision)
		{
			break;
		}
		//A - F
		cublasScopy(M, devA, 1, devTmpV1, 1);
		state = cublasGetError();
			
		cublasSaxpy(M, 1.0f, devF, 1, devTmpV1, 1);
		state = cublasGetError();

		cublasScopy(M, devTmpV1, 1, devF, 1);
		state = cublasGetError();

		cublasSgemv('N', M, M, gamma, devAnk, M, devF, 1, 0.0f, devTmpV1, 1);
		state = cublasGetError();
		cublasScopy(M, devTmpV1, 1, devF, 1);
		state = cublasGetError();
			
		cublasSgemm('N', 'T', M, M, M, 1.0f, devAnk, M, devAnk, M, 0.0f, devB, M);
		state = cublasGetError();

		cublasSaxpy(M, beta, devZ, 1, devF, 1);
		state = cublasGetError();

		cublasSgemv('T', M, M, 1.0f, devB, M, devZ, 1, 1.0f, devF, 1);
		state = cublasGetError();

		state = cublasGetVector(M, sizeof(float), devF, 1, TMP2.elements, 1);
		state = cublasGetMatrix(M, M, sizeof(float), devB, M, Ank.elements, M);
		Ank.AddToDiagonal(beta);

		SolveSLE(Ank, Z, TMP2, t.initialZ);

		iteration++;
	}

	Z *= -1.0f;

	GridData result(gp);
	memcpy(result.data, Z.elements, M * sizeof(float));

	return result;
}

__global__ void Dempfers(float *Z, float *dempf, int N)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= N)
		return;

	Z[index] *= dempf[index];
}

vector<GridData> LightMultilayerLinearisedMinimalError(MultilayerTask mTask, float alpha)
{
	GridParameters gp = mTask.GetGeneralGridParameters();
	int M = gp.NX * gp.NY;
	int L = mTask.GetLayersCount();

	Matrix A(M, 1), F(M, 1), TMP1(L * M, 1), Z(L * M, 1), Ze(L * M, 1);
	mTask.InitZ(Z);
	mTask.GetGeneralField().FillMatrix(F);

	float error;
	int iteration = 1;

	CudaDirectSolver dslvr(gp, false);
	CudaOperatorReplacer oper(mTask);

	float *devF;
	float *devA;
	float *devZ;
	float *devS;
	float *dempf;
					
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devF));
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devA));
	CheckCublas(cublasAlloc(L * M, sizeof(float), (void**)&devZ));
	CheckCublas(cublasAlloc(L * M, sizeof(float), (void**)&devS));
	CheckCublas(cublasAlloc(L * M, sizeof(float), (void**)&dempf));

	Matrix Demp(L * M, 1);
	for (int i = 0; i < L; i++)
	{
		memcpy(&Demp.elements[i * M], mTask[i].grid.data, M * sizeof(float));
		memcpy(&Ze.elements[i * M], mTask[i].exactSolution.data, M * sizeof(float));
	}

	float maxVal = -10000000.0f;
	for(int i = 0; i < L * M; i++)
	{
		Demp.elements[i] = pow(abs(Demp.elements[i]), 1.2f); 
		if (abs(Demp.elements[i]) > maxVal)
		{
			maxVal = Demp.elements[i];
		}
	}

	for(int i = 0; i < L * M; i++)
	{
		Demp.elements[i] /= maxVal;
	}
	
	//Нужно только раз
	CheckCublas(cublasSetVector(M, sizeof(float), F.elements, 1, devF, 1));
	CheckCublas(cublasSetVector(L * M, sizeof(float), Demp.elements, 1, dempf, 1));
	CheckCublas(cublasSetVector(L * M, sizeof(float), Z.elements, 1, devZ, 1));
		
	printf("Solving with Linearised Minimal Error Method\n");

	clock_t t0 = clock();
	

	while (true)
	{
		dslvr.SolveDirectMultilayerTask(mTask, Z, A);
		
		CheckCublas(cublasSetVector(M, sizeof(float), A.elements, 1, devA, 1));
		
		// A = A - F
		cublasSaxpy(M, -1.0f, devF, 1, devA, 1);
		CheckCublas(cublasGetError());

		// ||(A - F)||
		float n1 = cublasSnrm2(M, devA, 1);
		
		//---------------------Ошибка---------------------
		float err2 = cublasSnrm2(M, devF, 1);
		CheckCublas(cublasGetError());
		error = n1 / err2;

		printf("Iteration #%d.\tError = %f. Zerror = %f\n", iteration, error, (Ze - Z).Norm() / Ze.Norm());

		if (error < mTask[0].precision)
			break;
		//---------------------Ошибка---------------------

		// S = AnkT(A - F)
		oper.CalcATX(devZ, devA, devS);
		
		// ||S||
		float n2 = cublasSnrm2(L * M, devS, 1);

		float gamma = (alpha * n1 * n1) / (n2 * n2);
		
		int block_size = 256;
		int grid_size = L * M / block_size + (L * M % block_size != 0 ? 1 : 0);
		Dempfers<<<grid_size, block_size>>>(devS, dempf, L * M);
		CheckCuda(cudaGetLastError());
		
		cublasSaxpy(L * M, gamma, devS, 1, devZ, 1);

		CheckCublas(cublasGetVector(L * M, sizeof(float), devZ, 1, Z.elements, 1));

		iteration++;						
	}

	printf("Solving finished for %f sec\n", (double)(clock() - t0) / CLOCKS_PER_SEC);

	vector<GridData> result;
	Z *= -1.0;
	for (int i = 0; i < L; i++)
	{
		GridData layer(gp);
		memcpy(layer.data, &Z.elements[i * M], M * sizeof(float));
		result.push_back(layer);
	}

	CheckCublas(cublasFree(dempf));
	CheckCublas(cublasFree(devS));
	CheckCublas(cublasFree(devZ));
	CheckCublas(cublasFree(devA));
	CheckCublas(cublasFree(devF));

	return result;
}

vector<GridData> LightMultilayerLinearisedSpeedDescent(MultilayerTask mTask, float alpha)
{
	GridParameters gp = mTask.GetGeneralGridParameters();
	int M = gp.NX * gp.NY;
	int L = mTask.GetLayersCount();

	Matrix A(M, 1), F(M, 1), Z(L * M, 1), Ze(L * M, 1);
	mTask.InitZ(Z);
	mTask.GetGeneralField().FillMatrix(F);

	float error;
	int iteration = 1;

	CudaDirectSolver dslvr(gp, false);
	CudaOperatorReplacer oper(mTask);

	float *devF;
	float *devA;
	float *devZ;
	float *devS;
	float *devC;
	float *dempf;
					
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devF));
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devA));
	CheckCublas(cublasAlloc(L * M, sizeof(float), (void**)&devZ));
	CheckCublas(cublasAlloc(L * M, sizeof(float), (void**)&devS));
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devC));
	CheckCublas(cublasAlloc(L * M, sizeof(float), (void**)&dempf));

	Matrix Demp(L * M, 1);
	for (int i = 0; i < L; i++)
	{
		memcpy(&Demp.elements[i * M], mTask[i].grid.data, M * sizeof(float));
		memcpy(&Ze.elements[i * M], mTask[i].exactSolution.data, M * sizeof(float));
	}

	float maxVal = -10000000.0f;
	for(int i = 0; i < L * M; i++)
	{
		Demp.elements[i] = pow(abs(Demp.elements[i]), 1.2f); 
		if (abs(Demp.elements[i]) > maxVal)
		{
			maxVal = Demp.elements[i];
		}
	}

	for(int i = 0; i < L * M; i++)
	{
		Demp.elements[i] /= maxVal;
	}
	
	//Нужно только раз
	CheckCublas(cublasSetVector(M, sizeof(float), F.elements, 1, devF, 1));
	CheckCublas(cublasSetVector(L * M, sizeof(float), Demp.elements, 1, dempf, 1));
	CheckCublas(cublasSetVector(L * M, sizeof(float), Z.elements, 1, devZ, 1));
		
	printf("Solving with Linearised Speed Descent Method\n");

	clock_t t0 = clock();
	
	while (true)
	{
		dslvr.SolveDirectMultilayerTask(mTask, Z, A);

		CheckCublas(cublasSetVector(M, sizeof(float), A.elements, 1, devA, 1));
		
		// A = A - F
		cublasSaxpy(M, -1.0f, devF, 1, devA, 1);
		CheckCublas(cublasGetError());

		//---------------------Ошибка---------------------
		float err1 = cublasSnrm2(M, devA, 1);
		CheckCublas(cublasGetError());
		float err2 = cublasSnrm2(M, devF, 1);
		CheckCublas(cublasGetError());
		error = err1 / err2;
		
		printf("Iteration #%d.\tError = %f. Zerror = %f\n", iteration, error, (Ze - Z).Norm() / Ze.Norm());

		if (error < mTask[0].precision)
			break;
		//---------------------Ошибка---------------------
		
		// S = AnkT(A - F)
		oper.CalcATX(devZ, devA, devS);
		
		// ||S||
		float n1 = cublasSnrm2(L * M, devS, 1);

		// C = Ank * S
		oper.CalcAX(devZ, devS, devC);

		// ||C||
		float n2 = cublasSnrm2(M, devC, 1);

		float gamma = (alpha * n1 * n1) / (n2 * n2);
		
		int block_size = 256;
		int grid_size = L * M / block_size + (L * M % block_size != 0 ? 1 : 0);
		Dempfers<<<grid_size, block_size>>>(devS, dempf, L * M);
		CheckCuda(cudaGetLastError());
		
		cublasSaxpy(L * M, gamma, devS, 1, devZ, 1);

		CheckCublas(cublasGetVector(L * M, sizeof(float), devZ, 1, Z.elements, 1));

		iteration++;					
	}

	printf("Solving finished for %f sec\n", (double)(clock() - t0) / CLOCKS_PER_SEC);

	vector<GridData> result;
	Z *= -1.0;
	for (int i = 0; i < L; i++)
	{
		GridData layer(gp);
		memcpy(layer.data, &Z.elements[i * M], M * sizeof(float));
		result.push_back(layer);
	}

	CheckCublas(cublasFree(dempf));
	CheckCublas(cublasFree(devS));
	CheckCublas(cublasFree(devZ));
	CheckCublas(cublasFree(devA));
	CheckCublas(cublasFree(devF));

	return result;
}

vector<GridData> LightMultilayerLevenbergMarkvardt(MultilayerTask mTask, float alpha)
{
	//Definitions
	GridParameters gp = mTask.GetGeneralGridParameters();

	int L = mTask.GetLayersCount();
	int M = gp.NX * gp.NY;

	Matrix Ank(L * M, M), Ank1(M,L * M), B(M, M);
	Matrix A(M, 1), F(M, 1), TMP1(M, 1), TMP2(M, 1);
	Matrix Z(L * M, 1);

	//Initialising of Z vector
	for (int i = 0; i < L; i++)
	{
		for(int j = 0; j < M; j++)
		{
			Z.elements[i * M + j] = mTask[i].initialZ;
		}
	}

	//Form dempfer values
	Matrix Demp(L * M, 1);
	for (int i = 0; i < L; i++)
	{
		memcpy(&Demp.elements[i * M], mTask[i].grid.data, M * sizeof(float));
	}

	float maxVal = -10000000.0f;
	for(int i = 0; i < L * M; i++)
	{
		if (abs(Demp.elements[i]) > maxVal)
		{
			maxVal = Demp.elements[i];
		}
	}

	for(int i = 0; i < L * M; i++)
	{
		Demp.elements[i] /= maxVal;
	}

	//Initalising of general field matrix
	mTask.GetGeneralField().FillMatrix(F);
	//F *= -1.0f;
	
	//Error value declaration
	float precision = 100000000.0f;
	for (int i = 0; i < L; i++)
	{
		if (mTask[i].precision < 0)
		{
			throw new string("Invalid precision value");
		}

		if (mTask[i].precision < precision)
		{
			precision = mTask[i].precision;
		}
	}
		
	float *devF;
	float *devA;
	float *devZ;
	float *devTmp1;
	float *devTmp2;
	float *devAnk;
					
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devF));
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devA));
	CheckCublas(cublasAlloc(L * M, sizeof(float), (void**)&devZ));
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devTmp1));
	CheckCublas(cublasAlloc(M, sizeof(float), (void**)&devTmp2));
	CheckCublas(cublasAlloc(M * M, sizeof(float), (void**)&devAnk));

	int iteration = 1;
	clock_t t0 = clock();
	float beta = 0.01f;
	float gamma = 0.5f;

	int grid_size;
	int block_size;
	CudaDirectSolver dslvr(mTask.GetGeneralGridParameters(), true);
	
	printf("Solving multilayer task with Levenberg-Markvardt Method\n");

	while (true)
	{
		//dslvr.SolveDirectMultilayerTask(mTask, Z, A, Ank);
		dslvr.SolveDirectTask(Z, A, Ank, mTask[0]);

		CheckCublas(cublasSetVector(M, sizeof(float), A.elements, 1, devA, 1));
		CheckCublas(cublasSetVector(M, sizeof(float), F.elements, 1, devF, 1));
		CheckCublas(cublasSetVector(M, sizeof(float), Z.elements, 1, devZ, 1));
		CheckCublas(cublasSetMatrix(M, M, sizeof(float), Ank.elements, M, devAnk, M));

		//A = A - F;
		cublasSaxpy(M, -1.0f, devF, 1, devA, 1);
		CheckCublas(cublasGetError());		

		//Tmp1 = AT * A
		block_size = 32;
		grid_size = M / block_size + (M % block_size != 0 ? 1 : 0);
		//KernelCalcAZ<<<block_size, grid_size>>> (M, M, devA, devTmp1, devZ, mTask[0].geltaSigm, mTask[0].asimptHeight, gp.dX, gp.dY, gp.NX, gp.NY);


		float u = (A + F).Norm() / F.Norm();
		cout << "Iteration: " << iteration << " Error = " << u << endl;
		if (u < 0.001f || iteration == 7)
		{
			break;
		}
		iteration++;
	}

	printf("Solving finished for %f sec\n", (double)(clock() - t0) / CLOCKS_PER_SEC);

	vector<GridData> result;

	Z *= -1.0;

	for (int i = 0; i < L; i++)
	{
		GridData layer(gp);
		memcpy(layer.data, &Z.elements[i * M], M * sizeof(float));
		result.push_back(layer);
	}

	CheckCublas(cublasFree(devAnk));
	CheckCublas(cublasFree(devTmp2));
	CheckCublas(cublasFree(devTmp1));
	CheckCublas(cublasFree(devZ));
	CheckCublas(cublasFree(devA));
	CheckCublas(cublasFree(devF));

	return result;
}

