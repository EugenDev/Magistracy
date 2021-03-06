#include "ReverseSolver.hpp"

void Makeb(Matrix &Ank, Matrix &Fnk, Matrix &b);
void MakeB(Matrix &Ank, Matrix &mB);
void CheckTasks(vector<Task> tasks);

GridData LinearisedSpeedDescent(DirectSolverBase *directSolver, Task task, float alpha)
{
	GridParameters gp = task.grid.GetGridParameters();
	int NX = gp.NX;
	int NY = gp.NY;
	int m = NX * NY;
	float dX = gp.dX;
	float dY = gp.dY;

	Matrix A(m, 1), Ank(m, m), Z(m, 1), S(m, 1), C(m, 1);
	Z.Fill(task.initialZ);

	Matrix F(m, 1);
	task.grid.FillMatrix(F);

	float error;
	int iteration = 1;

	GridData exactGrid(task.exactSolution);
	Matrix exactSolution(m, 1);
	exactGrid.FillMatrix(exactSolution);
	
	printf("Solving with Linearised Speed Descent Method\n");

	clock_t t0 = clock();

	while (1)
	{
		directSolver->SolveDirectTask(Z, A, Ank, task);
		C = A - F;
		Makeb(Ank, C, S);

		float n1 = S.Norm();
		float n2 = (Ank * S).Norm();
		float gamma = (alpha * n1 * n1) / (n2 * n2);

		S *= gamma;

		Z += S;

		switch (task.residualType)
		{
		case RESIDUAL_TYPE_RIGHTHAND:
			error = (F - A).Norm() / F.Norm();
			break;
		case RESIDUAL_TYPE_EXACTSOLUTION:
			error = (exactSolution + Z).Norm() / exactSolution.Norm();
			break;
		default:
			throw "Unknown residual type";
		}

		printf("Iteration #%d.\tError = %f\n", iteration, error);

		if(error < task.precision)
			break;

		iteration++;
	}

	printf("Solving finished for %f sec\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
	
	Z *= -1.0f;
	GridData out(gp);
	memcpy(out.data, Z.elements, m * sizeof(float));

	return out;
}

GridData LinearisedMinimalError(DirectSolverBase *directSolver, Task task, float alpha)
{
	GridParameters gp = task.grid.GetGridParameters();
	int NX = gp.NX;
	int NY = gp.NY;
	int m = NX * NY;
	float dX = gp.dX;
	float dY = gp.dY;

	Matrix A(m, 1), Ank(m, m), Z(m, 1), S(m, 1);
	Z.Fill(task.initialZ);

	Matrix F(m, 1);
	task.grid.FillMatrix(F);

	float error;
	int iteration = 1;

	Matrix exactSolution(m, 1);
	task.exactSolution.FillMatrix(exactSolution);
	
	printf("Solving with Linearised Minimal Error Method\n");

	clock_t t0 = clock();

	while (true)
	{
		directSolver->SolveDirectTask(Z, A, Ank, task);
		
		Matrix C = F - A;
		Makeb(Ank, C, S);
		float n1 = C.Norm();
		float n2 = S.Norm();
		float gamma = (alpha * n1 * n1) / (n2 * n2);

		S *= gamma;

		Z -= S;

		switch (task.residualType)
		{
			case RESIDUAL_TYPE_RIGHTHAND:
				error = (F - A).Norm() / F.Norm();
				break;
			case RESIDUAL_TYPE_EXACTSOLUTION:
				error = (exactSolution + Z).Norm() / exactSolution.Norm();
				break;
			default:
				throw "Unknown residual type";
		}

		printf("Iteration #%d.\tError = %f.\n", iteration, error);

		if (error < task.precision)
			break;

		iteration++;
	}

	printf("Solving finished for %f sec\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
	Z *= -1.0f; 
	GridData out(gp);
	memcpy(out.data, Z.elements, m * sizeof(float));
	return out;
}

void G_LinearisedSpeedDescentMultilayer()
{
	/*GridData l12_field("grids\\model_data\\gr46x51_layer12_field.dat");
	Task l12_task = {
		l12_field
	};

	GridData layer1_field("grids\\model_data\\gr46x51_layer1_field.dat");
	GridData layer2_field("grids\\model_data\\gr46x51_layer2_field.dat");

	Task l1;
	l1.asimptHeight = 5.0f;
	l1.geltaSigm = 0.25f;
	l1.initialZ = 5.0f;
	l1.taskType = GRAVIMETRY;

	Task l2;
	l2.asimptHeight = 20.0f;
	l2.geltaSigm = 0.3f;
	l2.initialZ = 20.0f;
	l2.taskType = GRAVIMETRY;

	int NX = l12_task.grid.getNX();
	int NY = l12_task.grid.getNY();
	int M = NX * NY;
	float dX = l12_task.grid.getdX();
	float dY = l12_task.grid.getdY();

	Matrix Ank(M * 2, M), Ank1(M, M), Ank2(M, M);
	Matrix A(M, 1), A1(M, 1), A2(M, 1);
	Matrix Z(M * 2, 1), Z1(M, 1), Z2(M, 1);
	Z1.Fill(5.0f);
	Z2.Fill(20.0f);
	
	Matrix L1F(M, 1);
	Matrix L2F(M, 1);
	layer1_field.fillMatrix(L1F);
	layer2_field.fillMatrix(L2F);

	float max = -0x7fffffff;

	for(int i = 0; i < M; i++)
	{
		if (L1F.elements[i] > max)
			max = L1F.elements[i];

		if (L2F.elements[i] > max)
			max = L2F.elements[i];
	}

	for (int i = 0; i < M; i++)
	{
		L1F.elements[i] /= max;
		L2F.elements[i] /= max;
	}

	for(int i = 0; i < M; i++)
		Z.elements[i] = 5.0f;
	for(int i = M; i < 2 * M; i++)
		Z.elements[i] = 20.0f;

	Matrix F(M, 1), S(M, 1);
	l12_field.fillMatrix(F);

	Matrix S1(M, 1);
	Matrix S2(M, 1);

	DirectSolver dSlvr(l12_task.grid.getGridParameters());
	
	float alpha = 0.5;

	printf("Solving multilayer task with Linearised Speed Descent Method\n");
	clock_t t0 = clock();
	int iteration = 1;

	while (true)
	{
		dSlvr.SolveDirectTask(Z1, A1, Ank1, l1);
		dSlvr.SolveDirectTask(Z2, A2, Ank2, l2);

		A = A1 + A2;

		memcpy(Ank.elements, Ank1.elements, M * M * sizeof(float));
		memcpy(&Ank.elements[M * M], Ank2.elements, M * M * sizeof(float));

		S = Ank * (A - F);

		float n1 = (A - F).Norm();
		float n2 = S.Norm();
		float gamma = (alpha * n1 * n1) / (n2 * n2);

		S *= gamma;

		for(int i = 0; i < M; i++)
		{
			S.elements[i] *= L1F.elements[i];
			S.elements[i + M] *= L2F.elements[i];
		}

		Z += S;

		float error = (F - A).Norm() / F.Norm();

		printf("Iteration #%d.\tError = %f\n", iteration, error);
		iteration++;

		memcpy(Z1.elements, Z.elements, M * sizeof(float));
		memcpy(Z2.elements, &Z.elements[M], M * sizeof(float));

		if (error < 0.01)
			break;*/
	/*}

	printf("Solving finished for %f sec\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
	
	Z1 *= -1.0f;
	Z2 *= -1.0f;

	GridData layer1_out(l12_task.grid.getGridParameters());
	GridData layer2_out(l12_task.grid.getGridParameters());

	memcpy(layer1_out.data, Z1.elements, M * sizeof(float));
	memcpy(layer2_out.data, Z2.elements, M * sizeof(float));

	layer1_out.SaveToDATFile(string("grids\\model_data\\gr46x51_layer1_ans.dat"));
	layer2_out.SaveToDATFile(new s"grids\\model_data\\gr46x51_layer2_ans.dat");
	return;*/
}

#include <cublas.h>

void CublasSolveSLE1(float *A, float *Z, float *b, int m)
{
	clock_t time;
	printf(" Solving SLE ... \n");
	time = clock();

	int i;

	float *zpp = new float [m];

	float *devb;
	float *devA;
	float *devzpp;
	float *devrk;
	float *devArk;
	float *devzp;
	float *devpk;
	float *devApk;
	float *devTmp;
	
	cublasAlloc(m, sizeof(float), (void**)&devb);
	cublasSetVector(m, sizeof(float), b, 1, devb, 1);

	float normb = cublasSnrm2(m, devb, 1);

	for (i = 0; i < m; i++)
		zpp[i] = 10.0f;

	cublasAlloc(m, sizeof(float), (void**)&devzpp);
	cublasSetVector(m, sizeof(float), zpp, 1, devzpp, 1);
	
	cublasAlloc(m * m, sizeof(float), (void**)&devA);
	cublasSetMatrix(m, m, sizeof(float), A, m, devA, m);

	cublasAlloc(m, sizeof(float), (void**)&devrk);
	cublasScopy(m, devb, 1, devrk, 1);

	cublasSgemv('T', m, m, 1.0f, devA, m, devzpp, 1, -1.0f, devrk, 1);

	float normr = cublasSnrm2(m, devrk, 1);

	cublasAlloc(m, sizeof(float), (void**)&devArk);
	cublasSgemv('T', m, m, 1.0f, devA, m, devrk, 1, 0.0f, devArk, 1);

	float d = cublasSdot(m, devArk, 1, devrk, 1);

	cublasAlloc(m, sizeof(float), (void**)&devzp);
	cublasScopy(m, devzpp, 1, devzp, 1);

	cublasSaxpy(m, - (normr * normr / d), devrk, 1, devzp, 1);

	cublasAlloc(m, sizeof(float), (void**)&devpk);
	cublasAlloc(m, sizeof(float), (void**)&devApk);
	cublasAlloc(m, sizeof(float), (void**)&devTmp);

	int flag = 1;
	int iterations = 1;

	while (flag == 1)
	{
		cublasScopy(m, devb, 1, devrk, 1);
		cublasSgemv('T', m, m, 1.0f, devA, m, devzp, 1, -1.0f, devrk, 1);

		normr = cublasSnrm2(m, devrk, 1);

		cublasScopy(m, devzp, 1, devpk, 1);

		cublasSaxpy(m, -1.0f, devzpp, 1, devpk, 1);

		cublasSgemv('T', m, m, 1.0f, devA, m, devrk, 1, 0.0f, devArk, 1);
		cublasSgemv('T', m, m, 1.0f, devA, m, devpk, 1, 0.0f, devApk, 1);

		float dot1 = cublasSdot(m, devArk, 1, devpk, 1);
		float dot2 = cublasSdot(m, devrk, 1, devpk, 1);
		float dot3 = cublasSdot(m, devArk, 1, devrk, 1);
		float dot4 = cublasSdot(m, devApk, 1, devpk, 1);

		d = dot3 * dot4 - dot1 * dot1;

		float gamma = ((normr * normr) * dot4 - dot2 * dot1) / d;
		float beta = ((normr * normr) * dot1 - dot2 * dot3) / d;

		cublasScopy(m, devzp, 1, devzpp, 1);

		cublasSaxpy(m, -gamma, devrk, 1, devzp, 1);
		cublasSaxpy(m, beta, devpk, 1, devzp, 1);

		cublasScopy(m, devb, 1, devTmp, 1);

		cublasSgemv('T', m, m, 1.0f, devA, m, devzp, 1, -1.0f, devTmp, 1);

		double norm = cublasSnrm2(m, devTmp, 1);

		double error = norm / normb;

		printf("   Iteration:%d\terror:%f\n", iterations, error);

		if (error < 0.001)
			flag = 0;

		iterations++;
	}
	printf("SLE solved with %d iterations for %f sec\n", iterations, (double)(clock() - time) / CLOCKS_PER_SEC);

	cublasGetVector(m, sizeof(float), devzp, 1, Z, 1);

	cublasFree(devA);
	cublasFree(devb);
	cublasFree(devzp);
	cublasFree(devzpp);
	cublasFree(devArk);
	cublasFree(devApk);
	cublasFree(devrk);
	cublasFree(devpk);
	cublasFree(devTmp);

	return;
}

GridData G_Newton(Task task)
{
	GridParameters gp = task.grid.GetGridParameters();
	int NX = gp.NX;
	int NY = gp.NY;
	int m = NX * NY;
	float dX = gp.dX;
	float dY = gp.dY;

	Matrix A(m, 1), Ank(m, m), Z(m, 1);
	Z.Fill(task.initialZ);

	DirectSolver dSlvr(gp);
	float alpha = 0.05f;
	int iteration = 1;

	while (1)
	{
		cout << "Calculating A and Ank...\n";
		dSlvr.SolveDirectTask(Z, A, Ank, task);

		Ank.AddToDiagonal(alpha);

		if(iteration == 2)
			break;

		Matrix F(m, 1);
		memcpy(F.elements, task.grid.data, m * sizeof(float));

		Matrix Fnk(m, 1);

		cout << "Calculating Fnk...\n";
		Fnk = Ank * Z + Z * alpha - F + A;

		cout << "Calculating b...\n";
		Matrix b(m, 1);
		b.Fill(0.0f);
		Makeb(Ank, Fnk, b);

		cout << "Calculating B...\n";
		Matrix B(m, m);
		B.Fill(0.0f);
		MakeB(Ank, B);

		Matrix tmpV(m, 1);
		tmpV = B * Z - b;
		float newtonError = tmpV.Norm() / b.Norm();

		printf("Newton error on iteration #%d = %f\n", iteration, newtonError);

		if (newtonError < 0.002)
			break;
		
		//CublasSolveSLE1(B.elements, Z.elements, b.elements, m);
		SolveSLE(B, Z, b, 1);
		
		iteration++;
	}
	
	Z *= -1.0f;
	
	GridData out(gp);
	memcpy(out.data, Z.elements, m * sizeof(float)); 
	return out;
}

void MakeB(Matrix &mAnk, Matrix &mB)
{
	float *Ank = mAnk.elements;
	float *B = mB.elements;
	int m = mAnk.GetColsCount();

	int i, j, k;

	for (k = 0; k < m; k++)	
		for (i = 0; i < m; i++)
			for (j = i; j < m; j++)
				B[i * m + j] += Ank[k * m + i] * Ank[k * m + j];

	for (i = 0; i < m; i++)
		for (j = i; j < m; j++)
			B[j * m + i] = B [i * m + j];
}

void Makeb(Matrix &Ank, Matrix &Fnk, Matrix &b)
{
	b.Fill(0.0f);

	float *A = Ank.elements;
	float *F = Fnk.elements;
	float *B = b.elements;
	int m = Ank.GetColsCount();
	int n = Fnk.GetRowsCount();
	
	for (int i = 0 ; i < n; i++)
		for (int k = 0; k < m; k++)
			B[i] += A[k * m + i] * F[k];

	return;
}

vector<GridData> G_GradientMultilayer(MultilayerTask mTask, float alpha)
{	
	//Definitions
	GridParameters gp = mTask.GetGeneralGridParameters();

	int L = mTask.GetLayersCount();
	int M = gp.NX * gp.NY;

	Matrix Ank(L * M, M);
	Matrix A(M, 1), F(M, 1);
	Matrix Z(L * M, 1), S(L * M, 1), C(L * M, 1);

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
	
	//Initialising of direct solver
	CudaDirectSolver dslvr(gp, true);

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

	printf("Solving multilayer task with Linearised Minimal Error Method\n");
	int iteration = 1;
	clock_t t0 = clock();

	while (true)
	{
		dslvr.SolveDirectMultilayerTask(mTask, Z, A, Ank);

		C = A - F;

		S = Ank * C;

		float normS = S.Norm();
		float normC = C.Norm();

		S *= (normC * normC) / (normS * normS);

		for(int i = 0; i < L * M; i++)
		{
			S.elements[i] *= Demp.elements[i] * alpha;
		}

		Z += S;

		float error = (F - A).Norm() / F.Norm();
		printf("Iteration #%d.\tError = %f.\n", iteration++, error);

		if (error < precision)
			break;
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

	return result;
}

void CheckTasks(vector<Task> tasks)
{
	int taskCount = tasks.size();

	if (taskCount < 1)
	{
		throw "No input tasks";
	}

	for (int i = 1; i < taskCount; i++)
	{
		if (tasks[i - 1].grid.GetGridParameters() != tasks[i].grid.GetGridParameters())
		{
			throw "Not all tasks has equal grids";
		}
	}
	
	return;
}

void MakeMagic(Matrix &Ank, Matrix &C, Matrix &S)
{

	return;
}

GridData LevenbergMarkwardt(DirectSolverBase *directSolver, Task task, float alpha)
{
	GridParameters gp = task.grid.GetGridParameters();
	int NX = gp.NX;
	int NY = gp.NY;
	int m = NX * NY;
	float dX = gp.dX;
	float dY = gp.dY;

	Matrix A(m, 1), Ank(m, m), Z(m, 1), S(m, 1);
	Z.Fill(task.initialZ);

	Matrix F(m, 1);
	task.grid.FillMatrix(F);

	float error;
	int iteration = 1;

	Matrix exactSolution(m, 1);
	task.exactSolution.FillMatrix(exactSolution);
	
	printf("Solving with Levenberg-Markvardt Method\n");

	clock_t t0 = clock();

	while (true)
	{
		directSolver->SolveDirectTask(Z, A, Ank, task);
		
		Matrix C = F - A;
		Makeb(Ank, C, S);
		//MakeMagic(Ank, C, S
		float n1 = C.Norm();
		float n2 = S.Norm();
		

		S *= 1.0f;

		Z -= S;

		switch (task.residualType)
		{
			case RESIDUAL_TYPE_RIGHTHAND:
				error = (F - A).Norm() / F.Norm();
				break;
			case RESIDUAL_TYPE_EXACTSOLUTION:
				error = (exactSolution + Z).Norm() / exactSolution.Norm();
				break;
			default:
				throw "Unknown residual type";
		}

		printf("Iteration #%d.\tError = %f.\n", iteration, error);

		if (error < task.precision)
			break;

		iteration++;
	}

	printf("Solving finished for %f sec\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
	Z *= -1.0f; 
	GridData out(gp);
	memcpy(out.data, Z.elements, m * sizeof(float));
	return out;
}
