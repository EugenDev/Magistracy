#include "DirectSolverBase.hpp"
#define GRAVY_CONST 6.673f

DirectSolverBase::DirectSolverBase()
{
	gridParams.NX = gridParams.NY = 0;
	gridParams.MinX = gridParams.MinY =
		gridParams.dX = gridParams.dY = 0.f;
}

DirectSolverBase::DirectSolverBase(GridParameters params)
{
	this->gridParams = params;
	this->coeff = params.dX * params.dY;
	return;
}

GridData DirectSolverBase::SolveDirectTask(Task task)
{
	GridData out(this->gridParams);

	GridParameters params = task.grid.GetGridParameters();

	if (params.dX != gridParams.dX || params.dY != gridParams.dY ||
		params.NX != gridParams.NX || params.NY != gridParams.NY)
		throw "Direct solver: Inapropriate grid\n";

	switch (task.taskType)
	{
		case TASK_TYPE_GRAVIMETRY:
			SolveDirectGP(task.grid.data, out.data, task.geltaSigm, task.asimptHeight);
			break;

		case TASK_TYPE_MAGNITOMETRY:
			SolveDirectGP(task.grid.data, out.data, task.geltaJ, task.asimptHeight);
			break;

		default:
			break;
	}

	return out;
}

void DirectSolverBase::SolveDirectTask(const Matrix &Z, Matrix &ans, Matrix &opDer, Task task)
{
	if (gridParams.NX * gridParams.NY != Z.GetColsCount() * Z.GetRowsCount())
		throw("Direct solver: Inapropriate matrix\n");

	switch (task.taskType)
	{
		case TASK_TYPE_GRAVIMETRY:
			SolveDirectGP(Z.elements, ans.elements, opDer.elements, task.geltaSigm, task.asimptHeight);
			break;

		case TASK_TYPE_MAGNITOMETRY:
			SolveDirectMP(Z.elements, ans.elements, opDer.elements, task.geltaJ, task.asimptHeight);
			break;

		default:
			throw "Unknown task type\n";
			break;
	}

	return;
}

void DirectSolverBase::SolveDirectTask(const Matrix &Z, Matrix &ans, Task task)
{
	if (gridParams.NX * gridParams.NY != Z.GetColsCount() * Z.GetRowsCount())
		throw("Direct solver: Inapropriate matrix\n");

	switch (task.taskType)
	{
		case TASK_TYPE_GRAVIMETRY:
			SolveDirectGP(Z.elements, ans.elements, task.geltaSigm, task.asimptHeight);
			break;

		case TASK_TYPE_MAGNITOMETRY:
			SolveDirectMP(Z.elements, ans.elements, task.geltaJ, task.asimptHeight);
			break;

		default:
			throw "Unknown task type\n";
			break;
	}

	return;
}

void DirectSolverBase::SolveDirectMultilayerTask(MultilayerTask mTask, Matrix &Z,Matrix &A, Matrix &Ank)
{
	int L = mTask.GetLayersCount();
	int M = this->gridParams.NX * this->gridParams.NY;

	if (A.GetColsCount() * A.GetRowsCount() != M)
	{
		throw "Incorrect out matrix";
	}

	float *accumulator = new float[M];
	memset(accumulator, 0, M * sizeof(float));
	memset(A.elements, 0, M * sizeof(float));
	memset(Ank.elements, 0, M * M * L * sizeof(float));

	for (int i = 0; i < L; i++)
	{
		//TODO: Магнитометрия
		this->SolveDirectGP(&Z.elements[i * M], accumulator, &Ank.elements[M * M * i], 
			mTask[i].geltaSigm, mTask[i].asimptHeight);

		for (int j = 0; j < M; j++)
		{
			A.elements[j] += accumulator[j];
		}
	}

	delete [] accumulator;

	return;
}

void DirectSolverBase::SolveDirectMultilayerTask(MultilayerTask mTask, Matrix &Z,Matrix &A)
{
	int L = mTask.GetLayersCount();
	int M = this->gridParams.NX * this->gridParams.NY;

	if (A.GetColsCount() * A.GetRowsCount() != M)
	{
		throw "Incorrect out matrix";
	}

	float *accumulator = new float[M];
	memset(accumulator, 0, M * sizeof(float));
	memset(A.elements, 0, M * sizeof(float));

	for (int i = 0; i < L; i++)
	{
		//TODO: Магнитометрия
		this->SolveDirectGP(&Z.elements[i * M], accumulator, mTask[i].geltaSigm, mTask[i].asimptHeight);

		for (int j = 0; j < M; j++)
		{
			A.elements[j] += accumulator[j];
		}
	}

	delete [] accumulator;

	return;
}

GridData DirectSolverBase::SolveDirectMultilayerTask(MultilayerTask mTask)
{
	GridParameters gp = mTask.GetGeneralGridParameters();
	int L = mTask.GetLayersCount();
	int M = gp.NX * gp.NY;
	Matrix outMatrix(M, 1);
	Matrix input(M * L, 1);
	GridData result(gp);

	for(int i = 0; i < L; i++)
	{
		memcpy(&input.elements[i * M], mTask[i].grid.data, M * sizeof(float));
	}

	this->SolveDirectMultilayerTask(mTask, input, outMatrix);
	memcpy(result.data, outMatrix.elements, M * sizeof(float));

	return result;
}
