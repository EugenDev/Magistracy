#include "CudaDirectSolver.hpp"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#define GRAVY_CONST 6.673f

__global__ void OpCudaDirectSolverKernel(float *in, float *out, float *opOut, float dSigm, float H, float dX, float dY, int NX, int NY)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;

	int m = NX * NY;

	if (index >= m)
		return;
		
	float coeff = dX * dY * dSigm * GRAVY_CONST;

	float xs, ys, zs, d1, d2, sum;
	int k, l;

	l = index / NX;
	k = index % NX;

	sum = 0.0;
	for (int j = 0; j < NY; j++)
	{
		for(int i = 0; i < NX; i++)
		{
			xs = (k - i) * dX;
			ys = (l - j) * dY;
			zs = in[j * NX + i] * in[j * NX + i];
			d1 = sqrt( xs * xs + ys * ys + zs );
			d2 = sqrt( xs * xs + ys * ys + H * H );
			sum += (1.0f / d1 - 1.0f / d2);
			opOut[(l * NX + k) * m + j * NX + i] = in[j * NX + i] * coeff / (d1 * d1 * d1);
		}
	}
	out[index] = sum * coeff;

	return;
}

void CudaDirectSolver::SolveDirectGP(float *in, float *out, float *opOut, float dSigm, float H)
{
	if (!this->isOperatorNeeded)
	{
		throw new string("Can't solve direct task because memory for operator is not allocated");
	}
	
	int m = gridParams.NX * gridParams.NY;

	cudaMemcpy(dp_inputData, in, m * sizeof(float), cudaMemcpyHostToDevice);

	OpCudaDirectSolverKernel<<<grid, block>>>(
		dp_inputData, 
		dp_outputData, 
		dp_opOutData, 
		dSigm, 
		H, 
		gridParams.dX, 
		gridParams.dY,
		gridParams.NX,
		gridParams.NY);
	
	cudaError_t err = cudaGetLastError();

	cudaMemcpy(out, dp_outputData, m * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(opOut, dp_opOutData, m * m * sizeof(float), cudaMemcpyDeviceToHost);
	
	return;
}

__global__ void CudaDirectSolverKernel(float *in, float *out, float dSigm, float H, float dX, float dY, int NX, int NY)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;

	int m = NX * NY;

	if (index >= m)
		return;
	// ���������� �� ���������
	float coeff = dX * dY * dSigm * GRAVY_CONST;

	float xs, ys, zs, d1, d2, sum;
	int k, l;

	l = index / NX;
	k = index % NX;

	sum = 0.0;
	for (int j = 0; j < NY; j++)
	{
		for(int i = 0; i < NX; i++)
		{
			xs = (k - i) * dX;
			ys = (l - j) * dY;
			zs = in[j * NX + i] * in[j * NX + i];
			d1 = sqrt( xs * xs + ys * ys + zs );
			d2 = sqrt( xs * xs + ys * ys + H * H );
			sum += (1.0f / d1 - 1.0f / d2);
		}
	}
	out[index] = sum * coeff;

	return;
}

void CudaDirectSolver::SolveDirectGP(float *in, float *out, float dSigm, float H)
{
	int m = gridParams.NX * gridParams.NY;

	CheckCuda(cudaMemcpy(dp_inputData, in, m * sizeof(float), cudaMemcpyHostToDevice));
	
	CudaDirectSolverKernel<<<grid, block>>>(
		dp_inputData, 
		dp_outputData, 
		dSigm, 
		H, 
		gridParams.dX, 
		gridParams.dY,
		gridParams.NX,
		gridParams.NY);
	
	CheckCuda(cudaGetLastError());

	CheckCuda(cudaMemcpy(out, dp_outputData, m * sizeof(float), cudaMemcpyDeviceToHost));

	return;
}

__global__ void OpMCudaDirectSolverKernel(float *in, float *out, float *opOut, float dJ, float H, float dX, float dY, int NX, int NY)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;

	int m = NX * NY;

	if (index >= m)
		return;
		
	float coeff = dX * dY * dJ;

	float xs, ys, zs, d1, d2, sum;
	int k, l;

	l = index / NX;
	k = index % NX;

	sum = 0.0;
	for (int j = 0; j < NY; j++)
	{
		for(int i = 0; i < NX; i++)
		{
			xs = (k - i) * dX;
			ys = (l - j) * dY;
			zs = in[j * NX + i] * in[j * NX + i];
			d1 = sqrt( xs * xs + ys * ys + zs );
			d2 = sqrt( xs * xs + ys * ys + H * H );
			sum += (in[j * NX + i] / (d1 * d1 * d1) - H / (d2 * d2 * d2));
			opOut[(l * NX + k) * m + j * NX + i] = (1.0f / (d1 * d1 * d1) - 3 * zs / (d1 * d1 * d1 * d1 * d1)) * coeff;
		}
	}
	out[index] = sum * coeff;

	return;
}

void CudaDirectSolver::SolveDirectMP(float *in, float *out, float *opOut, float dJ, float H)
{
	if (!this->isOperatorNeeded)
	{
		throw new string("Can't solve direct task because memory for operator is not allocated");
	}
	
	int m = gridParams.NX * gridParams.NY;

	cudaMemcpy(dp_inputData, in, m * sizeof(float), cudaMemcpyHostToDevice);

	OpMCudaDirectSolverKernel<<<grid, block>>>(
		dp_inputData, 
		dp_outputData, 
		dp_opOutData, 
		dJ, 
		H, 
		gridParams.dX, 
		gridParams.dY,
		gridParams.NX,
		gridParams.NY);
	
	cudaError_t err = cudaGetLastError();

	cudaMemcpy(out, dp_outputData, m * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(opOut, dp_opOutData, m * m * sizeof(float), cudaMemcpyDeviceToHost);
	
	return;
}

__global__ void CudaMDirectSolverKernel(float *in, float *out, float dJ, float H, float dX, float dY, int NX, int NY)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;

	int m = NX * NY;

	if (index >= m)
		return;

	float xs, ys, zs, d1, d2, sum;
	int k, l;

	l = index / NX;
	k = index % NX;

	sum = 0.0;
	for (int j = 0; j < NY; j++)
	{
		for(int i = 0; i < NX; i++)
		{
			xs = (k - i) * dX;
			ys = (l - j) * dY;
			zs = in[j * NX + i] * in[j * NX + i];
			d1 = sqrt( xs * xs + ys * ys + zs );
			d2 = sqrt( xs * xs + ys * ys + H * H );
			sum += (in[j * NX + i] / (d1 * d1 * d1) - H / (d2 * d2 * d2));
		}
	}
	out[index] = sum * dX * dY * dJ;

	return;
}

void CudaDirectSolver::SolveDirectMP(float *in, float *out, float dJ, float H)
{
	int m = gridParams.NX * gridParams.NY;

	cudaMemcpy(dp_inputData, in, m * sizeof(float), cudaMemcpyHostToDevice);

	CudaMDirectSolverKernel<<<grid, block>>>(
		dp_inputData, 
		dp_outputData, 
		dJ, 
		H, 
		gridParams.dX, 
		gridParams.dY,
		gridParams.NX,
		gridParams.NY);
	
	cudaError_t err = cudaGetLastError();

	cudaMemcpy(out, dp_outputData, m * sizeof(float), cudaMemcpyDeviceToHost);

	return;
}
//TODO: ��� ���-��?!!?!?!
void CudaDirectSolver::CalculateCudaGrid()
{
	int m = gridParams.NX * gridParams.NY;

	block.x = 256;
	grid.x = m / block.x + (m % block.x != 0 ? 1 : 0);

	return;
}

CudaDirectSolver::CudaDirectSolver(GridParameters gp, bool operatorNeeded) : DirectSolverBase(gp)
{
	this->isOperatorNeeded = operatorNeeded;

	int m = gp.NX * gp.NY;

	CheckCuda(cudaMalloc((void**)&dp_inputData,  m * sizeof(float)));
	CheckCuda(cudaMalloc((void**)&dp_outputData, m * sizeof(float)));

	if (isOperatorNeeded)
	{
		CheckCuda(cudaMalloc((void**)&dp_opOutData,  m * m * sizeof(float)));
	}

	this->CalculateCudaGrid();
}

CudaDirectSolver::~CudaDirectSolver()
{
	CheckCuda(cudaFree(dp_inputData));
	CheckCuda(cudaFree(dp_outputData));

	if (isOperatorNeeded)
	{
		CheckCuda(cudaFree(dp_opOutData));
	}
}