#include "CudaOperatorReplacer.hpp"

__device__ float GetAElement(int index, float *in, float *dSigmas, float *heights, float dX, float dY, int NX, int NY, int LayersCount)
{
	float xs, ys, d;
	int M = NX * NY;
	unsigned long long row = index / (LayersCount * M);
	unsigned long long col = index % (LayersCount * M);

	int l = col / NX;
	int k = col % NX;
	int j = row / NX;
	int i = row % NX;

	int layer = col / M;

	float z = in[j * NX + i + layer * M];

	xs = (k - i) * dX;
	ys = (l - j) * dY;
	d = sqrt( xs * xs + ys * ys + z * z );
	return z * dX * dY * dSigmas[layer] * GRAVY_CONST / (d * d * d);
}

__device__ float GetATElement1(int index, float *in, float *dSigmas, float *heights, float dX, float dY, int NX, int NY, int LayersCount)
{
	float xs, ys, d;
	int M = NX * NY;
	unsigned long long row = (index / M) % M;
	unsigned long long col = index % M;
	unsigned long long layer = index / (M * M);

	int l = row / NX;
	int k = row % NX;
	int j = col / NX;
	int i = col % NX;

	float z = in[j * NX + i + layer * M];

	xs = (k - i) * dX;
	ys = (l - j) * dY;
	d = sqrt( xs * xs + ys * ys + z * z );
	return z * dX * dY * dSigmas[layer] * GRAVY_CONST / (d * d * d);
}

__device__ float GetATElement(int row, int col, float *in, float *dSigmas, float *heights, float dX, float dY, int NX, int NY, int LayersCount)
{
	float xs, ys, d;
	int M = NX * NY;
	int layer = row / M;

	row = row % M;

	int l = row / NX;
	int k = row % NX;
	int j = col / NX;
	int i = col % NX;

	float z = in[j * NX + i + layer * M];

	xs = (k - i) * dX;
	ys = (l - j) * dY;
	d = sqrt( xs * xs + ys * ys + z * z );
	return z * dX * dY * dSigmas[layer] * GRAVY_CONST / (d * d * d);
}

__global__ void KernelCalcATZ(int nRows, int nCols, float *Z, float *res, float *in, float *dSigmas, float *heights, float dX, float dY, int NX, int NY, int LayersCount)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= nCols)
		return;

	float sum = 0.0f;
	for (int i = 0 ; i < nRows; i++)
	{
		sum += Z[i] * GetATElement(index, i, in, dSigmas, heights, dX, dY, NX, NY, LayersCount);
	}
	res[index] = sum;
 }

__global__ void KernelCalcAZ(int nRows, int nCols, float *Z, float *res, float *in, float *dSigmas, float *heights, float dX, float dY, int NX, int NY, int LayersCount)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= nRows)
		return;

	float sum = 0.0f;
	for (int i = 0 ; i < nCols; i++)
	{
		//xz
		sum += Z[i] * GetATElement(index, i, in, dSigmas, heights, dX, dY, NX, NY, LayersCount);
	}
	res[index] = sum;

	return;
}

__global__ void KernelCalcTMatrix(int nRows, int nCols, float *matrix, float* Z, float *dSigmas, float *heights, float dX, float dY, int NX, int NY, int LayersCount)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x; 
	if (index >= nCols)
		return;

	for(int i = 0; i < nRows; i++)
	{
		matrix[index * nRows + i] = GetATElement(index, i, Z, dSigmas, heights, dX, dY, NX, NY, LayersCount);
	}
}

void CudaOperatorReplacer::CalcAX(float *Z, float *X, float *result)
{
	int block_size = 256;
	int grid_size = this->nRows / block_size + (this->nRows % block_size != 0 ? 1 : 0);
	KernelCalcAZ<<<grid_size, block_size>>>(nRows, nCols, X, result, Z, sigmas, heights, gp.dX, gp.dY, gp.NX, gp.NY, this->layersCount);
	CheckCuda(cudaGetLastError());
	return;
}

void CudaOperatorReplacer::CalcATX(float *Z, float *X, float *result)
{
	int block_size = 256;
	int grid_size = this->nCols / block_size + (this->nCols % block_size != 0 ? 1 : 0);
	KernelCalcATZ<<<grid_size, block_size>>>(nRows, nCols, X, result, Z, sigmas, heights, gp.dX, gp.dY, gp.NX, gp.NY, this->layersCount);
	CheckCuda(cudaGetLastError());
	return;
}

void CudaOperatorReplacer::GetMatrix(Matrix &Res, float* Z)
{
	if (this->nRows != Res.GetRowsCount()
		|| this->nCols != Res.GetColsCount())
	{
		throw new string("Inapropriate matrix");
	}

	
}

void CudaOperatorReplacer::GetTransposedMatrix(Matrix &Res, float* Z)
{
	if (this->nCols != Res.GetRowsCount()
		|| this->nRows != Res.GetColsCount())
	{
		throw new string("Inapropriate matrix");
	}
	float *tmpMat;
	CheckCuda(cudaMalloc((void**)&tmpMat, sizeof(float) * nRows * nCols));

	int block_size = 256;
	int grid_size = this->nCols / block_size + (this->nCols % block_size != 0 ? 1 : 0);
	KernelCalcTMatrix<<<grid_size, block_size>>>(nRows, nCols, tmpMat, Z, sigmas, heights, gp.dX, gp.dY, gp.NX, gp.NY, this->layersCount);
	CheckCuda(cudaGetLastError());

	CheckCuda(cudaMemcpy(Res.elements, tmpMat, sizeof(float) * nRows * nCols, cudaMemcpyDeviceToHost));
	CheckCuda(cudaFree(tmpMat));
}

CudaOperatorReplacer::CudaOperatorReplacer(Task t)
{
	this->isMultilayer = false;
	this->gp = t.grid.GetGridParameters();
	this->task = t;
	nRows = nCols = gp.NX * gp.NY;

	int L = this->layersCount = 1;
	float *heightsHost = new float[L];
	float *sigmasHost = new float[L];
	heightsHost[0] = t.asimptHeight;
	sigmasHost[0] = t.geltaSigm;	
	CheckCuda(cudaMalloc((void**)&heights, L * sizeof(float)));
	CheckCuda(cudaMalloc((void**)&sigmas,  L * sizeof(float)));
	CheckCuda(cudaMemcpy(heights, heightsHost, L * sizeof(float), cudaMemcpyHostToDevice));
	CheckCuda(cudaMemcpy(sigmas, sigmasHost, L * sizeof(float), cudaMemcpyHostToDevice));
	delete [] sigmasHost;
	delete [] heightsHost;
}

CudaOperatorReplacer::CudaOperatorReplacer(MultilayerTask mTask)
{
	this->isMultilayer = true;
	this->gp = mTask.GetGeneralGridParameters();
	this->mTask = &mTask;
	nRows = gp.NX * gp.NY;
	nCols = nRows * mTask.GetLayersCount();

	int L = this->layersCount = mTask.GetLayersCount();
	float *heightsHost = new float[L];
	float *sigmasHost = new float[L];
	for(int l = 0; l < L; l++)
	{
		heightsHost[l] = mTask[l].asimptHeight;
		sigmasHost[l] = mTask[l].geltaSigm;		
	}
	CheckCuda(cudaMalloc((void**)&heights, L * sizeof(float)));
	CheckCuda(cudaMalloc((void**)&sigmas,  L * sizeof(float)));
	CheckCuda(cudaMemcpy(heights, heightsHost, L * sizeof(float), cudaMemcpyHostToDevice));
	CheckCuda(cudaMemcpy(sigmas, sigmasHost, L * sizeof(float), cudaMemcpyHostToDevice));
	delete [] sigmasHost;
	delete [] heightsHost;
}

CudaOperatorReplacer::~CudaOperatorReplacer()
{
	CheckCuda(cudaFree(this->heights));
	CheckCuda(cudaFree(this->sigmas));
}
