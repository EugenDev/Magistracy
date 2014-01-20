#pragma once
#include "DirectSolver.hpp"

class CudaDirectSolver : public DirectSolverBase
{
private:
	//Indicates if memory for operator is allocated on GPU
	bool isOperatorNeeded;

	//Pointers to GPU memory
	float *dp_inputData;
	float *dp_outputData;
	float *dp_opOutData;

	//Cuda grid parameters
	int block;//dim3 block;
	int grid;//dim3 grid;

	//Methods for direct task solving
	void SolveDirectGP(float *in, float *out, float dSigm, float H);
	void SolveDirectGP(float *in, float *out, float *opOut, float dSigm, float H);

	//Method for calculating optimal grid and block size
	//I'm not sure i cane code it? but it will be here
	void CalculateCudaGrid();

public:
	CudaDirectSolver(GridParameters gp, bool isOperatorNeeded);
	~CudaDirectSolver();
};