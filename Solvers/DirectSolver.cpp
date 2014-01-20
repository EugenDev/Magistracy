#include "DirectSolver.hpp"
#define GRAVY_CONST 6.673f

DirectSolver::DirectSolver(GridParameters params) : DirectSolverBase(params)
{
}

//Вычислительное ядро для прямой задачи гравиметрии без оператора  
void DirectSolver::SolveDirectGP(float *in, float *out, float dSigm, float H)
{
	float sum, xs, ys, zs, d1, d2;
	int NX = gridParams.NX;
	int NY = gridParams.NY;
	float dX = gridParams.dX;
	float dY = gridParams.dY;
	
	int l, k, i, j;
	//TODO: Оптимизировать
	for(l = 0; l < NY; l++)
	{
		for(k = 0; k < NX; k++)
		{
			sum = 0.0;
			for (j = 0; j < NY; j++)
			{
				for(i = 0; i < NX; i++)
				{
					xs = (k - i) * dX;
					ys = (l - j) * dY;
					zs = in[j * NX + i] * in[j * NX + i];
					d1 = sqrt( xs * xs + ys * ys + zs );
					d2 = sqrt( xs * xs + ys * ys + H * H );
					sum += (1.0f / d1 - 1.0f / d2);
				}
			}
			out[l * NX + k] = sum * coeff * dSigm * GRAVY_CONST;
		}
	}

	return;
}

//Вычислительное ядро для прямой задачи гравиметрии с оператором
void DirectSolver::SolveDirectGP(float *in, float *out, float *opOut, float dSigm, float H)
{
	float sum, xs, ys, zs, d1, d2;
	int NX = gridParams.NX;
	int NY = gridParams.NY;
	float dX = gridParams.dX;
	float dY = gridParams.dY;

	int l, k, i, j;
	//TODO: Оптимизировать
	for(l = 0; l < NY; l++)
	{
		for(k = 0; k < NX; k++)
		{
			sum = 0.0;
			for (j = 0; j < NY; j++)
			{
				for(i = 0; i < NX; i++)
				{
					xs = (k - i) * dX;
					ys = (l - j) * dY;
					zs = in[j * NX + i] * in[j * NX + i];
					d1 = sqrt( xs * xs + ys * ys + zs );
					d2 = sqrt( xs * xs + ys * ys + H * H );
					sum += (1.0f / d1 - 1.0f / d2);
					opOut[(l * NX + k) * NX * NY + j * NX + i] = 
						in[j * NX + i] * coeff * dSigm * GRAVY_CONST / (d1 * d1 * d1);
				}
			}
			out[l * NX + k] = sum * coeff * dSigm * GRAVY_CONST;
		}
	}

	return;
}

//Вычислительное ядро для прямой задачи магнитометрии без оператора  
void DirectSolver::SolveDirectMP(float *in, float *out, float dJ, float H)
{
	float sum, xs, ys, xy, zs, d1, d2;
	int NX = gridParams.NX;
	int NY = gridParams.NY;
	float dX = gridParams.dX;
	float dY = gridParams.dY;
	
	int l, k, i, j;
	//TODO: Оптимизировать
	for(l = 0; l < NY; l++)
	{
		for(k = 0; k < NX; k++)
		{
			sum = 0.0;
			for (j = 0; j < NY; j++)
			{
				for(i = 0; i < NX; i++)
				{
					/*xs = (k - i) * dX;
					ys = (l - j) * dY;
					zs = in[j * NX + i] * in[j * NX + i];
					d1 = sqrt( xs * xs + ys * ys + zs );
					d2 = sqrt( xs * xs + ys * ys + H * H );
					sum += (1.0f / d1 - 1.0f / d2);*/
					xs = (k - i) * dX;
					ys = (l - j) * dY;
					xy = xs*xs + ys*ys;
					zs = in[j * NX + i] * in[j * NX + i];
					d1 = sqrt( xy + zs );
					d2 = sqrt( xy + H*H );
					sum = sum + (in[j * NX + i] / (d1 * d1 * d1) - H / (d2 * d2 * d2));
				}
			}
			out[l * NX + k] = sum * dX * dY * dJ;
		}
	}

	return;
}

//Вычислительное ядро для прямой задачи магнитометрии с оператором
void DirectSolver::SolveDirectMP(float *in, float *out, float *opOut, float dJ, float H)
{
	float sum, xs, ys, xy, zs, d1, d2;
	int NX = gridParams.NX;
	int NY = gridParams.NY;
	float dX = gridParams.dX;
	float dY = gridParams.dY;

	int l, k, i, j;
	//TODO: Оптимизировать
	for(l = 0; l < NY; l++)
	{
		for(k = 0; k < NX; k++)
		{
			sum = 0.0;
			for (j = 0; j < NY; j++)
			{
				for(i = 0; i < NX; i++)
				{
					/*xs = (k - i) * dX;
					ys = (l - j) * dY;
					zs = in[j * NX + i] * in[j * NX + i];
					d1 = sqrt( xs * xs + ys * ys + zs );
					d2 = sqrt( xs * xs + ys * ys + H * H );
					sum += (1.0f / d1 - 1.0f / d2);
					opOut[(l * NX + k) * NX * NY + j * NX + i] = 
						in[j * NX + i] * coeff * dSigm * GRAVY_CONST / (d1 * d1 * d1);*/
					xs = (k - i) * dX;
					ys = (l - j) * dY;
					xy = xs*xs + ys*ys;
					zs = in[j * NX + i] * in[j * NX + i];
					d1 = sqrt( xy + zs );
					d2 = sqrt( xy + H*H );
					sum = sum + (in[j * NX + i] / (d1 * d1 * d1) - H / (d2 * d2 * d2));
					opOut[(l * NX + k) * NX * NY + j * NX + i] = (1.0f / (d1 * d1 * d1) - 3 * zs / (d1 * d1 * d1 * d1 * d1)) * dX * dY * dJ;
				}
			}
			out[l * NX + k] = sum * dX * dY * dJ;
		}
	}

	return;
}
