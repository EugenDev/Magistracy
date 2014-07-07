#pragma once
#include <GridData.hpp>
#include "DirectSolver.hpp"
#include "CudaDirectSolver.hpp"
#include "DirectSolverBase.hpp"
#include "SLESolver.hpp"
#include <Matrix.hpp>
#include <Task.hpp>
#include <MultilayerTask.hpp>
#include <vector>

GridData LinearisedSpeedDescent(DirectSolverBase *directSolver, Task task, float alpha);
GridData LinearisedMinimalError(DirectSolverBase *directSolver, Task task, float alpha);

GridData M_LinearisedSpeedDescent(Task task, float alpha);
void	 G_LinearisedSpeedDescentMultilayer();
GridData G_Newton(Task task);
vector<GridData> G_GradientMultilayer(MultilayerTask mTask, float alpha);
//vector<GridData> G_LevMarMultilayer(MultilayerTask mTask, float alpha);

GridData LevenbergMarkwardt(DirectSolverBase *directSolver, Task task, float alpha);


void MakeB(Matrix &mAnk, Matrix &mB);
void Makeb(Matrix &Ank, Matrix &Fnk, Matrix &b);