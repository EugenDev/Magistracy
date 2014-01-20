#pragma once
#include <GridData.hpp>
#include "DirectSolver.hpp"
#include "SLESolver.hpp"
#include <Matrix.hpp>
#include <Task.hpp>
#include <MultilayerTask.hpp>
#include <vector>

GridData G_LinearisedSpeedDescent(Task task, float alpha);
GridData M_LinearisedSpeedDescent(Task task, float alpha);
void	 G_LinearisedSpeedDescentMultilayer();
GridData G_LinearisedMinimalError(Task task, float alpha);
GridData G_Newton(Task task);
vector<GridData> G_GradientMultilayer(MultilayerTask mTask, float alpha);