#pragma once
#include <cublas.h>
#include "CudaDirectSolver.hpp"
#include <Matrix.hpp>
#include <GridData.hpp>
#include <Task.hpp>
#include <time.h>
#include <vector>

GridData LightLevenbergMarkvardt(Task t, float alpha);
GridData LightLinearisedMinimalError(Task t, float alpha);
GridData LightLinearisedSpeedDescent(Task t, float alpha);
vector<GridData> LightMultilayerLevenbergMarkvardt(MultilayerTask mTask, float alpha);
vector<GridData> LightMultilayerLinearisedMinimalError(MultilayerTask mTask, float alpha);
vector<GridData> LightMultilayerLinearisedSpeedDescent(MultilayerTask mTask, float alpha);