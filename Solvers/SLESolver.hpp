#pragma once
#include <Matrix.hpp>
#include <Time.h>

void SolveSLE(const Matrix &A, Matrix &Z, const Matrix &b, float initialZ);
