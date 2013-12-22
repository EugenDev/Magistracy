#pragma once
#include <iostream>
#include <fstream>
using namespace std;
#include <string>
#include "Matrix.hpp"

#define dllExport _declspec(dllexport)

enum dllExport SurferFileFormat
{
	DAT_FILE_FORMAT,
	SURFER7_FILE_FORMAT,
	SURFER6_FILE_FORMAT
};

class dllExport GridParameters 
{
public:
	float dX;
	float dY;
	float MinX;
	float MinY;
	int NX;
	int NY;
	bool operator == (const GridParameters other)
	{
		bool isEqual = true;
		isEqual &= dX == other.dX;
		isEqual &= dY == other.dY;
		isEqual &= NX == other.NX;
		isEqual &= NY == other.NY;
		isEqual &= MinX == other.MinX;
		isEqual &= MinY == other.MinY;
		return isEqual;
	}
	bool operator != (const GridParameters other)
	{
		bool isEqual = false;
		isEqual |= dX == other.dX;
		isEqual |= dY == other.dY;
		isEqual |= NX == other.NX;
		isEqual |= NY == other.NY;
		isEqual |= MinX == other.MinX;
		isEqual |= MinY == other.MinY;
		return !isEqual;
	}
	GridParameters()
	{
		NX = NY = 0;
		dX = dY = 0.0f;
		MinX = MinY = 0.0f;
	}
};

//TODO: Many understandable throws
class dllExport GridData
{
private:
	//Grid Parameters
	GridParameters gp;
	//Load grid from Surfer 7 file
	void ReadDSRBStream(ifstream &stream);
	//Load grid from .grd file
	void ReadGRDFile(string fileName);
	//Load grid from .grd file
	void ReadDATFile(string fileName);
	int GetLinesCount(ifstream* file);
	//Save grid to Surfer 7 file
	void SaveToDSRBFile(string filePath);
	//Save grid to Surfer 6 file
	void SaveToDSBBFile(string filePath);
	//Save grid to .dat file
	void SaveToDATFile(string fileName);
	//UselessValues
	double zMin, zMax, rotation, blank;
public:
	//GridData
	float* data;
	//Print grid info
	void ShowInfo();
	//Grid parameters accessor
	GridParameters GetGridParameters();
	//Constructors
	GridData(void);
	GridData(string inputFileName);
	GridData(const GridData &idol);
	GridData(GridParameters gridParameters);
	GridData(GridParameters gridParameters, Matrix matrix);
	//Get data as matrix
	Matrix AsMatrix();
	//Assignment operator overload
	GridData& operator = (GridData &other);
	//Fill matrix with grid data
	void FillMatrix(Matrix &m);
	//Destructor
	~GridData(void);
	//Saving Grid To File
	void SaveToFile(string outputFilePath, SurferFileFormat fileFormat);
};
