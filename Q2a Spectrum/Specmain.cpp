/*
Term Project 02a		====================================
Hoon Seok Kim		====================================
STD #: 200452816	====================================
*/
#include <complex>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <bitset>
using namespace std;

// Dimensions
#define N 256		// # of rows/columns
#define MaxR 182	// Max circle radius from origin to corner

// Buffer for FFT values
complex<double> DFT_ori[N][N];	// array for Original DFT result
complex<double> DFT[N][N];		// array for Centered DFT results

//Dimension of the vector, max radius of circle inside the plane
int M = N/2;

// Frequency Spectrum Array
complex<double> FH[N/2];			// Radius 1-128, 128+=beyond the ring
complex<double> FH2[MaxR];			// Radius 1-182, Enlarge circle to corner
complex<double> FHc[MaxR];			// cummulative spectrum up to origin to corner 
complex<double> sum = (0.0,0.0);	// variable for sum of magnitude of DFT.

int main()
{
//------------------------------Reading input file-----------------------
//-----------------------------------------------------------------------

	//-------------------------------------------------------------------
	//***********************************Select Part*********************
	//-------------------------------------------------------------------
	// Reading Original DFT values.
//	string inFileName = "DFTsquare.txt";		// activate for 1st image
	string inFileName = "DFTcar.txt";			// activate for 2nd image
	//-------------------------------------------------------------------
	//***********************************End of Select*******************
	//-------------------------------------------------------------------


	ifstream inFile;
	inFile.open(inFileName.c_str());
	if (inFile.is_open())
	{
		for (int row = 0; row < N; row++)
		{
			for (int col = 0; col < N; col++)
			{
				// DFT file to DFT_ori[][] array
				inFile >> DFT_ori[row][col];
			}
		}
		inFile.close();
	}
	else
	{
		cerr << "Can't find input file" << inFileName << endl;
	}
//---------------------------end of reading file-------------------------
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
// Centering F(0,0) & calculate sum of DFT
//-----------------------------------------------------------------------

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{	
			// Shifting array so DFT[0] can come to center
			DFT[(i + (N/2)) % N][(j + (N/2)) % N] = DFT_ori[i][j];
			// obtain sum of all magnitudes of DFT values 
			sum += (abs(DFT_ori[i][j]));
		}
	}
//-----------------------------------------------------------------------
// End of Center DFT
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
// Gettign frequency histogram by radius 1 to 182 = 128*sqrt(2)
// Max distance from origin to corner coordinate plane is apporx 182
//-----------------------------------------------------------------------

	// Increase Di from 1 to 182
	for (int i = 1; i <= MaxR; i++)
	{
		// Sweeping row dimension(x axis) value up to Di,
		// If x(r) > Di(i), it is out of circle
		// So we don't need to incread r more than i
		for (int r = 1; r <= i; r++)
		{
			// Sweeping column dimension(y axis) value up to Di,
			// If y(c) > Di(i), it is out of circle
			// So we don't need to incread c more than i
			for (int c = 1; c <= i; c++)
			{
				// Check if DFT values are in Di circle
				// sqrt( r^2 + c^2 ) < radius+1 from circle
				if (pow((r*r + c*c), 0.5) < i+1)
				{
					// x(r), y(c) cannot exceed M, or cannot be smaller than 0
					if (r >= 0 && r <= M && c >= 0 && c <= M)
					{
						// Apply x(r), y(c) to all quadrants 
						FHc[i] += (abs(DFT[r + 127][c + 127]))			// quadrant IV (x+,y-)
								+ (abs(DFT[r + 127][128 - c]))			// quadrant III(x-,y-)
								+ (abs(DFT[128 - r][c + 127]))			// quadrant I(x+,y+)
								+ (abs(DFT[128 - r][128 - c]));			// quadrant II(x-,y+)
					}
				}
			}
		}
	}
	
	// convert cummulative FHc to FH upto radius = M
	for (int z = 1; z < M; z++)
	{
		FH[z-1] = FHc[z] - FHc[z - 1];
	}
	// Add all remained DFT values to FH[last radius]
	FH[M - 1] = sum - FHc[M-1];
	// convert cummulative FHc to FH2 upto radius = origin to corner
	for (int z = 1; z <= 182; z++)
	{
		FH2[z - 1] = FHc[z] - FHc[z - 1];
	}
	// Print total sum of DFT values
	cout << "sum: " << sum << endl;
//-----------------------------------------------------------------------
// End of Frequency Spectrum
//-----------------------------------------------------------------------


//------------------------------saving output file-----------------------
//-----------------------------------------------------------------------

	//---------------------------------------------------------------
	//*******************************Select Part*********************
	//---------------------------------------------------------------
	// file out for FH array
//	ofstream myfile("FHsquare.txt");		// activate for 1st image
//	ofstream myfile2("FHsquaretoend.txt");	// activate for 1st image
//	ofstream myfile3("FHsquarecumm.txt");	// activate for 1st image

	ofstream myfile("FHcar.txt");			// activate for 2nd image
	ofstream myfile2("FHcartoend.txt");		// activate for 2nd image
	ofstream myfile3("FHcarcumm.txt");		// activate for 2nd image
	//---------------------------------------------------------------
	//*******************************End of Select*******************
	//---------------------------------------------------------------

	// file out for FH array
	if (myfile.is_open())
	{
		for (int item = 0; item < M; item++)
		{
			// take only real part since imaginary is 0
			myfile << FH[item].real() << endl;
		}
		myfile.close();
	}
	else cout << "Unable to open file";
	// file out for FH2 array
	if (myfile2.is_open())
	{
		for (int item = 0; item < MaxR-1; item++)
		{
			// take only real part since imaginary is 0
			myfile2 << FH2[item].real() << endl;
		}
		myfile2.close();
	}
	// file out for FHc array
	if (myfile3.is_open())
	{
		for (int item = 1; item < MaxR; item++)
		{
			// take only real part since imaginary is 0
			myfile3 << FHc[item].real() << endl;
		}
		myfile3.close();
	}
	return 0;
}