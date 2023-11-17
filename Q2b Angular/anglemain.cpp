/*
Term Project 02b		====================================
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

// Dimensions X,Y(Columns, rows)
#define N 256

// Buffer for FFT images
complex<double> DFT[N][N];		// Centered DFT result
complex<double> DFT_ori[N][N];	// Original DFT result

//Dimension of the vector: 360` / K
#define K 80

// Angular histogram
complex<double> AH[K];		// angular spectrum
complex<double> AHc[K];		// cumulative angular spectrum

// Define pi and d-theta by K
double pi = 2 * asin(1);
double dtheta = 2 * pi / K;

// variant for sum of each quadrant of DFT
complex<double> Q1, Q2, Q3, Q4;
complex<double> sum;

//==========================================================
// Sub Function ============================================
//==========================================================

// function for obtaining sum of quadrant I II III IV of DFT
// from DFT[][]
void sumQ()
{
	// initialize Q1-Q4 & sum
	Q1 = (0.0, 0.0);
	Q2 = (0.0, 0.0);
	Q3 = (0.0, 0.0);
	Q4 = (0.0, 0.0);
	sum = (0.0, 0.0);
	// Q1 = DFT Sum of quadrant I
	for (int r = 0; r < (N / 2); r++)
	{
		for (int c = (N/2); c < N; c++)
		{
			Q1 += abs(DFT[r][c]);
		}
	}
	// Q2 = DFT Sum of quadrant II
	for (int r = 0; r < (N / 2); r++)
	{
		for (int c = 0; c < (N/2); c++)
		{
			Q2 += abs(DFT[r][c]);
		}
	}
	// Q3 = DFT Sum of quadrant III
	for (int r = (N/2); r < N; r++)
	{
		for (int c = 0; c < (N / 2); c++)
		{
			Q3 += abs(DFT[r][c]);
		}
	}
	// Q4 = DFT Sum of quadrant IV
	for (int r = (N / 2); r < N; r++)
	{
		for (int c = (N/2); c < N; c++)
		{
			Q4 += abs(DFT[r][c]);
		}
	}
	// Total DFT values for the entire plane
	for (int r = 0; r < N; r++)
	{
		for (int c = 0; c < N; c++)
		{
			sum += abs(DFT[r][c]);
		}
	}
}
//==========================================================
// End of sub function definition===========================
//==========================================================


//=========================================================================
//===== Main Starts =======================================================
//=========================================================================
int main()
{
//-------------------------------------------------------------------------
//------------------------------File reading-------------------------------

	//-------------------------------------------------------------------
	//***************** Selection ***************************************
//	string inFileName = "DFTsquare.txt";		// activate only for 1st image
	string inFileName = "DFTcar.txt";			// activate only for 2nd image
	//***************** End of Selection ********************************
	//-------------------------------------------------------------------

	ifstream inFile;
	inFile.open(inFileName.c_str());

	if (inFile.is_open())
	{
		for (int row = 0; row < N; row++)
		{
			for (int col = 0; col < N; col++)
			{	// Load original DFT to DFT_ori[][]
				inFile >> DFT_ori[row][col];
			}
		}
		inFile.close();
	}
	else
	{
		cerr << "Can't find input file" << inFileName << endl;
	}
//------------------------------End reading file---------------------------
//-------------------------------------------------------------------------
	

//-------------------------------------------------------------------------
//------------------------------ Center DFT values ------------------------

	// Center DFT_ori[][] and save it to DFT[][]
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			DFT[(i + N/2) % N][(j + N/2) % N] = DFT_ori[i][j];
		}
	}
//------------------------------End centering DFT--------------------------
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//------------------------------ Angular spectum calculation --------------

	// Sum of DFT for each quadrant I ~ IV (Q1 ~ Q4) & total
	sumQ();

	// calculate i*d_theta ( 0 degree to 360 with K steps )
	for (int i = K-1; i >= 0; i--)
	{
		// When theta == 2 pi, add all values
		if ((dtheta)*(i + 1) >= (2 * pi))
		{
			// Sum of all DFT values ( all of Quadrant I II III IV )
			AHc[i] = Q1 + Q2 + Q3 + Q4;
		}

		// When 3/2 pi < theta < 2 pi
		else if ( ((dtheta)*(i + 1)) < (2*pi) && ((dtheta)*(i + 1)) > (1.5 * pi))
		{
			// Sum of quadrant (I II III) of DFT 
			AHc[i] = Q1 + Q2 + Q3;

			// Sum of quadrant (IV) of DFT depends on angle
			for (int col = (N/2); col < N; col++)
			{
				// x*tan(theta)=> minus
				for (int row = (int)(127 - (col - 127)*tan((dtheta)*(i + 1))); row < N; row++)
				{
					// ignore if x*tan(theta) values are out of row range(128~255) 
					if (row < N && row >= (N/2))
					{
						AHc[i] += abs(DFT[row][col]);
					}
				}
			}
		}

		// When theta is 3/2 pi, sum of quadrnat (I,II,III) of DFT
		else if ( ((dtheta)*(i + 1)) == (1.5*pi) )
		{
			// Sum of quadrant (I II III) of DFT 
			AHc[i] = Q1 + Q2 + Q3;
		}

		// When pi < theta < 3/2 pi
		else if ( (((dtheta)*(i + 1)) < (1.5*pi)) && ((dtheta)*(i + 1)) > pi )
		{
			// Sum of quadrant (I && II) of DFT 
			AHc[i] = Q1 + Q2;
			
			// Sum of quadrant (III) of DFt depends on angle
			for (int col = 0; col < (N/2); col++)
			{
				for (int row = (int)(127 + (128 - col)*tan((dtheta)*(i + 1))); row >= (N/2); row--)
				{
					if (row < N && row >= (N/2))
					{
						AHc[i] += abs(DFT[row][col]);
					}
				}
			}
		}

		// When pi == theta
		else if (((dtheta)*(i + 1)) == pi)
		{
			// Sum of quadrant (I && II) of DFT 
			AHc[i] = Q1 + Q2;
		}

		// When 0.5*pi < theta < pi
		else if ( (((dtheta)*(i + 1)) > (0.5*pi) ) && ((dtheta)*(i + 1)) < pi)
		{
			// Sum of quadrant (I) of DFT 
			AHc[i] = Q1;

			// Sum of quadrant (II) of DFt depends on angle
			for (int col = 0; col < (N/2); col++)
			{
				for (int row = (int)(128 + (128 - col)*tan((dtheta)*(i + 1))); row >= 0; row--)
				{
					if (row < (N/2) && row >= 0)
					{
						AHc[i] += abs(DFT[row][col]);
					}
				}
			}
		}

		// When 0.5 pi == theta
		else if (((dtheta)*(i + 1)) == 0.5*pi)
		{
			// Sum of quadrant (I) of DFT 
			AHc[i] = Q1;
		}

		// When 0 < theta < 0.5 * pi
		else if ((dtheta)*(i + 1) < (0.5 * pi) && ((dtheta)*(i + 1)) > 0)
		{
			// Sum of quadrant (I) of DFt depends on angle
			for (int col = (N/2); col < N; col++)
			{
				for (int row = (int)(128 - (col - 127)*tan((dtheta)*(i + 1))); row < (N/2); row++)
				{
					if (row < (N / 2) && row >= 0)
					{
						AHc[i] += abs(DFT[row][col]);
					}
				}
			}
		}

		// When theta == 0
		else
		{
			AHc[i] = AHc[i];
		}
	}

	// obtain angular spectrum AH[0 ~ K-1] from cumulative AH[0~K-1] 		
	for (int i = 0; i < K; i++)
	{
		if (i > 0)
		{
			AH[i] = AHc[i] - AHc[i - 1];
		}
		else
		{	// AH[0] = AHc[0]
			AH[i] = AHc[i];
		}
	}
	// Compare AHc[K-1] and sum if those are the same.
	cout << "AHC last: " << AHc[K - 1] << " sum: " << sum;

//------------------------------End centering DFT--------------------------
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//------------------------ Save outputs -----------------------------------

	//-------------------------------------------------------------------
	//***************** Selection ***************************************
	// file out for AH array
//	ofstream myfile("AHsquare.txt");		// activate only for 1st image
	ofstream myfile("AHcar.txt");			// activate only for 2nd image
	//***************** End of Selection ********************************
	//-------------------------------------------------------------------
	if (myfile.is_open())
	{
		for (int item = 0; item < K; item++)
		{	// save real part 'cause imaginary part is 0
			myfile << AH[item].real() << endl;
		}
		myfile.close();
	}
	else cout << "Unable to open file";
	return 0;
}