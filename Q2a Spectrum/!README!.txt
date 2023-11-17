--------------------------------------
Inputs from Q1: Original DFT values 
--------------------------------------
DFTcar.txt 
DFTsquare.txt
---------------

---------------------------------------------------
Outputs: Car folder, Square folder, barchart.xlsx
--------------------------------------------------------------------------------------------------------------------------
Car folder:
DFTcar.txt	==> FHcar.txt (FH[0] ~ FH[127], FH[127] contains magnitudes of DFT values on the corner)
		==> FHcar.png (Image for FH[] bar histogram)

		==> FHcartoend.txt (FH2[0] ~ FH2[181], Frequency spectrum radius extended to corner)
		==> FHcartoend.png (Image for FH2[] bar histogram)

		==> FHcarcumm.txt (FHc[1]~ FHc[182], cumulative bar histogram
		==> FHcarcumm.png (Image for FHc[] bar histogram)
Square folder:
DFTsquare.txt	==> FHsquare.txt (FH[0] ~ FH[127], FH[127] contains magnitudes of DFT values on the corner)
		==> FHsquare.png (Image for FH[] bar histogram)

		==> FHsquaretoend.txt (FH2[0] ~ FH2[181], Frequency spectrum radius extended to corner)
		==> FHsquaretoend.png (Image for FH2[] bar histogram)

		==> FHsquarecumm.txt (FHc[1]~ FHc[182], cumulative bar histogram
		==> FHsquarecumm.png (Image for FHc[] bar histogram)

DFTcar/DFTsquare.txt	==> barchart.xlsx ( All numerical values and histograms are included)
--------------------------------------------------------------------------------------------------------------------------