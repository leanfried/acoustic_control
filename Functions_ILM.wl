(* ::Package:: *)

(* ::Title:: *)
(*Functions for analysis of inverted light microscope images*)


(* ::Text:: *)
(*Last updated 06/01/16 by Leanne Friedrich*)
(*used by analysis_ILM.nb*)


(* ::Section:: *)
(*ILM Functions*)


(* ::Subsection:: *)
(*initialization*)


<<Functions_master.wl;
ILMscale = (1000./772); (*microns/pixel - this corresponds to the lowest magnification on our nikon microscope*)


(* ::Subsection:: *)
(*checking data*)


(* ::Text:: *)
(*getImages finds all images in a given directory and returns a list of them and a graphic*)
(*ILM interface collects the information from getImages and displays it in an interface*)


getImages[d_, bs_, MinInt_]:= Module[{images, ims},
Catch[
images = FileNames["ILM*.tif*", d];
If[Length[images]>0,
	ims = ImageResize[Import[#], 600]&/@images;
	ims = {Show[ImageAdjust@#, ImageSize->Large], 
					Show[getParticles[ImageAdjust@ImageCrop[#, ImageDimensions[#]*{0.6, 1}], bs, MinInt], ImageSize->Medium]}&/@ims;
	Return[Column[{d,
			Grid[Transpose@Flatten[{{StringTake[#, -7]&/@images}, Transpose@ims, {Button["Delete file", DeleteFile[#]]&/@images}},1]]}]]
	,
	Return[Column[{d, "No images"}]];
]]]


ILMInterface:=Manipulate[
	If[i<Length[directories]&& DirectoryQ[directories[[i]]], 
		getImages[directories[[i]], 0.85, 0.3]
	]
, {{i, 2, "Directory index"}, 1, Length[directories],1}

, ContinuousAction->False, 
ControlType->InputField];


(* ::Subsection:: *)
(*analysis*)


(* ::Text:: *)
(*getAllImageStats is called on directories and gets composite stats for the first 7 images in that folder, given critical segmentation intensity thresholds*)
(**)
(*this function needs getImageStats*)


(* ::Text:: *)
(*usually MaxInt = 0.85, MinInt = 0.3*)


getAllImageStats[d_, MaxInt_, MinInt_]:=Module[{images, x, retTable},
	Print[d]; (*print out the name of the folder*)
	images = FileNames["ILM*tif*", d];  (*get a list of all the images in the folder*)
	If[Length[images]>6, (*if there are at least 7 images, proceed*)
		images = images[[1;;7]]; (*take the first 7 images*)
		x = getImageStats[#, MaxInt, MinInt]&/@images; (*get stats for each image*)

		(*x will be a 7x2 matrix, where each row is for an image, the first column is stats about that image, and the second column is the intensity profile for that image*)

(*		Print[Row[x[[;;,2]]]];*)
(*		Export[StringJoin[d, StringJoin[{"\\ILM_", FileNameTake[d, -1],"_fits.tiff"}]], Row[x[[;;,2]]]];*)
(*Print[x[[;;,3]]];*)
		Export[StringJoin[d, StringJoin[{"\\ILM_", FileNameTake[d, -1],"_profiles.csv"}]], x[[;;,2]]]; (*this will save all of our intensity profiles in one file*)
		retTable = Append[
						Prepend[x[[;;,1]], 
								{"Mean center (\[Micro]m)", "Mean STDEV (\[Micro]m)", "STDEV center (\[Micro]m)", "STDEV STDEV (\[Micro]m)", "Composite center (\[Micro]m)", "Composite STDEV (\[Micro]m)"}
						], 
						Mean[x[[;;,1]]]
					]; (*retTable will be a compilation of all of the stats we measured in getImageStats. The last row is an average of the stats for all of the images*)
(*		Print[Grid[retTable]];*)
		Export[StringJoin[d, StringJoin[{"\\ILM_", FileNameTake[d, -1],"_stats.csv"}]] , retTable]; (*this exports all of our stats into one file*)
	];
]


(* ::Text:: *)
(*getImageStats gets stats for one image, given segmentation parameters*)
(**)
(*getImageStats needs getParticles and getSD (in masterFunctions)*)


getImageStats[file_, MaxInt_, MinInt_]:=Module[
			{cropHeight, im1, imFiltered, imageIntensities, imDims, numRows, numCols, ps,
				compositeStats, widths, flatImageIntensities, retTable, retImage},
Catch[
	im1 = Import[file]; (*get the image*)
	im1 = ImageAdjust@ImageCrop[im1, ImageDimensions[im1]*{0.6, 1}]; (*crop the image to the middle 60% and adjust the contrast*)
	imFiltered = getParticles[im1, MaxInt, MinInt]; (*segment the particles out of the image*)

	Export[StringJoin[StringSplit[file, "."][[1]], "_filt.bmp"], ColorNegate[AlphaChannel[imFiltered]]]; (*take the segmented binary image that we stored in the alpha channel of imFiltered and save it*)

	imageIntensities = ImageData[AlphaChannel[imFiltered]]; (*turn the binary image we stored in the alpha channel of imFiltered and turn it into a matrix of 0 and 1s*)
	flatImageIntensities = Total[imageIntensities];	(*take the sum of all of the rows, so you just have one intensity list from left-to-right across the image*)
	flatImageIntensities = flatImageIntensities/Max[flatImageIntensities]; (*normalize by the maximum intensity*)
	flatImageIntensities = Transpose[{Range[Length[flatImageIntensities]]*ILMscale, flatImageIntensities}]; (*turn this into a nx2 matrix of horizontal location within the image v. intensity*)

(*	imageIntensities = ImageData[imFiltered];
	imDims = Dimensions[imageIntensities][[1;;2]];
	numRows = imDims[[1]];
	numCols = imDims[[2]];
(*	cropHeight = Floor[0.2*numRows];
	widths = DeleteCases[findWidth[#, imageIntensities]&/@Range[cropHeight, numRows-cropHeight], {_,0,0}];*)
	flatImageIntensities = {#[[2]], imageIntensities[[#[[1]], #[[2]], 2]]}&/@Flatten[Table[{i,j}, {i, 1, numRows}, {j, 1,numCols}],1];
Print[flatImageIntensities];*)
	compositeStats = getSD[flatImageIntensities]; (*get the mean and standard deviation of this intensity profile*)
	retTable = Flatten[{(*Mean[widths[[;;,2;;]]], StandardDeviation[widths[[;;,2;;]]], *){0,0}, {0,0},compositeStats}];
(*ps = 2;
(*	retImage =Show[imFiltered, Graphics[{Red, AbsolutePointSize[ps], Point[#[[{2,1}]]]}]&/@widths, 
						Graphics[{Blue, AbsolutePointSize[ps], Point[{#[[2]]+#[[3]],#[[1]]}]}]&/@widths, 
						Graphics[{Blue, AbsolutePointSize[ps], Point[{#[[2]]-#[[3]],#[[1]]}]}]&/@widths]; *)*)
	Return[{retTable, (*Column[{(*Grid[Transpose[Prepend[{retTable}, {"Mean center", "Mean STDEV", "STDEV center", "STDEV STDEV", "Composite center", "Composite STDEV"}]]]*)
								Show[retImage, ImageSize->Large]}], *)flatImageIntensities}]
]]


(* ::Text:: *)
(*There's a bunch of commented out code in getParticles because I tried to do fancy shit, but it didn't work out, so I turned to using an aggressive segmentation technique instead.*)
(**)
(*getParticles finds dark particles on a light background. It essentially just uses MorphologicalBinarize, which finds all pixels above a high threshold (MaxInt) and all pixels attached to those pixels above a low threshold (MinInt)*)
(*It will store this segmented information in the alpha channel of the image, which is the transparency channel. As a result, ILMinterface can use this output to show you your original image, but set the background to transparent. getImageStats will just use the alpha channel as binary information (1 if this pixel contains a particle, 0 if this pixel does not contain a particle)*)


getParticles[im1_, MaxInt_, MinInt_]:=Module[{largestGroupInit, largestGroup, whiteScale, imFiltered, m1, badNumbers, f3, width, height, m2, intensities, m3},
(*	largestGroupInit = Times@@ImageDimensions[im1]; (*total # of pixels*)
	largestGroup = largestGroupInit; 
	whiteScale = 0.3; (*factor to filter out background*)
	imFiltered = ImageAdjust@RemoveBackground[im1 ,{White, whiteScale}]; (*take out white background*)
	m1 = MorphologicalComponents[imFiltered]; (*find connected components*)
	badNumbers = Select[Tally[Flatten[m1,1]], #[[2]]>largestGroupInit*MaxInt&][[2;;,1]]; (*get the large segments*)
	If[Length[badNumbers]>0, 
		f3 = ImageData[imFiltered][[;;,;;,1]];
		width = Dimensions[m1][[1]];
		height =  Dimensions[m1][[2]];
		m2 = Flatten[Table[{f3[[i,j]], m1[[i,j]]}, {i,width}, {j,height}],1]; (*table of intensities and alpha channel*)
		intensities = Transpose@{badNumbers, getIntensity[m2, #]&/@badNumbers}; (*get the average intensity of the large segments*)
		badNumbers = Select[intensities, #[[2]]>MinInt&][[;;,1]]; (*if they're too white, keep them in the list of bad numbers*)
		m3 =Replace[ m1,((#->0)&/@badNumbers),2]; (*take the bad segments out of the alpha channel*)
		m3 = m3/.x_/;x>0->1;
		imFiltered = SetAlphaChannel[imFiltered, Image[m3]]; (*insert the new alpha channel*)
	];*)
imFiltered = SetAlphaChannel[im1, MorphologicalBinarize[ColorNegate@im1, {MaxInt, MinInt}]];
	imFiltered

]


(* ::Subsection:: *)
(*Post checking*)


getILMProfile[dir_]:=Module[{x,y, xx},
x = ToExpression@Import[FileNames["*ILM*profiles.csv", dir][[1]]];
y =Import[FileNames["*ILM*stats.csv", dir][[1]]];
xx = Total[x[[;;,;;,2]]];
Row[{
ListLinePlot[Transpose[{x[[1,;;,1]], xx/Max[xx]}], PlotStyle->Black, FrameLabel->{"Width (micron)", "Intensity"}, ImageSize->300]
, "\t",
Grid[Transpose[Prepend[Transpose[y[[;;, {5,6}]]], {"", "Image 1", "Image 2", "Image 3", "Image 4", "Image 5", "Image 6", "Image 7", "MEAN"}]]]
}]]




(* ::Subsubsection::Closed:: *)
(*unused functions*)


(*findWidth[row_, imageIntensities_]:=Module[{row1},
	row1 = {#, imageIntensities[[row, #, 2]]}&/@Range[Length[imageIntensities[[row]]]];
	Flatten[{{Length[imageIntensities]-row}, getSD[row1]}]
];*)


(*getIntensity[m2_, number_]:=Module[{g1}, 
	g1 = Select[m2, #[[2]]==number&];
	Total[g1[[;;,1]]]/ Length[g1]
]*)
