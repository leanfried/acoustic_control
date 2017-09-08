(* ::Package:: *)

(* ::Section:: *)
(*Guppy functions*)


<<Functions_master.wl;
video = FileNames["guppy*CSph*.avi", BASEDIR, Infinity]; 
	(*all guppy videos*)
GS = 1.15; (*this is the scale for guppy videos, in px/micron. You should measure this scale using imageJ and change this if necessary*)


(* ::Subsection:: *)
(*analyze guppy videos*)


(* ::Subsubsection:: *)
(*takes an intensity line profile from point {x1,y1} to {x2, y2} on an image using npts points and background fitting parameter bParam*)


(* ::Text:: *)
(*used by: CalibrateProfile*)


(* ::Text:: *)
(*ImageProfile gets an intensity profile between two points on an image, removes the background, and rescales all values so 1 is the max value*)
(*ImageProfileRAW gets an intensity profile between two points on an image, and rescales all values to a range from 0 to 1*)
(**)
(*the current version of the analysis function uses ImageProfileRAW*)


ImageProfileRAW[img_,{x1_,y1_},{x2_,y2_},npts_,  crop_]:=Module[{profile, pback, profile2, pmax},
(*img = entire frame
	{x1, y1} = left point, in pixels
	{x2, y2} = right point, in pixels
	npts = the number of points to collect
	bParam = the parameter for estimating the background NOT USED
	crop = percentage from left to crop*)
Catch[
	profile = Table[{N@((i-1)/(npts-1)) Sqrt[(x1-x2)^2+(y1-y2)^2],
						ImageValue[img,{x1+(x2-x1) (i-1)/(npts-1),y1+(y2-y1) (i-1)/(npts-1)}]}
				,{i,1,npts}];
		(*table column 1 = position in pixels, column 2 = intensity*)
	profile = profile[[Max[1, Floor[crop*npts]];;]];
	profile2 = profile;
	profile2[[;;,2]] = Rescale[profile2[[;;,2]]];
(*	(*diagnostic*)Print[ListLinePlot[{profile, pback, profile2} , PlotRange\[Rule]{{0, 350*GS}, {0,1}}]];*)
	Return[profile2];
]];


(* ::Subsubsection:: *)
(*old imageProfile*)


ImageProfile[img_,{x1_,y1_},{x2_,y2_},npts_, bParam_, crop_]:=Module[{profile, pback, profile2, pmax},
(*img = entire frame
	{x1, y1} = left point, in pixels
	{x2, y2} = right point, in pixels
	npts = the number of points to collect
	bParam = the parameter for estimating the background
	crop = percentage from left to crop*)
Catch[
	profile = Table[{N@((i-1)/(npts-1)) Sqrt[(x1-x2)^2+(y1-y2)^2],
						ImageValue[img,{x1+(x2-x1) (i-1)/(npts-1),y1+(y2-y1) (i-1)/(npts-1)}]}
				,{i,1,npts}];
		(*table column 1 = position in pixels, column 2 = intensity*)
	profile = profile[[Max[1, Floor[crop*npts]];;]];
	pback = profile; (*pback = background of intensity profile*)
	pback[[;;,2]] = EstimatedBackground[profile[[;;, 2]],bParam]; (*find background*)
	profile2 = profile; (*profile2 is the profile to return*)
	profile2[[;;,2]] = profile[[;;,2]] - pback[[;;,2]];(*subtract background*)
	profile2[[;;,2]] = profile2[[;;,2]]/Max[profile2[[;;,2]]]; (*normalize*)

(*	(*diagnostic*)Print[ListLinePlot[{profile, pback, profile2} , PlotRange\[Rule]{{0, 350*GS}, {0,1}}]];*)

	Return[profile2];
]];


(* ::Subsubsection:: *)
(*finds left edge of channel*)


(* ::Text:: *)
(*used by: CalibrateVideo*)


(* ::Text:: *)
(*straightedges->0.1*)


getLine[video_, diagnosticMode_, start_, slope1_, slope2_, seParam_]:=Module[{frame1, smallIm, line, p1, vector, unit, p2, newLeft, lw, id3,
													line2, p12, p22, scale, numFrames, smallIm0,  smallIm2, pixelRegion, lineWidth, id, id2},
(*video = video,
	length is fractional distance from the top of the channel across which to take a cross-section,
	diagnosticMode is boolean - whether or not to print pictures and other diagnostics,
	start = first frame*)
Catch[
(*find channel*)

numFrames = Import[video, "FrameCount"]; (*number of frames in video*)
frame1 = ColorConvert[Import[video, {"Frames", Floor[numFrames*start]}], "Grayscale"]; (*first frame we will measure*)
scale = 4; (*use a larger number here to compress the image more - if this number is too big, you might not find the line*)
line = {}; (*line is the list of possible lines to use*)
While[Length[line]<1, (*while we don't have any viable lines, keep running this*)
	smallIm0 = ImageResize[frame1, 2588/scale]; (*compress the image image*)
	smallIm = HistogramTransform@smallIm0; (*normalize intensities across the frame*)
	pixelRegion = 7; (*pixel distance over which to detect edges*)
	While[Length[line]<1 && pixelRegion>0, (*if no lines are found, decrease the contrast radius*)
		smallIm2 = ImageMultiply[EdgeDetect[smallIm0,pixelRegion, Method->{"ShenCastan", "StraightEdges"->seParam}],smallIm]; (*find edges and multiply by the intensity of the image at each point*)
		id = ImageData[smallIm2]; (*Mathematica stores images using a special data structure - this converts the image into a matrix of numbers*)
		id2 = Table[If[id[[i,j]]>0,{i,j}, {0,0}], {i, Dimensions[id][[1]]}, {j, Dimensions[id][[2]]}]; (*get a list of indices where EdgeDetect found an edge*)
		id3 = DeleteCases[#, {0,0}]&/@id2; (*only take the indices where EdgeDetect found an edge*)
		lw = MinMax[id3[[;;,;;,2]]];  (*find the leftmost and rightmost index where EdgeDetect found an edge*)
		lineWidth = (#2-#1)&@@lw; (*lineWidth is the distance between the leftmost and rightmost index*)
		If[diagnosticMode, Print[1.*lineWidth/(Dimensions[id][[2]])]]; 
		If[0.14<lineWidth/(Dimensions[id][[2]])<0.17, 
						(*if the line width is smaller than an arbitary width - I chose 0.14 and 0.17 because they worked, but you can change this - 
							then treat the found edges as the channel region and find the left edge from that region. otherwise, use line detection (Hough transform)*)
			line = {{{Min[DeleteCases[0, Flatten[id3[[-5;;-1, lw[[1]];;lw[[2]], 2]]]]], 0}, {Min[DeleteCases[0, Flatten[id3[[1;;5, lw[[1]];;lw[[2]], 2]]]]], Dimensions[id][[2]]}}}*scale; 
				(*draw a line between the leftmost point in the top 5 lines and the leftmost point in the bottom 5 lines*)
			,
			line = (ImageLines[smallIm2,0.125(*, MaxFeatures\[Rule]100*)]*scale); (*find lines*)
			line = Select[line, 200<#[[1,1]]&]; (*only take lines to the right of 200 px*)
			line = Select[line, slope1<#[[1,1]] - #[[2,1]]<slope2&]; (*only take lines that are nearly vertical*)
			pixelRegion = pixelRegion-1 (*if this doesn't work, detect edges over a smaller region*)
			];
	];
	scale = scale+1 (*if no lines are found, scale the image to a smaller size*)
];

If[diagnosticMode, 
(*diagnostic*) Print[smallIm2](* ; Print[smallIm]*)
(*diagnostic*) Print[Show[frame1, Graphics[{Thick, Orange, Line@line}], ImageSize->Large]];
];

line = MinimalBy[line, First][[1]]; (*take line that is farthest to the left*)

vector = {line[[1,2]] - line[[2,2]], -line[[1,1]] + line[[2,1]]}; (*direction of cross-section*)
unit = Normalize[vector];

Return[{line,unit}];
]]


(* ::Subsubsection:: *)
(*get positions from left edge and height*)


(* ::Text:: *)
(*used by: CalibrateProfile*)


getp[line_,unit_,  length_]:=Module[{ p1, p2},Catch[
	p1 = line[[1]] + length*(line[[2]] - line[[1]]); (*left point*)
	p2 = p1 + unit*(350)*GS; (*right point*)
	Return[{p1, p2}];	
]]


(* ::Subsubsection:: *)
(*trim a frame to include just area of interest*)


(* ::Text:: *)
(*used by: calibrateVideo*)


trimVideo[frame1_, line_, diagnosticMode_]:=Module[{line2, p12, p22, newLeft, frame},
Catch[

frame = Show[frame1, Graphics[{Red, Thick, Line@line}]];
(*If[diagnosticMode, Print[Show[frame, ImageSize\[Rule]Large]]];*)
frame = ImageTrim[frame, {line[[1]], 0.4*(line[[2]]-line[[1]])+line[[1]],line[[1]]+{400,0}, 0.4*(line[[2]]-line[[1]])+line[[1]]+{400,0} }]; 
	(*trim to width of channel*)
frame = ImageTake[frame, ImageDimensions[frame][[2]]*Norm[0.4*(line[[2]]-line[[1]])]];  (*trim image*)

Return[{frame}]]
]


(* ::Subsubsection::Closed:: *)
(*calibrateProfile old junk*)


(*	profile = ConstantArray[0, {npts, 2}];
	numFrames = Length[Import[video, {"Frames"}]]; (*total length of video, in frames*)
	lowBound =  Floor[numFrames*start]; (*index of first frame to take*)
	upBound =  Floor[numFrames*finish]; (*index of last frame to take*)
	numFramesToGet = 3; (*average 5 frames*)
	stepSize = (upBound - lowBound)/numFramesToGet;
	Do[
		index = Floor[lowBound + (j-1)*stepSize];
		frame1 =  ColorConvert[Import[video, {"Frames", index}], "Grayscale"]; (*import frame*)
		Do[
			{p1, p2} = getp[line, unit, length]; (*find cross-sectional region*)
			p = ImageProfile[frame1, p1, p2,npts, bParam]; (*take profile of frame*)
			profile[[;;,2]] = profile[[;;,2]] + p[[;;,2]]; (*add to sum*)
		, {length, 0.1, 0.4, 0.05}]
	,{j, numFramesToGet}];
	profile [[;;,2]] = profile[[;;,2]]/Max[profile[[;;,2]]]; (*normalize*)	
	profile[[;;,1]] = p[[;;,1]]; (*put distances into column 1 of profile*)

	num = 20; 
		(*distance to take 1st derivative - increase for less noise, decrease for more precision*)
	maxPos = Ordering[profile[[1;;100,2]], -1][[1]]; 
		(*far left peak corresponding to reflection off channel wall*)
	profile = Transpose[Append[Transpose[profile], EstimatedBackground[profile[[;;, 2]], num]]]; 
		(*find the close background*)
	profile = Transpose[Append[Transpose[profile], Join[profile[[num;;-1,3]] -profile[[1;;-num, 3]], ConstantArray[0, num-1]]]]; 
		(*1st derivative*)
	profile= Transpose[Append[Transpose[profile], Join[profile[[num;;-1, 4]] - profile[[1;;-num, 4]],ConstantArray[0, num-1]]]];  
		(*2nd derivative*)
	newMin = FirstPosition[profile[[maxPos;;]], {x_ ,noise_ ,in_, d1_, d2_}/;d1>0 && in<0.5 && d2>0][[1]] + maxPos;
		(*first position after peak where intensity<30%, slope>0, curvature>0*)
	newMin = Min[newMin, 250];	(*crop at 25% if the minimum is too far to the right*)*)
(*	profile2 = profile[[newMin;;, {1,2}]]; (*crop*)
	profile2[[;;,2]] = profile2[[;;,2]]/Max[profile2[[;;,2]]]; (*normalize*)*)


(* ::Subsubsection:: *)
(*originally calibrated a line profile to remove left reflection*)


(* ::Text:: *)
(*needs: ImageProfile, getp*)


CalibrateProfile[video_,line_,unit_, npts_,  pic_, start_]:=Module[
			{frame1, profile, pmax, maxPos, earlyMax, newMin, numFrames, p, lateMax, stage1, 
			stage2, num, lowBound, upBound, stepSize, index, numFramesToGet, profile2, p1, p2, points},
Catch[
	numFrames = Length[Import[video, {"Frames"}]]; (*total length of video, in frames*)	
	frame1 =  ColorConvert[Import[video, {"Frames", Floor[start*numFrames]}], "Grayscale"]; (*import frame*)
	{p1, p2} = getp[line, unit, 0.2]; (*find cross-sectional region*)
	profile2 = ImageProfileRAW[frame1, p1, p2, npts,  0.25]; (*crop fffrom the left*)
	If[pic,
(*		stage2 = ListLinePlot[{profile2[[;;, 1;;2]](*, profile[[;;, {1,3}]], profile[[;;, {1,4}]], profile[[;;, {1,5}]]*)} , 
							PlotRange->{{0, 350*GS}, {0,1}}, PlotStyle->{Red, Thickness[0.02]}, Frame->None, Axes->None, Ticks->None];
		*)Return[{profile2, frame1}];
	,
		Return[profile2]
	];
]];


(* ::Subsubsection:: *)
(*sets up vision parameters to analyze guppy videos*)


(* ::Text:: *)
(*needs: getLine, CalibrateProfile, getSD, trimVideo*)


calibrateVideo[video_, start_, diagnosticMode_, slope1_, slope2_, seParam_]:=Module[
							{line, unit, profile, frame1},
Catch[
	If[!StringQ[video], Return["Video is invalid"]];
	{line, unit} = getLine[video, diagnosticMode, start, slope1, slope2, seParam]; (*find left edge*)
	{profile, frame1} = CalibrateProfile[video,line,unit, 1000, True, start]; (*calibrate profile region*)
	(*frame1 = ImageTrim[frame1, Flatten[{line, {p2}},1]];*)
(*	{m,v} = getSD[profile]; (*get stats on profile*)*)
	(*frame1 = ImageCompose[frame1,  Graphics[{Thick, Orange, Line@line2, Blue, Line[{p12, p22}], Red, PointSize[0.05](*, Point[{p12 + unit*m}]*)}]];
	*)
(*	{p1, p2} = getp[line, unit, 0.2];*)
	{frame1} = trimVideo[frame1, line, diagnosticMode]; (*trim frame for exporting sample*)
(*	If[diagnosticMode, Print[frame1]];*)
(*	frame1 = ImageCompose[frame1, profilePlot, Scaled[{0.5, 0.5}]];  (*put profile onto image*)*)
	If[diagnosticMode, Print[Show[frame1, ImageSize->Medium]]];
	Return[{profile, line, frame1, unit}];
]]


(* ::Subsubsection:: *)
(*analyzes and exports profile diagnostics from guppy videos*)


(* ::Text:: *)
(*usually [video, 0.25, 0.75, -20, 120, 0.1]*)


analyzeAndExport[video_, start_, finish_, slope1_, slope2_, seParam_]:=Module[
				{numFrames, lowBound, upBound, results,profile, profTotal, newSteps, unit, length,
					p1adj, p1, p2, index, frame1, frame2 , stepSize, numFramesToGet, line, h1, h2, hstep, numhsteps,
				resultsGrid, resultsGrid2, profiles},
Catch[
	If[!StringQ[video], Return["Video is invalid"]];
	Print[video]; (*print name of video*)
	{profile, line, frame1, unit} = calibrateVideo[video, start, False, slope1, slope2, seParam]; (*calibrate profile location*)
	p1adj = unit*profile[[1,1]]; (*p1adj is the adjustment made to the left bound*)
	numFrames = Length[Import[video, {"Frames"}]]; (*number of frames*)
	lowBound =  Floor[numFrames*start]; (*first frame*)
	upBound =  Floor[numFrames*finish]; (*last frame*)
	numFramesToGet = 8;
	stepSize = (upBound - lowBound)/(numFramesToGet-1); 

	h1 = 0.05; (*first height*)
	h2 = 0.4; (*last height*)
	hstep = 0.01; (*height step size*)
	numhsteps = Floor[(h2-h1)/hstep+1];

	results = ConstantArray[0, {numFramesToGet*numhsteps, 4}]; (*mean and stdev results*)
	profTotal = ConstantArray[0, {Length[profile],2}]; (*composite profile*)
	profTotal[[;;,1]] = profile[[;;,1]]; (*copy over distances*)
	newSteps = Length[profile]; (*number of steps to take in intensity plot*)
	PrintTemporary[Row[{Dynamic[k], "/", numFramesToGet}]]; (*monitor*)
	
	profiles = ConstantArray[0, {numFramesToGet*numhsteps, Length[profile]}];

	Do[
		index = Floor[lowBound + (k-1)*stepSize]; (*frame number*)
		frame2= ColorConvert[Import[video, {"Frames", index}], "Grayscale"]; (*frame*)

		Do[
			length = h1 + (l-1)*hstep;
			{p1, p2} = getp[line,unit, length]; (*find cross-sectional region*)
			profile[[;;,2]] = ImageProfileRAW[frame2, p1 + p1adj, p2,newSteps, 0][[;;,2]]; (*take profile of frame/height*)
			profTotal[[;;,2]] = profTotal[[;;,2]] + profile[[;;,2]]; (*add to sum*)
			profiles[[(k-1)*numhsteps+l, ;;]]=profile[[;;,2]]; (*store this profile in our big matrix 'o' profiles*)
			results[[(k-1)*numhsteps+l, ;;]] = Flatten[{{index, length}, getSD[profile]}]; (*store mean and SD results*)
		, {l, numhsteps}]
	,{k,numFramesToGet}];
	profTotal[[;;,2]] = profTotal[[;;,2]]/ Max[profTotal[[;;,2]]]; (*normalize composite intensity profile*)
(*Print[ListLinePlot[Transpose[{profile[[;;,1]], #}]&/@profiles, PlotRange\[Rule]{{0, 350*GS}, {0,1}}]];*)

	resultsGrid = {{"Mean center", "Mean STDEV"}, Mean[results[[;;, {3,4}]]], 
					{"STDEV center", "STDEV STDEV"}, StandardDeviation[results[[;;, {3,4}]]], 
					{"Composite Center", "Composite STDEV"}, getSD[profTotal]};
	resultsGrid2 = {{"Mean center", "Mean STDEV", "STDEV center", "STDEV STDEV","Composite Center", "Composite STDEV"}, 
					Flatten[{Mean[results[[;;, {3,4}]]], StandardDeviation[results[[;;, {3,4}]]],  getSD[profTotal]}, 1]};

	Export[StringJoin[StringSplit[video, ".avi"][[1]], "_XSlocation_RAW.csv"], Flatten[{line, {unit}}, 1]]; (* - the location of the channel*)
	Export[StringJoin[StringSplit[video, ".avi"][[1]], "_profile_RAW.tiff"], frame1]; (* - the composite profile superimposed on a frame cropped to the channel width*)
	Export[StringJoin[StringSplit[video, ".avi"][[1]], "_data_RAW.csv"], resultsGrid2]; (* - stats about the video*)
	Export[StringJoin[StringSplit[video, ".avi"][[1]], "_meanProfile_RAW.csv"], Prepend[profTotal, {"Position", "Intensity"}]]; (* - the composite intensity profile*)
	Export[StringJoin[StringSplit[video, ".avi"][[1]], "_means_RAW.csv"], Prepend[results, {"Frame index", "Height", "Mean", "STDEV"}]]; (* - stats for each line profile*)
	Export[StringJoin[StringSplit[video, ".avi"][[1]], "_profiles_RAW.csv"], profiles]; (* - all captured line profiles*)
	Return[Row[{Show[frame1, ImageSize->150], 
			Grid[	resultsGrid], 
			ListLinePlot[profTotal , PlotRange->{{0, 350*GS}, {0,1}}, FrameLabel->{"Width (px)", "Intensity"}, PlotStyle->Thick, ImageSize->Medium]}, "\t"]];
]
];


(* ::Subsubsection:: *)
(*gets all files related to a specific guppy video*)


getGuppy[dir_]:=Module[{mp, stats, prof},
	If[Length[FileNames["guppy*", dir]]<1,
			Print["No guppy file"];
		,
		mp = Import[FileNames["*meanProfile_RAW.csv", dir][[1]]];
		stats = Import[FileNames["*data_RAW.csv", dir][[1]]];
		prof = Import[FileNames["*profile_RAW.tiff", dir][[1]]];
		Column[{dir, Row[{Show[prof, ImageSize->150], 
			"\t",
			ListLinePlot[mp, PlotRange->{{0, 350}, {0,1}},PlotStyle->Black,ImageSize->300], 
				Grid[Transpose@stats]}]}]]


]


(* ::Subsubsection:: *)
(*takes all files from the same run and puts them in a folder together*)


testDirectory[num_]:=Module[{x, newDirectory, filesToMove},
x = FileNameSplit[video[[num]]];
x[[-1]] = StringSplit[video[[num]], {"guppy_", ".avi"}][[2]];
newDirectory =  FileNameJoin[x];
Print[newDirectory];
If[StringTake[x[[-2]], 1]=="V", 
If[!DirectoryQ[newDirectory], CreateDirectory[newDirectory]];
If[DirectoryQ[newDirectory], 
filesToMove = FileNames[StringJoin["*", StringSplit[video[[num]], {"guppy_", ".avi"}][[2]], "*"], "*", Infinity];
RenameFile[#, FileNameJoin[Append[x, FileNameTake[#]]]]&/@filesToMove;
];
]
];


(* ::Input::Plain:: *)
(**)


(* ::Subsection:: *)
(*post processing*)


(* ::Subsubsection::Closed:: *)
(*old post processing stuff*)


getRow[file_, pressures_]:=Module[{x, y, r, z},
x = parseFileName[file];
y = Flatten[Import[file, "csv"]];
y = y/GS; (*this converts from pixel to micron*)
z = Select[pressures, #[[1]]==ToExpression@x[[2]] && #[[2]]==ToExpression@x[[3]] && #[[3]]==ToExpression@x[[7]]&][[1, 4;;]];
r = Join[x, y, z]
]


getRowHA[file_, pressures_]:=Module[{x, y, r, z},
x = parseFileNameHA[file];
y = Import[file, "csv"][[2]];
y = y/GS; (*this converts from pixel to micron*)
z = Select[pressures, #[[1]]==ToExpression@x[[2]] && #[[2]]==ToExpression@x[[3]] && #[[3]]==ToExpression@x[[7]]&][[1, 4;;]];
r = Join[x, y, z]
]


parseFileName[file_]:=Module[{x, lf, sio, a , CSph, V, P, S, t},
x = StringSplit[FileNameTake[file], "_"];
lf = StringTake[x[[2]], 3;;];
sio = StringTake[x[[3]], 4;;];
a = StringTake[x[[4]], 2;;];
CSph = StringTake[x[[5]], 5;;];
V = StringTake[x[[6]], 2;;];
P = StringTake[x[[7]], 2;;];
S = StringTake[x[[8]], 2;;];
t = StringTake[StringJoin@Riffle[x[[9;;-2]], "_"], 2;;];

{lf, sio, a, CSph, V, P, S, t}
]


parseFileNameHA[file_]:=Module[{x, lf, sio, a , CSph, V, P, S, t},
x = StringSplit[FileNameSplit[file][[-1]], "_"];
lf = StringTake[x[[2]], 3;;];
sio = StringTake[x[[3]], 4;;];
a = StringTake[x[[4]], 2;;];
CSph = StringTake[x[[5]], 5;;];
V = StringTake[x[[6]], 2;;];
P = StringTake[x[[7]], 2;;];
S = StringTake[x[[8]], 2;;];
t = StringTake[StringJoin[Riffle[x[[9;;-3]], "_"]], 2;;];

{lf, sio, a, CSph, V, P, S, t}
]


(* ::Subsubsection:: *)
(*get file info out of a name of a guppy video*)


parseGuppyName[file_]:=Module[{x, lf, sio, a , CSph, V, P, S, t},
x = StringSplit[FileNameSplit[file][[-1]], "_"];
lf = StringTake[x[[2]], 3;;];
sio = StringTake[x[[3]], 4;;];
a = StringTake[x[[4]], 2;;];
CSph = StringTake[x[[5]], 5;;];
V = StringTake[x[[6]], 2;;];
P = StringTake[x[[7]], 2;;];
S = StringTake[x[[8]], 2;;];
t = StringSplit[StringTake[StringJoin[Riffle[x[[9;;-1]], "_"]], 2;;], ".avi"][[1]];

{lf, sio, a, CSph, V, P, S, t}
]


(* ::Subsubsection::Closed:: *)
(*old stat compilation*)


reCompileStats:=Module[{pressures},
	pressures = Import["C:\\Users\\Leanne\\Box Sync\\Rheology\\Pressures.txt", "csv"];
files = FileNames["*CSph120*data.txt*", "*", Infinity]; 
	allStats = Prepend[getRow[#, pressures]&/@files, Join[header, pressures[[1,4;;]]] ]; (*this was for initial compilation*)
	Export["Processed guppy videos HA.txt", allStats, "csv"];
	allStats = Import["Processed guppy videos HA.txt", "csv"];
]


reCompileStatsHA:=Module[{pressures},
	pressures = Import["C:\\Users\\Leanne\\Box Sync\\Rheology\\Pressures.txt", "csv"];
files = FileNames["*CSph120*data_HA.txt*", "*", Infinity]; 
	allStats = Prepend[getRowHA[#, pressures]&/@files, Join[header, pressures[[1,4;;]]] ]; (*this was for initial compilation*)
	Export["Processed guppy videos HA.txt", allStats, "csv"];
	allStats = Import["Processed guppy videos HA.txt", "csv"];
]


(* ::Subsubsection::Closed:: *)
(*old stat gathering*)


(* ::Text:: *)
(*get all stats for a file*)


getStats[file_]:=Module[{files, pressure, folder, guppyDiagnostics, meanProfile, calibrationImage,means},
folder =  FileNameTake[file,{1,-2}];
(*files = FileNames["*",folder];*)
pressure = Import[ FileNames["*data*csv*",folder][[1]], "csv"];
guppyDiagnostics = Flatten[Import[ FileNames["*guppy*data*",folder][[1]], "csv"]]/GS;
meanProfile = Import[ FileNames["*guppy*meanProfile*",folder][[1]], "csv"];
meanProfile[[;;, 1]] = meanProfile[[;;,1]]/GS;
means = Import[ FileNames["*guppy*means*",folder][[1]], "csv"];
means[[2;;]] = means[[2;;]]/GS;
calibrationImage = Import[ FileNames["*guppy*profile.tiff",folder][[1]]];
{Grid@Transpose@{header[[{1, 6,8}]], parseFileName[file][[{1, 6,8}]]},

ListLinePlot[pressure[[;;,{1,2}]], ImageSize->Medium, Frame->True, FrameLabel->{"Time (s)", "Pressure (mBar)"}, PlotRange->{All, {0, 7000}}],
Row[{Grid[means, Frame->All], 
Grid[Transpose@{{"Mean center (\[Micro]m)", "Mean STDEV (\[Micro]m)", "STDEV center (\[Micro]m)", "STDEV STDEV (\[Micro]m)", "Composite center (\[Micro]m)", "Composite STDEV (\[Micro]m)"}, guppyDiagnostics}]
}],ListLinePlot[meanProfile, ImageSize->325, FrameLabel->{"Distance (\[Micro]m)", "Frequency"},PlotLabel->"Mean profile",  Frame->True, PlotRange->{{0, 350}, {0, 1}}],
Show[calibrationImage, ImageSize->280]
}
]


(* ::Text:: *)
(*get all stats for height scanned files*)


getStatsHA[file_]:=Module[{files, pressure, folder, guppyDiagnostics, meanProfile, calibrationImage,means},
folder =  FileNameTake[file,{1,-2}];
(*files = FileNames["*",folder];*)
pressure = Import[ FileNames["*data*csv*",folder][[1]], "csv"];
guppyDiagnostics = Import[ FileNames["*guppy*data_HA.txt",folder][[1]], "csv"][[2]]/GS;
meanProfile = Import[ FileNames["*guppy*meanProfile_HA*",folder][[1]], "csv"];
meanProfile[[;;, 1]] = meanProfile[[;;,1]]/GS;
means = Import[ FileNames["*guppy*means_HA*",folder][[1]], "csv"];
means[[2;;]] = means[[2;;]]/GS;
calibrationImage = Import[ FileNames["*guppy*profile_HA.tiff",folder][[1]]];
{Grid@Transpose@{header[[{1, 6,8}]], parseFileName[file][[{1, 6,8}]]},

ListLinePlot[pressure[[;;,{1,2}]], ImageSize->Medium, Frame->True, FrameLabel->{"Time (s)", "Pressure (mBar)"}, PlotRange->{All, {0, 7000}}],
Row[{
Grid[Transpose@{{"Mean center (\[Micro]m)", "Mean STDEV (\[Micro]m)", "STDEV center (\[Micro]m)", "STDEV STDEV (\[Micro]m)", "Composite center (\[Micro]m)", "Composite STDEV (\[Micro]m)"}, guppyDiagnostics}]
}],ListLinePlot[meanProfile, ImageSize->325, FrameLabel->{"Distance (\[Micro]m)", "Frequency"},PlotLabel->"Mean profile",  Frame->True, PlotRange->{{0, 350}, {0, 1}}],
Show[calibrationImage, ImageSize->280],
Grid[means, Frame->All]
}
]


(* ::Text:: *)
(*interface for plots*)


(* ::Text:: *)
(*get plot, for use with plot interface*)


(* ::Text:: *)
(*findFiles finds all files for a given combo of sio/acetone/voltage/speed*)


FindFiles[sio_, a_, v_, s_]:=Module[{f,g, o},
f =FileNames[StringJoin["*guppy*SiO", ToString[sio], "*A", ToString[a], "*V", ToString[v],"*S", ToString[s], "_*data*"], "*", Infinity];
o = Ordering[parseFileName[#][[8]]&/@f];
f = f[[o]];
Grid@Transpose@(getStatsHA[#]&/@f)
]


(* ::Text:: *)
(*findFilesInterface finds files for a given combination of silica, acetone, voltage, and speed, and displays their data*)


findFilesInterface:=Manipulate[
FindFiles[sio4, acetone, voltage, speed]
,{{sio4, 5}, ALLSIOVALS}
,{{acetone, 8}, ALLAVALS}
, {{voltage, 25}, ALLVVALS}
,{{speed, 5}, ALLSVALS}
]


(* ::Text:: *)
(*finds the index of videos for a given combination of silica, acetone, voltage, and speed*)


FindVideoInList[sio_, a_, v_, s_]:=Module[{str},
str = StringJoin["*guppy*SiO", ToString[sio], "*A", ToString[a], "*V", ToString[v],"*S", ToString[s], "_*"];
Position[video ,y_String/; StringMatchQ[y,str]]
];


(* ::Text:: *)
(*interface to find indices of videos in the list*)


FindVideoInListInterface:=Manipulate[Module[{indices},
indices = First/@FindVideoInList[sio1,a1, v1, s1];
Grid[Transpose[{indices,video[[indices]]}], Frame->All, FrameStyle->Gray]]
 , {{sio1,5}, Append[ALLSIOVALS, "*"]}
, {{a1, 8}, Append[ALLAVALS,"*"]}
, {{v1, 25}, Append[ALLVVALS,"*"]}
, {{s1, 5}, Append[ALLSVALS,"*"]}];
