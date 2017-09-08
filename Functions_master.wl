(* ::Package:: *)

(* ::Title:: *)
(*Functions for analyzing Epon 828 - conductive glass sphere print lines*)


(* ::Section:: *)
(*Initialization: establish global variables*)


BASEDIR = "C:\\Users\\Leanne\\Documents\\Research\\Data\\Epoxy Map";
	(*BASEDIR is the directory where all of your files are stored*)
PRESFILE = "C:\\Users\\Leanne\\Box Sync\\Rheology\\Epon 828\\Pressures.txt";
If[FileExistsQ["G R ILM K compiled plus.xls"], 
	ALLSTATSMORE = Import["G R ILM K compiled plus.xls", "xls"][[1]];
	(*all stats with manually entered data*)
	HEADER = ALLSTATSMORE[[1]];
		(*header to allstats*)
	COMPILED = True;
	,
	COMPILED = False;
];
directories = Select[FileNames["*LF*CSph*P*", BASEDIR, Infinity], DirectoryQ[#]&];
	(*samples*)


<<ANOVA` (*needs anova packages*)
Needs["ErrorBarPlots`"] (*need to run this line to make plots with error bars*)


columnIndex[v_]:=Module[{}, Transpose@{Range[Length[v]], v}]


columnIndex[HEADER]
	(*list of indices and the variables they correspond to - run this to figure out what numbers to use for your associations*)


(* ::Text:: *)
(*Associations are basically switch statements or hash tables. Calling nameAssoc[2] would return "Silica w%"*)


nameAssoc = <|2-> "silica w%" , 3-> "acetone w%", 4->"particle w", 5-> "amplitude", 6->"pressure", 7-> "speed", 
						21->"shear rate", 25->"flowing viscosity", 61->"static viscosity", 13->"silica v%",14->"acetone v%" ,15->"particle v%", 62->"speed * viscosity"|>; 
viableColors = Select[columnIndex[CountDistinct/@Transpose[ALLSTATSMORE[[2;;]]]], #[[2]]<7&][[;;,1]];
	(*column to category*)
unitAssoc = <|2->" w% silica", 3-> " w% acetone", 4->" mg", 5-> " V",6-> " mBar", 7 -> " mm/s", 21->"Hz", 25->"Pa s", 61->"Pa s",
							13->" v% silica", 14->" v% acetone",15->" v% particles", 62->"mm * Pa"|>; 
	(*column to unit, used for graphing*)
nmAssoc = <|2->sil, 3->ace, 4->par,5->amp, 6->pre, 7->spe, 21->shr, 25->vis, 61->svi, 13->silvol, 14->acevol,15->parvol, 62->spevis|>; 
	(*column to variable name, used by ANOVA*)
short2LongAssoc = <|"sil"->"silica", "ace"->"acetone", "amp"->"amplitude","par"->"Particle w", "spe"->"speed", "silvol"->"silica volpct", 
						"acevol"->"acetone volpct", "parvol"->"particle volpct", "spevis"->"speed * viscosity"|>; 
	(*variable name to category, used by ANOVA*)
If[COMPILED, 
depAssoc = <|#->HEADER[[#]]&/@Range[26, 60]|>;
	(*column to dependent variable*)
]


fixP[s_]:=StringReplace[s, {"120 mg"->"1.7 v% particles", "700 mg"->"9 v% particles"}]


(* ::Text:: *)
(*These global variables are lists of all of the possible values *)


If[COMPILED, 
ALLSIOVALS = DeleteDuplicates[ALLSTATSMORE[[2;;, 2]]]; (*possible silica values*)
ALLAVALS = DeleteDuplicates[ALLSTATSMORE[[2;;, 3]]];(*possible acetone values*)
ALLPVALS = DeleteDuplicates[ALLSTATSMORE[[2;;, 4]]]; (*particle values*)
ALLVVALS = DeleteDuplicates[ALLSTATSMORE[[2;;, 5]]]; (*possible voltage values*)
ALLSVALS = DeleteDuplicates[ALLSTATSMORE[[2;;, 7]]]; (*possible speed values*)
];
SetOptions[ListLinePlot, Frame->True, FrameStyle->Black, AspectRatio->1, LabelStyle->{20, Black}, ImageSize->400];
SetOptions[ListPlot, Frame->True, FrameStyle->Black, AspectRatio->1, LabelStyle->{20, Black}, ImageSize->400];


(* ::Text:: *)
(*Here I made a mistake and did not check the magnification before taking some of my keyence images, so the scale was different and I had to do this hacky stuff. *)


ks200 = 0.901;
ks150 = 0.75*0.901;
ks100 = 0.5*0.901; (*px/micron*)
keyenceScale[d_, ace_, v_]:=Switch[d
		,"26", If[v=="25", ks150, ks200]
		,"28", ks150
		,"31", If[ace!="24", ks200, If[v=="25", ks150, ks200]]
		,"36", ks200
		,"45", ks150
		,"59", ks100
		,"66", ks100];
GS = 1.15;


getKS[file_]:=keyenceScale@@(parseDirectoryName[DirectoryName[file]][[{1,3,5}]]);


(* ::Subsubsection:: *)
(*stats for histograms where column 1 is position and column 2 is intensity*)


(* ::Text:: *)
(*justgetMean just gets the mean from a frequency plot*)


justgetMean[p_]:=Module[{n}, 	
	n = Total[p[[1 ;; All,2]]];	
	If[n>0,
		getMean[p, n]  (*mean location*), Print["n<0", 0]
]]


(* ::Text:: *)
(*getMean gets the mean location for histogram, given a total intensity: p is the profile in format {position, intensity} and n is the total intensity across all columns*)


getMean[p_, n_] := Total[p[[1 ;; All,1]]*p[[1 ;; All,2]]]/n


(* ::Text:: *)
(*getSD gets the mean and standard deviation of a frequency plot*)


getSD[p_] := Module[{n, m, fm2}, 
Catch[
	n = Total[p[[1 ;; All,2]]]; (*total intensity across all columns*)
	If[n>0,
		m = getMean[p, n];  (*mean location*)
	    fm2 = Total[p[[1 ;; All,2]]*(p[[1 ;; All,1]] - m)^2];  (*standard deviation step 1*)
		Return[N@{m, Sqrt[fm2/n]}];  (*return center and standard deviation*)
	,
		Return[{0,0}];
	]
]
]; 


(* ::Subsubsection:: *)
(*general tools*)


(* ::Text:: *)
(*parseDirectoryName splits a folder name for a sample into a vector of stats*)


parseDirectoryName[file_]:=Module[{x, lf, sio, a , CSph, V, P, S, t},
x = StringSplit[FileNameSplit[file][[-1]], "_"];
lf = StringTake[x[[1]], 3;;];
sio = StringTake[x[[2]], 4;;];
a = StringTake[x[[3]], 2;;];
CSph = StringTake[x[[4]], 5;;];
V = StringTake[x[[5]], 2;;];
P = StringTake[x[[6]], 2;;];
S = StringTake[x[[7]], 2;;];
t = StringTake[StringJoin[Riffle[x[[8;;]], "_"]], 2;;];

{lf, sio, a, CSph, V, P, S, t}
]


(* ::Section:: *)
(*compiling and exporting spreadsheets*)


(* ::Text:: *)
(*getAllDirectoryStats gets all available information about a given directory*)


getAllDirectoryStats[d_, pressures_]:=Module[{guppy, dInfo, pInfo, ILM, r, gStats, ILMG, ILMTrue, k1, k2,ILMstdev,ks,
														guppyTrue, ps, k , kTrue, ILMK, dInfo2, area, width},
Check[
	dInfo = parseDirectoryName[d]; (*independent variables stored in the file name*)
	ks = keyenceScale[dInfo[[1]], dInfo[[3]], dInfo[[5]]]; (*keyence scale bar for this sample*)
(*	dInfo2 = ToExpression[dInfo[[1;;7]]]; (*independent variables to store in the file*)*)

	(*-----------------------guppy------------------------*)
	guppy = FileNames["*guppy*data_RAW*", d];
	If[Length[guppy]>0, 
		guppyTrue = True;
		gStats = Import[guppy[[1]], "csv"][[2]];
		gStats = gStats/GS; (*convert from pixel to micron*)
		gStats = ToExpression[gStats];
		,
		guppyTrue = False;
		gStats = ConstantArray["", 6]
	];
	(*--------------------------rheology--------------------------*)
	pInfo = Select[pressures, #[[1]]==ToExpression@dInfo[[2]] && #[[2]]==ToExpression@dInfo[[3]] && #[[3]]==ToExpression@dInfo[[7]]&];
	If[Length[pInfo]>0, 
			pInfo = pInfo[[1, 4;;]];
			pInfo = ToExpression[pInfo];
		, 
			pInfo = ConstantArray["", 10]
	];
	(*--------------------------ILM--------------------------------*)
	ILM = FileNames["*ILM*stats.csv*", d];
	If[Length[ILM]>0, 
		ILMTrue = True;
		ILM = Import[ILM[[1]], "csv"][[-1]];
		ILM = ToExpression[ILM][[5;;6]];
		ILMstdev =  ILM[[2]];
		,
		ILMTrue = False;
		ILM = ConstantArray["", 2]
	];
	(*--------------------------keyence----------------------------------*)
	k = FileNames["*K*stat_f.txt", d];
	If[Length[k]>0,
		kTrue = True;
		(*mean profile*)
		k = Mean[(Import[#, "csv"][[;;,2]])&/@k];
		width = k[[5]];
		area = k[[3]];

		(*all profiles*)
		k1 = Import/@FileNames["*K*allprofsStats.csv", d];
		k2 = {Mean[k1[[;;,1,1]]], Mean[k1[[;;,1,4]]], Mean[k1[[;;,1,5]]], Mean[k1[[;;,2,1]]], Mean[k1[[;;,2,4]]], Mean[k1[[;;,2,5]]]}*{1.,1./ks, 1./ks, 1, 1./ks, 1./ks};
		k = Join[k, k2];
		,
		kTrue = False;
		k = ConstantArray["", 16];
	];
	(*---------------------------ILM, guppy------------------------------*)
	If[ILMTrue && guppyTrue, 
(*		ps = (Log10[dInfo2[[7]]]+1)/dInfo2[[7]]*dInfo2[[6]]/pInfo[[8]];*)
		ILMG = {ILMstdev/gStats[[6]](*, ILMstdev/ps, ILMstdev/(ps*gStats[[6]])*)};
		,
		ILMG = {""(*, "", ""*)};
	];
	(*-----------------------------ILM, keyence------------------------------*)
	If[ILMTrue && kTrue,
		If[guppyTrue, 
			ILMK = {ILMstdev/width, ILMstdev/area, ILMstdev/(width*gStats[[6]]/350), ILMstdev/(area*gStats[[6]])};
			,
			ILMK = {ILMstdev/k[[2]], "", "", ""};
		];
		, ILMK = ConstantArray["", 4];
	];
	(*----------------put it all together-----------------------*)			
	r = Join[dInfo,pInfo, gStats,  ILM, ILMG, k, ILMK];
Return[r];
,
Print[d];]

]


compileStatswithILM:=Module[{pressures, directories, t, ALLSTATS},
	pressures = Import[PRESFILE, "csv"];
	directories = Select[FileNames["*CSph*P*S*t*", BASEDIR, Infinity], DirectoryQ[#]&];
	ALLSTATS = Prepend[getAllDirectoryStats[#, pressures]&/@directories, 
							{"Notebook page (LF)", 
								"Silica w%", "Acetone w%", "Particle w%", 
								"Amplitude (V)", "Pressure (mBar)", "Speed (mm/s)", "Timestamp (YY_MM_DD_HH_MM_SS)",
								"stress constant","stress exponent",
								"viscosity constant","viscosity exponent","apparent shear rate","wall shear rate","shear stress (Pa)",
								"\[CapitalDelta]\.94P (mBar) stress approximation","\[CapitalDelta]\.94P (mBar) viscosity approximation","viscosity (Pa s)",
								"Guppy Mean center (\[Micro]m)", "Guppy Mean STDEV (\[Micro]m)", "Guppy STDEV center (\[Micro]m)", 
								"Guppy STDEV STDEV (\[Micro]m)","Guppy Composite Center (\[Micro]m)", "Guppy Composite STDEV (\[Micro]m)", 	
								"ILM composite center (\[Micro]m)", "ILM composite STDEV (\[Micro]m)", 
								"Out/In Composite STDEV", 
								"Height (\[Micro]m)", "Full Width (\[Micro]m)", "Full area", "Full width/height", "Full area/height",
								"Full width half max", "Full area half max", "FWHM/height", "FAHM/height", "STDEV based width","width (stdev)/height",
								"Mean height", "Mean width", "Mean width/height", 
								"STDEV height","STDEV width", "STDEV height/width", 
								 "Out STDEV/K width", "Out STDEV/Area", "Out STDEV/K width/In STDEV", "Out STDEV/Area/In STDEV"}
							 ] ;
	Export["G R ILM K compiled RAW "<>DateString[{"YearShort", "_", "Month", "_", "Day", "_", "Hour", "_", "Minute"}]<>".csv", ALLSTATS];
]


(* ::Section:: *)
(*interfaces*)


(* ::Subsubsection:: *)
(*get stats for one combo*)


getProfiles[folder_]:=Module[{dInfo, ks, width, m, nameDisplay, middle, numProfiles, middles, meanProfile, ILMprofile, Kprofiles},
dInfo = parseDirectoryName[folder];
ks = keyenceScale[dInfo[[1]], dInfo[[3]], dInfo[[5]]];
width = 1200;
m = width/2;
nameDisplay = {Grid@Transpose@{HEADER[[{1, 6,8}]], dInfo[[{1, 6,8}]]}};
If[Length[FileNames["*guppy*", folder]]>0, 
		meanProfile = Import[ FileNames["*guppy*meanProfile_RAW.csv",folder][[1]]][[2;;]];
		meanProfile[[;;, 1]] = meanProfile[[;;,1]]/GS;
		middle = justgetMean[meanProfile];
		meanProfile[[;;,1]] = (#+(m-middle))&/@meanProfile[[;;,1]];
	,
		meanProfile = {};
];
If[Length[FileNames["*ILM*filt*", folder]]>0,
	ILMprofile = ToExpression@Import[FileNames["*ILM*profiles.csv", folder][[1]]][[2;;]];
	ILMprofile = Transpose[{ILMprofile[[1,;;,1]], Rescale[Total[ILMprofile[[;;,;;,2]]]]}];
	middle = justgetMean[ILMprofile];
	ILMprofile[[;;,1]] = (#+(m-middle))&/@ILMprofile[[;;,1]];
	,
	ILMprofile = {};
];
If[Length[FileNames["*K*", folder]]>0,
	Kprofiles = Transpose[{1./ks*Range[1181], Import[#, "csv"][[;;,1]]}]&/@FileNames["*K*profile.txt", folder];
	numProfiles = Length[Kprofiles];
	middles = ConstantArray[0, numProfiles];
	Do[
		middles[[i]] = justgetMean[Kprofiles[[i]]];
		Kprofiles[[i, ;;, 1]] = (#+(m-middles[[i]]))&/@Kprofiles[[i, ;;, 1]];
		Kprofiles[[i, ;;, 2]] = Rescale[Kprofiles[[i, ;;,2]]];
	,{i, numProfiles}];
	Kprofiles = Transpose[{Kprofiles[[1,;;,1]], Rescale[Mean[Kprofiles[[;;,;;,2]]]]}];
	,
	Kprofiles = {};
];

ListLinePlot[{ ILMprofile, meanProfile,Kprofiles}, PlotLegends->{"Particles in line","Particles in channel",  "Line profile"}, AspectRatio->1/3, ImageSize->300, FrameLabel->{"Width (\[Micro]m)", "Intensity"}, PlotRange->{{0, width}, {0,1}}]
]


getILMStatsHA[folder_, gP_, gG_, gI_, gK_, gILMPL_, gILMF_, gcomp_]:=Module[{dInfo,files, pressure, guppyDiagnostics, meanProfile, calibrationImage,means,  ILMprofile, profs, Kyprofile,ks,
												ILMcalibrationImage, ILMstats,nameDisplay, pDisplay, gDisplay, iDisplay, Kprofiles, Kstats, kDisplay, ILM1, ILM2, ps, profstats, ip, cdisplay,
														getILM, k, k1, k2, numProfiles},
(*files = FileNames["*",folder];*)
dInfo = parseDirectoryName[folder];
ks = keyenceScale[dInfo[[1]], dInfo[[3]], dInfo[[5]]];
ip = {{50, 5}, {40, 5}};
nameDisplay = {Grid@Transpose@{HEADER[[{1, 6,8}]], dInfo[[{1, 6,8}]]}};
If[gP, 
	pressure = Import[ FileNames["*data*csv*",folder][[1]], "csv"];
	pDisplay = {ListLinePlot[pressure[[;;,{1,2}]], ImageSize->Medium, Frame->True, FrameLabel->{"Time (s)", "Pressure (mBar)"}, PlotRange->{All, {0, 7000}}, AspectRatio->1/2]};
	,
	pDisplay = {""};
];
If[gG && Length[FileNames["*guppy*", folder]]>0, 
		guppyDiagnostics = Import[ FileNames["*guppy*data_RAW.csv",folder][[1]]][[2]]/GS;
(*		meanProfile = Import[ FileNames["*guppy*meanProfile_RAW.csv",folder][[1]]];
		meanProfile[[;;, 1]] = meanProfile[[;;,1]]/GS;*)
		calibrationImage = Import[ FileNames["*guppy*profile_RAW.tiff",folder][[1]]];
		gDisplay ={Style[ "Guppy", Large]
(*					, 
					ListLinePlot[meanProfile, ImageSize->325, FrameLabel->{"Distance (\[Micro]m)", "Frequency"},
							PlotLabel->"Mean profile",  Frame->True, PlotRange->{{0, 350}, {0, 1}}]*)
					,
					Show[calibrationImage, ImageSize->280]
					,
					Grid[Transpose@{{"Mean center (\[Micro]m)", "Mean STDEV (\[Micro]m)", "STDEV center (\[Micro]m)", 
								"STDEV STDEV (\[Micro]m)", "Composite center (\[Micro]m)", "Composite STDEV (\[Micro]m)"}, guppyDiagnostics}]
				};
	,
	gDisplay = ConstantArray["", 3];
];
If[gI && Length[FileNames["*ILM*filt*", folder]]>0,
	getILM = True;
(*	ILMcalibrationImage = Import[FileNames["*ILM*fits.tiff", folder][[1]]];
	ILMcalibrationImage = ImageCrop[ILMcalibrationImage, ImageDimensions[ILMcalibrationImage]*{1/7, 1}];*)
	If[gILMF, ILM1 = Import[FileNames["*ILM*filt.bmp", folder][[2]]]; ILM1 = ImageResize[ImageCrop[ILM1, {926, 1200}], 400];, ILM1=Graphics[]];
	If[gILMPL, ILM2 = Import[FileNames["*ILM*L*L*tif*", folder][[2]]]; ILM2 = ImageAdjust@ImageResize[ImageCrop[ILM2, {926, 1200}], 400], ILM2 = Graphics[]];
(*	ILMprofile = ToExpression@Import[FileNames["*ILM*profiles.csv", folder][[1]]];*)
	ILMstats = Import[FileNames["*ILM*stats.csv", folder][[1]]][[{1, -1},-1]];
	iDisplay = {Style["ILM", Large]
				,
				Row[ILMstats]
				,
				If[gILMF, Show[ILM1, ImageSize->400, ImagePadding->ip], ""]
				,
				If[gILMPL, Show[ILM2, ImageSize->400, ImagePadding->ip], ""]
(*				,
				ListLinePlot[Transpose[{ILMprofile[[1,;;,1]], Total[ILMprofile[[;;,;;,2]]]/Max[Total[ILMprofile[[;;,;;,2]]]]}], PlotStyle->Black, 
										AspectRatio->225/1200, ImageSize->400, PlotRange->{{100, 1300}, {0, 1}}, ImagePadding->ip, FrameLabel->{"Width (\[Micro]m)", "Intensity"}]*)
				};
	,
	iDisplay = ConstantArray["", 4];
	getILM = False;
];
If[gK && Length[FileNames["*K*", folder]]>0,
(*	Kprofiles = Import[#]&/@FileNames["*K*profile.tiff", folder];*)
		k = (Import[#, "csv"][[;;,2]])&/@FileNames["*K*stat_f.txt", folder];
		numProfiles = Length[k];

		k1 = Import/@FileNames["*K*allprofsStats.csv", folder];
		k2 = ({1.,1./ks, 1./ks, 1, 1./ks, 1./ks}*Flatten[#,1])&/@k1[[;;,{1,2},{1,4,5}]];



(*	Kstats = (Import[#, "csv"][[{1,4},1]]*{1., 1./ks})&/@FileNames["*K*stats.txt", folder];
	profstats = Flatten[((Import[#, "csv"][[{1,2},{1,4,5}]])*ConstantArray[{1., 1./ks, 1./ks},2]), 1]&/@FileNames["*K*allprofsStats.csv", folder];*)
	Kstats = Flatten/@Transpose[{k, k2}]; 
	profs = FileNames["*K*allprofs.csv", folder];
	ps = {{Red}, {Red, Dashed}, {Black}, {Black, Dashed}, {Blue}, {Blue, Dashed}}[[1;; 2*numProfiles]];
	Kyprofile = Flatten[Transpose[({1., 1./ks}*#)&/@Import[#, "csv"][[;;,{1,4}]]]&/@profs,1];
(*	Kyprofile = Transpose[{1./ks*Range[801], #}]&/@Kyprofile;*)
	Kyprofile = ListLinePlot[Kyprofile, ImageSize->(400-ip[[1,1]]-ip[[1,2]])*800/1181+ip[[1,1]]+ip[[1,2]], FrameLabel->{"Length (px)", "Height (thick) / Width (dotted) (\[Micro]m)"},
										AspectRatio->1, ImagePadding->ip,PlotStyle->ps, PlotRange->{{(*(801/ks-801)/2, 801/ks-(801/ks-801)/2*)0,800}, {0, 800}}];
	Kprofiles = Transpose[{1./ks*Range[1181], Import[#, "csv"][[;;,1]]}]&/@FileNames["*K*profile.txt", folder];
	Kprofiles = ListLinePlot[Kprofiles, ImageSize->400, AspectRatio->225/1200, FrameLabel->{"Width (\[Micro]m)", "Height (\[Micro]m)"},
										PlotRange->{{(1181/ks-1181)/2, 1181/ks-(1181/ks-1181)/2}, {0, 225}}, PlotStyle->Partition[ps,2][[;;,1]],ImagePadding->ip];
	kDisplay = {Style["Keyence", Large]
				,
				Kprofiles (*x cross-section*)
				,
				Kyprofile (*y cross-section*)
				,
				Grid[Transpose@Prepend[Append[Kstats, Mean[Kstats]], {"Height (\[Micro]m)", "Full Width (\[Micro]m)", "Full area", "Full width/height", "Full area/height",
								"Full width half max", "Full area half max", "FWHM/height", "FAHM/height", "STDEV based width",
								"STDEV based width/height", "Mean height", "Mean width", "Mean width/height", 
								"STDEV height","STDEV width", "STDEV height/width"}], 
																		ItemStyle->{{Bold, Black, Red, Blue}[[1;;Length[profs]+1]], None}]
				,
				If[getILM,
						Row[{Style["Focused width/line width\t", Bold], ILMstats[[2]]/Mean[Kstats[[;;,10]]]}], ""]};
	,
	kDisplay = ConstantArray["", 5];
];
If[gcomp, cdisplay = {getProfiles[folder]}, cdisplay = {""};];
Join[nameDisplay, pDisplay, gDisplay, iDisplay, kDisplay, cdisplay]

]


(* ::Text:: *)
(*srch fixes search terms so they're appropriate strings*)


srch[s_]:=If[s==6.5, ToString[s], StringReplace[ToString[s], "."->""]];


(* ::Text:: *)
(*findDirectories finds all directories for a given combo of sio/acetone/voltage/speed and gets their stats*)


FindDirectories[sio_, a_, v_, s_,pt_ ,gP_, gG_, gI_, gK_, gILMPL_, gILMF_, gcomp_]:=Module[{f,g, o, str},
str = StringJoin["*SiO", srch[sio], "_A", srch[a], "_CSph", srch[pt], "_V", srch[v],"*S", srch[s], "_t*"];
f =Select[FileNames[str, BASEDIR, Infinity], DirectoryQ[#]&];
If[Length[f]>0,
o = Ordering[parseDirectoryName[#][[8]]&/@f]; (*puts the files in time order*)
f = f[[o]];
Grid@Transpose@(getILMStatsHA[#, gP, gG, gI, gK, gILMPL, gILMF, gcomp]&/@f)
	,
	str]
]


(* ::Text:: *)
(*findFilesInterface finds files for a given combination of silica, acetone, voltage, and speed, and displays their data*)


ss = 30;


findDirectoriesInterface:=Manipulate[

FindDirectories[sio4, acetone, voltage, speed,pt,  getPressure, getGuppy, getILM, getKeyence, getILMPla, getILMFilt, getComposite]
,Row[{Control[{{sio4, 5}, ALLSIOVALS}]
,Spacer[ss]
,Control@{{acetone, 8}, ALLAVALS}
,Spacer[ss]
,Control@{{voltage, 25}, ALLVVALS}
,Spacer[ss]
,Control@{{speed, 5}, ALLSVALS}
,Spacer[ss]
,Control@{{pt, 120}, ALLPVALS}
,Spacer[ss]
,Control@{{getPressure, False}, {True, False}}
,Spacer[ss]
,Control@{{getGuppy, False}, {True, False}}
,Spacer[ss]
,Control@{{getILM, False}, {True, False}}
,Spacer[ss]
,Control@{{getILMPla, False}, {True, False}}
,Spacer[ss]
,Control@{{getILMFilt, False}, {True, False}}
,Spacer[ss]
,Control@{{getKeyence, False}, {True, False}}
,Spacer[ss]
,Control@{{getComposite, False}, {True, False}}
}]
]


(* ::Subsubsection:: *)
(*keyence profiles*)


(* ::Text:: *)
(*getArea gets the areas of all keyence profiles in a folder*)


getArea[folder_]:=Module[{k},
If[Length[FileNames["*K*", folder]]>0,
		k = ({#, Import[#, "csv"][[3,2]]})&/@FileNames["*K*stat_f.txt", folder]
,
{}
]
]


(* ::Text:: *)
(*makeStr makes a set of strings to use with FileNames to search for files in this folder*)


makeStr[sio_, a_, v_, s_,pt_ ]:=Module[{variables, varied, str},
Catch[
variables = {sio, a, pt, v, s};
varied = Position[variables, {__,__}];
Switch[Length[varied]
,1,
	str = {StringJoin["*SiO", srch[sio], "_A", srch[a], "_CSph", srch[pt], "_V", srch[v],"*S", srch[s], "_t*"]};
	Return[str];
,2,
	varied = varied[[1,1]];
	str = StringJoin["*SiO", If[varied==1, srch[#],srch[sio]], "_A", If[varied==2, srch[#],srch[a]], "_CSph", If[varied==3, srch[#],srch[pt]], "_V", If[varied==4, srch[#],srch[v]],"*S", If[varied==5, srch[#],srch[s]], "_t*"]&/@variables[[varied]];
	Return[str];
,_,
Return["Too many variables"]];
]]


(* ::Text:: *)
(*getSimilarProfiles finds the set of profiles with the most similar areas*)


getSimilarProfiles[sio_, a_, v_, s_,pt_ ]:=Module[{str, f, KAreas, groups, mostSimilar, dInfo, variables, varied, numVals},
Catch[
str = makeStr[sio, a, v, s, pt]; (*strings to search*)
If[Length[str]>0, 
varied = variedIndex[{sio, a, pt, v, s}];
f = Select[FileNames[str, BASEDIR, Infinity], DirectoryQ[#]&]; (*directories*)
numVals = Length[(DeleteDuplicates/@Transpose[(parseDirectoryName/@f)[[;;, {2,3,4,5,7}]]])[[varied]]];
If[numVals>1 || varied<1,
	KAreas = Flatten[DeleteCases[getArea/@f,{}],1]; (*file names and areas*)
	Switch[Length[str]
	,1, (*if makeStr only returned one directory type to query - this happens when the search set only has one value per variable*)
		If[Length[KAreas]>1, {KAreas[[;;,1]], StandardDeviation[KAreas[[;;,2]]]/Mean[KAreas[[;;,2]]]}, {KAreas[[;;,1]], 0}] (*return all k files in the eligible directories*)
	,_, (*if makeStr returned several directory types to query*)
		varied = Switch[varied, 1, "_SiO", 2, "_A", 3, "_CSph", 4, "_V", 5, "_s"]; (*string to split at*)
		groups = Tuples[GatherBy[KAreas, StringSplit[StringSplit[FileBaseName[#[[1]]],varied][[2]], "_"][[1]]&]]; (*all combinations of files that encompass all of the variables queried*)
		groups = {#[[;;,1]], StandardDeviation[#[[;;,2]]]/Mean[#[[;;,2]]]}&/@groups;
		mostSimilar = SortBy[groups,  Last][[1;;Min[10, Length[groups]]]]; (*sort groups of files by similarity to each other*)
		mostSimilar[[;;,1]] = putFilesInOrder/@(mostSimilar[[;;,1]]);(*sort the files to put the varied variable in order*)
		mostSimilar	
	]
,
{{{""}, 1000000}}]
,
	{{{""},1000000}}
]]]


putFilesInOrder[files_]:=SortBy[files, ToExpression[parseDirectoryName[DirectoryName[#]][[2;;5]]]&]


getVariedVariable[f_]:=Module[{dInfo, variables, varied},
	dInfo = parseDirectoryName/@f; (*directory stats from folders*)
	variedIndex[DeleteDuplicates/@Transpose[dInfo[[;;, {2,3,4,5,7}]]]]
]


variedIndex[variables_]:=Module[{varied},
	varied = Position[variables, {__,__}]; (*the index of the varied variable*)
	Switch[Length[varied],1, -1,2, varied[[1,1]],_, 0] (*if 1, then there aren't enough variables. if 2, then just enough variables. if more, then too many variables*)
]


getXY[file_, middle_]:=Module[{ks, profile, center},
ks = getKS[file];
profile = Import[file, "csv"][[;;,1]]; (*profiles as a list of heights*)
center = (Length[profile]/2-Ordering[profile,-1][[1]]); (*maximal points = centers*)
Transpose[{(Range[1, Length[profile]]+center)/ks+middle-Length[profile]/2/ks, profile}]
]


plotTogether[l_]:=Module[{width, middle, files, profiles, numvars, wh, folders, dInfo, variables, varied, v2u, unit, title, legend, colors},
(*CONSTANTS*)
width = 800;
middle = width/2;

(*DATA*)
files = StringReplace[#, "stat_f"->"profile"]&/@(l[[1]]); (*profile files to import*)
profiles = getXY[#, middle]&/@files;
numvars = Length[files];
wh = (Import[#, "csv"][[11,2]])&/@(l[[1]]); (*widths/heights*)

(*LABELING*)
folders = DirectoryName/@(l[[1]]); (*folders that the profiles rest in*)
dInfo = parseDirectoryName/@folders; (*directory stats from folders*)
variables = DeleteDuplicates/@Transpose[dInfo[[;;, {2,3,4,5,7}]]]; (*list of queried variables*)
varied = variedIndex[variables];
v2u=<|-1->"", 0->"", 1->2, 2->3, 3->4, 4->5, 5->7|>;
unit = v2u[varied];
title = StringJoin[Riffle[ToString[variables[[#, 1]]]<>unitAssoc[v2u[#]]&/@DeleteCases[Range[1,5], varied],", "]];
title = fixP[title];
If[varied>0, 
		legend = ToString[variables[[varied,#]]]<>unitAssoc[unit]<>", width/height = "<>ToString[wh[[#]]]&/@Range[numvars]
		, 
		legend = StringTake[FileBaseName/@files, 3;;]
];
colors = ColorData["DarkRainbow"][#]&/@(Range[0, numvars-1]/(numvars-1));
legend = fixP/@legend;

ListLinePlot[Append[profiles, {{0,0},{middle-175,0}, {middle-174.9999, 150}, {middle+175, 150}, {middle+175.00001, 0}, {width,0}}] , 
						PlotRange->{{0, width}, {0, 175}},
						ImageSize->800,
						AspectRatio->175/width, 
						FrameLabel->{"Width (\[Micro]m)", "Height (\[Micro]m)"}, 
						PlotLegends->Append[legend, "Channel"], 
						PlotLabel->title,
						PlotStyle->Append[colors, Black]]
]


findProfilesSwitch[x_, default_]:=Switch[Length[x], 0, default, 1, x, _, {x}];


FindProfiles:=Manipulate[
Module[{tuples, idx},
	update;
	idx = variedIndex[{sio1, a1, v1, s1, p1}];
	Column[{
	Button["Update", update = 1-update],
	Switch[idx
		,0, "Too many variables"
		,_,
		tuples = Tuples[{findProfilesSwitch[sio1, ALLSIOVALS], findProfilesSwitch[a1, ALLAVALS],findProfilesSwitch[v1, ALLVVALS],findProfilesSwitch[s1, ALLSVALS],findProfilesSwitch[p1, ALLPVALS]}];
		profs = Flatten[getSimilarProfiles[#[[1]], #[[2]], #[[3]], #[[4]],#[[5]]]&/@tuples,1];
		profs = Select[Sort[profs, #1[[2]]<#2[[2]]&][[1;;Min[10, Length[profs]]]], #[[2]]<0.2&];
		pics = Column[{plotTogether[#], 
					Row[{"\tstdev area/mean area = ", #[[2]]}]
				}]&/@profs;
		If[Length[profs]>0, 
			Row[Column[{pics[[#]], Button["Export "<>makeFileName[profs[[#]]], Export[makeFileName[profs[[#]]],pics[[#]]]]}]&/@Range[1, Length[pics]]]
			,
			""]
		]
	
	}]
]
, {{sio1, {5.}}, Append[Join[{#}&/@ALLSIOVALS, {ALLSIOVALS, {5., 6., 7., 8., 9., 10.},{5., 6., 6.5, 8., 9., 10.}, {5., 6.,7., 9., 10.}, {5., 7., 9.}, {6., 8., 10.}}], 0->"Any value"]}
, {{a1, 0}, Append[DeleteCases[Subsets[ALLAVALS], {}], 0->"Any value"]}
, {{v1, {50.}}, Append[DeleteCases[Subsets[ALLVVALS],{}], 0->"Any value"]}
, {{s1, {1.}}, Append[DeleteCases[Subsets[ALLSVALS],{}], 0->"Any value"]}
, {{p1, ALLPVALS}, Append[DeleteCases[Subsets[ALLPVALS],{}], 0->"Any value"]}
, {{update, 0}, ControlType->None}
, {{profs2, 0}, ControlType->None}
, TrackedSymbols:>{update}
]


makeFileName[group_]:=BASEDIR<>"\\plots\\profiles\\"<>cleanFiles[group[[1]]]<>".pdf";


cleanFiles[files_]:=StringRiffle[DeleteDuplicates[Flatten[StringSplit[#, "_"]&/@(StringReplace[FileBaseName[#], {"K_"->"", "_stat_f"->""}]&/@files)]], "_"]


(* ::Subsubsection:: *)
(*plot*)


(x = StringJoin[BASEDIR, "\\plots\\scatter plots\\", 
								Evaluate[StringReplace[#, {"/"->" over "}]]]; If[!DirectoryQ[x], CreateDirectory[x]];)&/@(depAssoc[#]&/@Keys[depAssoc]);


(* ::Text:: *)
(*get filtered data*)


PIgetGrid[sio1_, a1_, v1_, s1_,p1_, x_, c_, onlyCont_, maxArea_, minArea_, minwh_, maxwh_]:=Module[{sio2, a2, v2, s2, xPlot, selectCont, thisGrid, runName, p2, var},Catch[
 sio2 =sio1; a2 = a1; v2 = v1; s2 = s1; p2=p1;
	(*establish a local instance of manipulated variables*)
(*Switch[x,2, sio2 = ALLSIOVALS; ,3, a2 =ALLAVALS; ,5,v2 =ALLVVALS;,7, s2=ALLSVALS;];
Switch[c,2, sio2 = ALLSIOVALS; ,3, a2 = ALLAVALS; ,5,v2 =ALLVVALS;,7, s2= ALLSVALS;];*)
	(*establish which data to select*)
var = Switch[x, 2, sio2, 3, a2, 5, v2, 7, s2, 4, p2];
If[Length[var]<2, var = Switch[x, 2, ALLSIOVALS, 3, ALLAVALS, 5, ALLVVALS, 4, ALLPVALS]]; 
	(*determine the variable to plot on the x axis*)
If[onlyCont, selectCont = {0.}, selectCont = {0.,1.}];
	(*if we only want the continuous lines, select data where Noncontinuous = 0*)
thisGrid = Select[ALLSTATSMORE, (#[[57]]=="" || maxArea>#[[57]]>minArea )&& MemberQ[sio2*1.,#[[2]]] && 
						MemberQ[p2*1.,#[[4]]] && MemberQ[ a2*1., #[[3]]] && MemberQ[ v2*1., #[[5]]] && 
						MemberQ[ s2*1., #[[7]]] && MemberQ[selectCont, #[[9]]]&(* && (#[[37]]=="" || minwh<#[[37]]/#[[37]]<maxwh)&*)];
	(*select the data*)
If[Length[thisGrid]>0,  (*if we successfully selected data, do more stuff. if not, say nothing selected*)
	If[x!=c && c!=0,  (*if the color and x axis are different, split the data by the color variable. else, plot everything in black*)
		thisGrid = SplitBy[Sort[thisGrid,#1[[c]]<#2[[c]]&], Part[#,c]&];
		,
		thisGrid = {thisGrid}
	];
	,
	thisGrid = {};
];
runName=PIgetRunNames[sio2, a2, v2, s2,p2, maxArea, minArea, onlyCont];
Return[{thisGrid, runName, sio2, a2, v2, s2, p2}];
]]


(* ::Text:: *)
(*get best fit lines*)


PIgetLine[plotIndex_, thisGrid_, xPlot_, showMeans_, model_]:=Module[{lines, numSubs, vars},
Catch[
numSubs = Length[thisGrid];
lines = ConstantArray[0, numSubs];
	vars = Cases[#, {_?NumberQ, _?NumberQ}]&/@thisGrid[[;;,;;,{xPlot, plotIndex}]]; (*remove all empty spaces*)
	If[showMeans,
			vars = (meansNstderrs[#][[;;,1]])&/@vars;
	];
	If[Length[vars]>0 && CountDistinct[vars[[1,;;,1]]]>1, 

		Do[If[Length[vars[[j]]]>2,
			Switch[model
			,"linear",
				lines[[j]] = LinearModelFit[vars[[j]], {1,nmAssoc[xPlot]}, nmAssoc[xPlot]]; (*determine a linear fit*)
			,"quadratic",
				lines[[j]] = Check[NonlinearModelFit[vars[[j]], a+c*var^2, {a,c},var],"cannot fit model"];
			,"^1/2",
				lines[[j]] = Check[NonlinearModelFit[vars[[j]], a+c*(var)^(1/2), {a,c},var],"cannot fit model"];
			,"^-1",
				lines[[j]] = Check[NonlinearModelFit[vars[[j]], a+c*(var)^-1, {a,c}, var], "cannot fit model"];
			]
			,
			lines[[j]] = LinearModelFit[{{Min[thisGrid[[;;,;;,xPlot]]],-100}, {Mean[MinMax[thisGrid[[;;,;;,xPlot]]]],-101},{Max[thisGrid[[;;,;;,xPlot]]],-102}}, {1, var}, var]; 
			];
			, {j, 1, numSubs}];
		,
		lines = {};
	];
Return[lines];
]
]


(*PIgetLines[plotIndices_, thisGrid_, xPlot_, showMeans_]:=Module[{vars, lines, numSubs},Catch[
numSubs = Length[thisGrid];
lines = ConstantArray[0, {Length[plotIndices], numSubs}];
		(*create an array to hold the line functions*)
Do[
	lines[[i]] = PIgetLine[plotIndices[[i]], thisGrid, xPlot, showMeans];
,{i, Length[plotIndices]}];
Return[lines];
]]*)


(* ::Text:: *)
(*get names for saving files*)


PIgetRunNames[sio2_, a2_, v2_, s2_,p2_, maxArea_, minArea_, onlyCont_]:=Module[{sioName, aName, vName, sName, pName,runName},
sioName = Switch[Length[sio2], 1, ToString[sio2[[1]]], _, "all"];
aName = Switch[Length[a2], 1, ToString[a2[[1]]], _, "all"];
vName = Switch[Length[v2], 1, ToString[v2[[1]]], _, "all"];
sName = Switch[Length[s2], 1, ToString[s2[[1]]], _, "all"];
pName = Switch[Length[p2], 1, ToString[p2[[1]]], _, "all"];
runName = StringJoin["sio", sioName, "_a", aName, "_v", vName, "_s", sName, "_p", pName, "_", fileFriend[maxArea],"_", fileFriend[minArea],If[onlyCont, "_cont_","_"]];
runName
]


(* ::Text:: *)
(*standard error*)


standardError[data_]:=If[Length[data]>1,StandardDeviation[data]/Sqrt[Length[data]],0]


(* ::Text:: *)
(*gathers data into subsets by column 1*)


col1Sort[testData_]:=Module[{d}, d = Cases[testData,  {_?NumberQ, _?NumberQ}]; GatherBy[d, #[[1]]&]];


(* ::Text:: *)
(*meansNstderrs returns a table of data in the format {{col1, mean(col2)}, ErrorBar[standardError(col2)]} for use with errorlistplot*)


meansNstderrs[testData_]:=Module[{d}, d=col1Sort[testData]; ({#[[{1,2}]], #[[3]]})&/@Transpose[{d[[;;, 1,1]], Mean/@(d[[;;, ;;,2]]),  ErrorBar/@standardError/@(d[[;;,;;,2]])}]];


(* ::Text:: *)
(*get plot, for use with plot interface*)


PIgetImage[depvar_, thisGrid_, xPlot_, runName_, c_, showMeans_, fixedRange_, plotMode_, getLines_, io_, model_]:=Module[
					{xaxis,plotRanges, l, d, image, lines, colors, plotPreferences,allDataMeans, allDataErrs, f, pr, showLines, grid, numSubs, zaxis},
image = Switch[plotMode
	,0, PIgetImageSc[depvar, thisGrid, xPlot, c, showMeans, fixedRange, getLines, model]
	,1, PIgetImage3D[depvar, thisGrid, xPlot, c, showMeans, fixedRange, io]
	,2, PIgetImagecont[depvar, thisGrid, xPlot, c, fixedRange, io]
	];

image

]


depAssoc


(* ::Text:: *)
(*plot range association - update this as needed*)


prAssoc = <|24->{51, 70}, 33->{0.14, 0.216917}/0.216917, 40->{0.21, 0.39}, 42->{1.2, 2.0}, 58->{0.35, 1}|>;


(* ::Text:: *)
(*get contour plot*)


PIgetImagecont[depvar_, thisGrid_, xPlot_, c_, fixedRange_,io_]:=Module[
					{plotPreferences, xaxis, zaxis, grid, d, f, pr, image},

(*depvar is the column index of the dependent variable*)
(*thisGrid is the grid of all data, 3 levels deep*)
(*xPlot is the index of the independent variable*)
(*c is the color variable number*)
(*yAxes is the list of names of dependent variables*)
(*showMeans is a bool that indicates whether to plot all points or jsut the means and standard error*)
(*fixedRange is a bool that indicates whether to fix the plot range or just plot all*)

plotPreferences = { ImageSize->Medium,LabelStyle->{FontSize->14,FontColor->Black},ColorFunction->"RedBlueTones", PlotLegends->Automatic, TicksStyle->Black,  InterpolationOrder->io};

	(*values of independent variables for legend on color plots*)
xaxis = nameAssoc[xPlot];
zaxis = nameAssoc[c];

grid = Flatten[thisGrid,1];
d = {#[[1,1]], #[[1,2]], Mean[#[[;;,3]]]}&/@SplitBy[grid[[;;, {xPlot,c, depvar}]], {#[[1]],  #[[2]]}&]; (*data to plot*)
f = ListContourPlot; (*function to plot with*)

If[fixedRange,
	pr = {All, All, Lookup[prAssoc,depvar,All]};
	,
	pr = All];

(*image is a plot with a legend that includes ANOVA data for each color/subset*)
image = f[d,FrameLabel->{xaxis, zaxis}, PlotLabel->depAssoc[depvar], plotPreferences, PlotRange->pr];
image

]


(* ::Text:: *)
(*get 3D plot*)


PIgetImage3D[depvar_, thisGrid_, xPlot_, c_, showMeans_, fixedRange_, io_]:=Module[
					{plotPreferences, xaxis, zaxis, grid, d, f, pr, image},

(*depvar is the column index of the dependent variable*)
(*thisGrid is the grid of all data, 3 levels deep*)
(*xPlot is the index of the independent variable*)
(*c is the color variable number*)
(*yAxes is the list of names of dependent variables*)
(*showMeans is a bool that indicates whether to plot all points or jsut the means and standard error*)
(*fixedRange is a bool that indicates whether to fix the plot range or just plot all*)

plotPreferences = { ImageSize->Medium,Boxed->True, BoxStyle->Black, LabelStyle->{Black, Medium}, TicksStyle->Black, ColorFunction->"RedBlueTones"};


	(*values of independent variables for legend on color plots*)
xaxis = nameAssoc[xPlot];
zaxis = nameAssoc[c];

grid = Flatten[thisGrid,1];
If[showMeans, 
		(*if yes, show the mean*)
		d = {#[[1,1]], #[[1,2]], Mean[#[[;;,3]]]}&/@SplitBy[grid[[;;, {xPlot,c, depvar}]], {#[[1]],  #[[2]]}&]; (*data to plot*)
		f = ListPlot3D; (*function to plot with*)
		plotPreferences = Append[plotPreferences, InterpolationOrder->io];
		,
		(*else show all points*)
		d = grid[[;;, {xPlot,c, depvar}]]; (*data to plot*)
		f = ListPointPlot3D; (*function to plot with*)
		plotPreferences = Append[plotPreferences, PlotStyle->PointSize[0.02]];
];

If[fixedRange,
	pr = {All, All, Lookup[prAssoc,depvar,All]};
	,
	pr = All];


(*image is a plot with a legend that includes ANOVA data for each color/subset*)
image = f[d,AxesLabel->{xaxis, zaxis, depAssoc[depvar]}, plotPreferences, PlotRange->pr];
image

]


(* ::Text:: *)
(*get scatter plot*)


nameAssoc


PIgetImageSc[depvar_, thisGrid_, xPlot_, c_, showMeans_, fixedRange_, getLines_, model_]:=Module[{numSubs, colors, plotPreferences,showLines, l, xaxis, zaxis, d, f, pr, lines, image, shapes, ac, header},
Catch[
(*depvar is the column index of the dependent variable*)
(*thisGrid is the grid of all data, 3 levels deep*)
(*xPlot is the index of the independent variable*)
(*c is the color variable number*)
(*showMeans is a bool that indicates whether to plot all points or jsut the means and standard error*)
(*fixedRange is a bool that indicates whether to fix the plot range or just plot all*)
(*getLines is a bool that indicates whether to get best fit lines*)
numSubs = Length[thisGrid];
If[numSubs>1, colors = ColorData["DarkRainbow"][#]&/@(Range[0, numSubs-1]/(numSubs-1));, colors = {Black}];
shapes = {"\[FilledCircle]", "\[FilledSquare]", "\[FilledUpTriangle]", "\[FilledDiamond]", "\[FivePointedStar]", "\[CircleTimes]"}[[1;;numSubs]];

plotPreferences = { PlotStyle->colors,PlotMarkers->({#, 13}&/@shapes), ImageSize->250, AspectRatio->1, (*PlotLegends->Placed[l, Right], *) 
						Frame->True, FrameStyle->Black, ImagePadding->{{50,10}, {50,10}}, LabelStyle->{FontFamily->"Arial",FontSize->12,FontColor->Black}, PlotRangePadding->Scaled[0.08], TicksStyle->Black};

Quiet[If[getLines, lines = PIgetLine[depvar, thisGrid, xPlot, showMeans, model];, lines = {}]];
showLines = getLines;
If[lines=={} , showLines = False; , showLines = True;];

If[c!=xPlot && c!=0, l = StringJoin[ToString[#],unitAssoc[c]]&/@thisGrid[[;;, 1, c]];, l={""}];
	(*values of independent variables for legend on color plots*)
xaxis = nameAssoc[xPlot]; (*x label*)
zaxis = nameAssoc[c]; (*y label*)
	If[showMeans, 
		(*if yes, show the mean with error bars*)
		d = meansNstderrs/@thisGrid[[;;,;;, {xPlot, depvar}]]; (*data to plot*)
		f = ErrorListPlot; (*plot function*)
		,
		(*else show all points*)
		d = thisGrid[[;;,;;, {xPlot, depvar}]]; (*data to plot*)
		If[xPlot==25, f = ListLogLinearPlot,f = ListPlot] (*plot function*)
	];



If[fixedRange,
	pr = {All, Lookup[prAssoc,depvar,All]};
	,
	If[showMeans, pr = {MinMax@d[[;;,;;,1,1]], MinMax[{d[[;;,;;,1,2]]-d[[;;,;;,2,1]], d[[;;,;;,1,2]]+d[[;;,;;,2,1]]}]}, pr = All]]; (*manually calculate the plot range to fit error bars*)

ac = 3;
header = {"","","",ar@"R^2"};
header = If[model=="linear", Append[header,  ar@"p"], header];
(*image is a plot with a legend that includes ANOVA data for each color/subset*)
image = Column[{Show[f[d,FrameLabel-> {xaxis, depAssoc[depvar]}, plotPreferences, PlotRange->pr]
			, 
			If[showLines, Plot[#[var]&/@lines, {var, Min[flattenNclear[thisGrid[[;;,;;,xPlot]]]], Max[flattenNclear[thisGrid[[;;,;;,xPlot]]]]}, Evaluated->True, PlotStyle->colors],  Graphics[]]] (*plot*)	
			,
			If[showLines, (*Style[Grid[Prepend[Transpose[{Graphics[{#,PointSize[1], Point[{0,0}]}, ImageSize->9]&/@colors, StringReplace[#, "."->""]&/@l, 
						#["RSquared"]&/@lines}], {"","",(*"", *)  "R^2"(*, "p"*)}]], FontFamily->"Arial"] (*legend*)*)
						Grid[Prepend[{Style[shapes[[#]], colors[[#]]],
								ar[StringReplace[l[[#]], "."->""]],
								Style[SetPrecision[Normal[lines[[#]]],ac], FontFamily->"Arial"],
								If[!StringQ[Normal[lines[[#]]]], ar[SetPrecision[lines[[#]]["RSquared"],ac]], ""],
								If[model=="linear", ar[SetPrecision[(lines[[#]]["ANOVATablePValues"])[[1]],ac]], ""]} &/@Range[1, numSubs], header]] (*legend*)
(*						Transpose[Prepend[#[{"RSquared"(*, "ANOVATablePValues"*)}](*, Normal[#]*)]&/@lines]]], {"","",(*"", *)  "R^2"(*, "p"*)}]], FontFamily->"Arial"]*)
					,
					(*If[Length[thisGrid]>1, Style[Grid[{Graphics[{#,PointSize[1], Point[{0,0}]}, ImageSize->9]&/@colors, l}], FontFamily->"Arial"], ""]*)
					If[c!=xPlot && c!=0, Column[Style[shapes[[#]]<>" "<>StringReplace[l[[#]], "."->""], colors[[#]], FontFamily->"Arial"]&/@Range[1, numSubs]]] (*legend*)
			]
		}, Center];
Return[image];
]
]


flattenNclear[l_]:=DeleteCases[Flatten[l],""];


ar[x_]:=Style[ToString[x], FontFamily->"Arial"];


(* ::Text:: *)
(*plots to start with*)


plotIndices = {};


(* ::Text:: *)
(*interface that creates plots*)


plotInterface:=Manipulate[DynamicModule[{disp, fileNames,xPlot, runName, numSubs, colors, plotPreferences,
											vs,  yAxes, indexVector, display, sio2, a2, v2, s2, list, factors, 
											interactions, jjj, getSiO, getAcetone, getAmplitude, getSpeed, p2, getParts},

{thisGrid, runName,  sio2, a2, v2, s2, p2}=PIgetGrid[sio1, a1, v1, s1, p1, x, c, onlyCont, maxArea, minArea, minwh, maxwh];
sio1 = sio2; a1 = a2; v1 = v2; s1 = s2; p1 = p2;
xPlot = x;
Column[{
Row[{Button["Filter 0.8-1.5", maxArea = 1.5; minArea=0.8; onlyCont = True;],Button["Filter 0-10", maxArea = 10; minArea=0; onlyCont=True;]}],
TogglerBar[Dynamic[plotIndices], (#->depAssoc[#])&/@Keys[depAssoc], Appearance->"Row"] (*buttons with dependent variables to plot*), 
If[thisGrid!={},
	display = ConstantArray["", {Length[plotIndices], 2}];
	(*--------------PLOT MODE-------------------*)
	If[showPlot,

		indexVector = Range[1, Length[plotIndices]];
		disp = Quiet[PIgetImage[#, thisGrid, xPlot,  runName, c, showMeans, fixedRange, plotMode, showLines, io, model]&/@plotIndices];
			(*get the plots for each dependent variable*)
		fileNames = StringJoin[BASEDIR, "\\plots\\scatter plots\\", 
								Evaluate[fileFriend[depAssoc[#]]], "\\",
								Evaluate@runName, 
								Evaluate@ToString@nmAssoc[xPlot],"_", 
								Evaluate@ToString@nmAssoc[c],
								If[showMeans, "_mean", "_pt"],
								If[fixedRange, "_fr", ""], 
								Switch[plotMode, 0, "_scat", 1, StringJoin["_3d_", ToString[io]], 2, StringJoin["_cont", ToString[io]]], 
								If[showLines, "_lines", ""],  
								Evaluate[fileFriend[depAssoc[#]]]
					]&/@plotIndices;
		(*fileNames = names of files*)
		display[[;;, 1]] = Column[{disp[[#]], Button["Export plot", 
										Export[StringJoin[fileNames[[#]], ".pdf"], disp[[#]], ImageResolution->200 ]; 
										Export[StringJoin[fileNames[[#]], ".tiff"], disp[[#]], ImageResolution->200 ];
								]}]&/@indexVector
	(*column of graphs with buttons to display*)
	];
	
	(*---------------------------ANOVA MODE----------------------*)
	If[showANOVA, 
		If[Length[sio1]<2 || CountDistinct[Flatten[thisGrid[[;;,;;,2]],1]]<2, getSiO = 0, getSiO=2];
		If[Length[a1]<2 || CountDistinct[Flatten[thisGrid[[;;,;;,3]],1]]<2, getAcetone = 0, getAcetone=3];
		If[Length[v1]<2 || CountDistinct[Flatten[thisGrid[[;;,;;,5]],1]]<2, getAmplitude= 0, getAmplitude=5];
		If[Length[s1]<2 || CountDistinct[Flatten[thisGrid[[;;,;;,7]],1]]<2, getSpeed = 0, getSpeed=7];
		If[Length[p1]<2 || CountDistinct[Flatten[thisGrid[[;;,;;,4]],1]]<2, getParts = 0, getParts=4];
		list = Select[{getSiO, getAcetone, getAmplitude, getSpeed, getParts}, #>0&]; (*list of indices for factors to include in ANOVA*)
		factors = nmAssoc/@list; (*factors to include in ANOVA*)
		If[Length[list]>0,
			jjj = ANOVA[DeleteCases[Flatten[thisGrid,1][[;;, Flatten[{list,{#}}, 1]]], Append[ConstantArray[_, Length[list]], ""]],factors, nmAssoc/@list]&/@plotIndices;
			,
			jjj = ConstantArray[{ANOVA->"", CellMeans->{{{"No factors selected"}}}}, Length[plotIndices]];
		];
		display[[;;,2]] = Column[{
			Style[StringJoin["ANOVA of ", depAssoc[plotIndices[[#]]]], "Section"]
			,
			Row[{SetAccuracy[(ANOVA/.(jjj[[#]])), 4]}]
			,
			Row[{}]
			,
			Row[Grid[#, Frame->All, FrameStyle->Directive[Gray]]&/@SetAccuracy[SplitBy[(CellMeans/.jjj[[#]])[[1]],
				StringJoin[StringTake[StringSplit[ToString[#[[1]]], {"] "}],3]] &], 4]  , "\t"]
		}]&/@Range[Length[plotIndices]]
	];
	If[showPlot&&showANOVA, Grid[display], If[showPlot, Row[display[[;;,1]]], If[showANOVA, Column[display[[;;,2]]], "Select plot and/or ANOVA"]]]
(*Grid[display]*)
,
"Nothing selected"]
}]]

,{{showPlot, True}, {True, False}}
,{{showANOVA, True}, {True, False}}
,Delimiter
(*select subsets based on independent variables*)
,{{sio1, ALLSIOVALS}, Join[{#}&/@ALLSIOVALS, {ALLSIOVALS, {5., 6., 7., 8., 9., 10.},{5., 6., 6.5, 8., 9., 10.}, {5., 6.,7., 9., 10.}, {5., 7., 9.}, {6., 8., 10.}}]}
, {{a1, ALLAVALS}, DeleteCases[Subsets[ALLAVALS], {}]}
, {{v1, ALLVVALS}, DeleteCases[Subsets[ALLVVALS],{}]}
, {{s1, ALLSVALS}, DeleteCases[Subsets[ALLSVALS],{}]}
, {{p1, {120}}, DeleteCases[Subsets[ALLPVALS],{}]}
, Delimiter
(*independent variables to plot*)
,{{x,2, "x axis"}, Normal[nameAssoc]}
,{{c,5, "color"}, Append[Normal[KeyTake[nameAssoc, viableColors]], 0->"None"]}
, Delimiter
(*plot parameters*)
(*, {{colorScheme, "DarkRainbow"}, ColorData["Gradients"][[{4,5,6,8,14, 17, 20,26,32,36,37,40,41,45,51}]]}*)
, {{showMeans, True, "Show means with Err bars"}, {True, False}}
, {{fixedRange, False, "Fix plot range"}, {True, False}}
, {{showLines, False, "Show best fit lines"}, {True, False}}
, {{model, "linear", "Fit"}, {"linear", "quadratic", "^1/2", "^-1"}}
,{{plotMode, 0, "Plot mode"}, {0->"Scatter", 1->"3D", 2->"Contour"}}
, {{io, 0, "Interpolation order"}, {0,1}, ControlType->Setter}
(*,Delimiter
,{{getSiO, 2,"Silica"},{0,2}, ControlType->Checkbox}
,{{getAcetone, 3, "Acetone"}, {0,3}, ControlType->Checkbox}
,{{getAmplitude, 5, "Amplitude"}, {0,5}, ControlType->Checkbox}
,{{getSpeed, 7, "Speed"}, {0,7}, ControlType->Checkbox}*)
,Delimiter
, {{onlyCont, True, "Only continuous lines"}, {True, False}}
, {{minArea, 0, "Minimum cross-sectional area"}}
, {{maxArea, 10, "Maximum cross-sectional area"}}
,{{minwh, 0, "Minimum STDEV height/height"}, 0, 0.6}
,{{maxwh, 0.6, "Maximum STDEV height/height"}, 0.01, 0.6}
(*, ControlPlacement->Left*)
,ContinuousAction->False
]


(* ::Subsubsection::Closed:: *)
(*anova*)


(* ::Text:: *)
(*anova interface*)


anovaInterface:=Manipulate[Module[{sio3, a3, v3, s3, grid2, list, factors, jjj},
sio3 =sio1; a3 = a1; v3 = v1; s3 = s1;
If[Length[sio3]<2, getSiO = 0];
If[Length[a3]<2, getAcetone = 0];
If[Length[v3]<2, getAmplitude= 0];
If[Length[s3]<2, getSpeed = 0];
grid2 =Select[ALLSTATS, MemberQ[sio3,#[[2]]] && MemberQ[ a3, #[[3]]] && MemberQ[ v3, #[[5]]] && MemberQ[ s3, #[[7]]]&];
list = Select[{getSiO, getAcetone, getAmplitude, getSpeed}, #>0&];
factors = If[interactions, Flatten[{nmAssoc/@list, {All}}, 1], nmAssoc/@list];
If[Length[list]>0,
jjj = ANOVA[DeleteCases[grid2[[;;, Flatten[{list,{dependent}}, 1]]], Append[ConstantArray[_, Length[list]], ""]],factors, nmAssoc/@list];
,
jjj = {ANOVA->"", CellMeans->{{{"No factors selected"}}}};
];
anovaGrid = SetAccuracy[(ANOVA/.jjj), 6];
Column[{
Row[{"ANOVA of ", depAssoc[dependent]}]
,
Row[{anovaGrid}]
,
Row[{}]
,
Row[Grid[#, Frame->All, FrameStyle->Directive[Gray]]&/@SetAccuracy[SplitBy[(CellMeans/.jjj)[[1]],
			StringJoin[StringTake[StringSplit[ToString[#[[1]]], {"] "}],3]] &], 5]  , "\t"]
,
Row[{}]
,
If[showGrid, Grid@Prepend[grid2, HEADER], ""]}]]
, {{sio1,ALLSIOVALS}, Append[{#}&/@ALLSIOVALS, ALLSIOVALS]}
, {{a1, ALLAVALS}, Append[{#}&/@ALLAVALS,ALLAVALS]}
, {{v1, ALLVVALS}, Append[{#}&/@ALLVVALS,ALLVVALS]}
, {{s1, ALLSVALS}, Append[{#}&/@ALLSVALS,ALLSVALS]}
, Delimiter
, {{getSiO, 2, "Silica"}, {2->True, 0->False}}
,{{getAcetone, 3, "Acetone"}, {3->True, 0->False}}
,{{getAmplitude, 5, "Amplitude"}, {5->True, 0->False}}
,{{getSpeed, 7, "Speed"}, {7->True, 0->False}}
,{{interactions, False, "Interactions"}, {True, False}}
,{{dependent, 10}, Normal[depAssoc]}
, Delimiter
,{{showGrid,False, "Show Grid"}, {True, False}}]


(* ::Subsubsection::Closed:: *)
(*contours*)


(* ::Text:: *)
(*make a contour plot based on two factor block averages*)


twoFactorContour[x_, index_, V_]:=Module[{y, plot, xlabel, ylabel, depName, dirName, V1},
V1 = ToString[V];
y = {StringSplit[StringSplit[ToString[#], "["][[2]], "]"][[1]], StringSplit[StringSplit[ToString[#], "["][[3]], "]"][[1]], #[[2]]}&/@x; (*y is grid of values we're plotting {ind1, ind2, dependent}*)
xlabel = short2LongAssoc[StringTake[StringSplit[ToString[x[[1]]], "["][[1]], 2;;]]; (*strips the variable names from the table and finds the corresponding string name*)
ylabel = short2LongAssoc[StringSplit[StringSplit[ToString[x[[1]]], " "][[2]], "["][[1]]];
depName = depAssoc[index]; (*dependent string*)
plot = ListContourPlot[ToExpression@y, FrameLabel->{xlabel, ylabel},PlotLegends->Automatic,
	 ImageSize->Medium, FrameStyle->Black, PlotRangePadding->Scaled[0.0],LabelStyle->{FontSize->14,FontColor->Black},  PlotLabel->depName, InterpolationOrder->0];
dirName = StringJoin["E:\\documents\\Epoxy map data\\SiC epoxy map data\\map\\plots\\contours\\", fileFriend@depName];
If[!DirectoryQ[dirName], CreateDirectory[dirName]];
Export[StringJoin[dirName, "\\", fileFriend@depName, "_", xlabel, "_", ylabel, "_v", V1, ".pdf"], plot];
Export[StringJoin[dirName, "\\", fileFriend@depName, "_", xlabel, "_", ylabel, "_v", V1, ".bmp"], plot];	
plot
]


fileFriend[s_]:=StringReplace[ToString[s], {"/"->" over ", "."->","}]


(* ::Text:: *)
(*get all contour plots for one dependent variable*)


getContours[index_]:=Module[{gridVall, list, factors, jjj, grids, x, minsize},
gridVall  = ALLSTATSMORE[[2;;]];
	(*get all stats*)
list = {2,3,5,7}; 
	(*list of ink parameters*)
factors = Flatten[{nmAssoc/@list, {All}}, 1]; 
	(*variable names*)
jjj = ANOVA[DeleteCases[gridVall[[;;, Flatten[{list,{index}}, 1]]], {_,_,_,_, ""}],factors, nmAssoc/@list]; 
	(*anova results with interactions*)
x = {Length[ALLAVALS],Length[ALLSVALS],Length[ALLSIOVALS],Length[ALLVVALS]};
minsize = MinMax@Reap[Do[Do[Sow[x[[i]]*x[[j]]], {i, 1, j-1}], {j, 1,4}]][[2,1]];
grids = Select[SplitBy[(CellMeans/.jjj)[[1]],
			StringJoin[StringTake[StringSplit[ToString[#[[1]]], {"] "}],3]] &] , minsize[[2]]>=Length[#]>=minsize[[1]]&];
	(*take the cell means for two-factor interactions*)
Print@Row[Prepend[twoFactorContour[#, index, "all"]&/@grids, Style["all V", Large]], "\t\t"];
	(*print plots for all two-factor interactions*)
Do[
	gridVall =Select[ALLSTATSMORE, MemberQ[ {thisV}, #[[5]]]&];
		(*get just the stats for this voltage*)
	list = {2,3,7};
		(*get sio, acetone, speed*)
	factors = Flatten[{nmAssoc/@list, {All}}, 1];
		(*variable names*)
	jjj = ANOVA[DeleteCases[gridVall[[;;, Flatten[{list,{index}}, 1]]], {_,_,_,""}],factors, nmAssoc/@list];
		(*anova stats*)
	grids = Select[SplitBy[(CellMeans/.jjj)[[1]],
			StringJoin[StringTake[StringSplit[ToString[#[[1]]], {"] "}],3]] &] , minsize[[2]]>=Length[#]>=minsize[[1]]&];
		(*only the two-factor interactions*)
	Print@Row[Prepend[twoFactorContour[#, index, thisV]&/@grids, Style[StringJoin["V=", ToString[thisV]], Large]], "\t\t"];
,{thisV, ALLVVALS}];
]


(* ::Subsubsection:: *)
(*covariance plots*)


(* ::Text:: *)
(*covariance plots allow you to plot dependent variables against each other to see if there is a correlation between dependent variables*)


covPlot[thisGrid_, d1_, d2_, c_,s_, cf_, ps_, showLines_]:=Module[{g1, gminmax, f, line, reportTable, shapes, svals, plots, lineStyles, g2, lines, numVals},
(*thisGrid is a grid of data*)
(*d1 is the x-axis dependent variable*)
(*d2 is the y-axis dependent variable*)
(*c is the color variable*)
(*cf is the color function (e.g. "RainbowColors")*)
	If[s>0, g1 = Cases[thisGrid[[;;, {d1, d2, c, s}]],  {_?NumberQ, _?NumberQ, _?NumberQ, _?NumberQ}];,g1 = Cases[thisGrid[[;;, {d1, d2, c}]],  {_?NumberQ, _?NumberQ, _?NumberQ}];]  (*only plot points that have values*)
	If[c==25, f=Log, f=1.*#&]; (*use a log scale for viscosity*)
	gminmax = MinMax[f/@g1[[;;,3]]]; (*min and max color value*)
	g1[[;;,3]] = Rescale[f/@g1[[;;,3]]]; (*rescale the colors onto a scale from 0 to 1*)
	line = LinearModelFit[g1[[;;,{1,2}]], {1, xx}, xx];
	reportTable = {{"Correlation", Correlation[g1[[;;,1]], g1[[;;,2]]]}, 
					{"Linear Fit", Normal[line]}, 
					{"Linear fit R^2", line["RSquared"]}, 
					{"Linear fit p", line["ANOVATablePValues"][[1]]}};
		
	If[s>0, 
		svals = DeleteDuplicates[g1[[;;,4]]];
		numVals = Length[svals];
		g2 = GatherBy[g1, Last];
		lines = LinearModelFit[#[[;;,{1,2}]], {1,xx}, xx]&/@g2;
		reportTable = Prepend[Transpose[Join[Transpose[reportTable], {Correlation[g2[[#, ;;,1]], g2[[#, ;;, 2]]], Normal[lines[[#]]], lines[[#]]["RSquared"], lines[[#]]["ANOVATablePValues"][[1]]}&/@Range[numVals]]], Join[{"", "All"}, svals]];
		reportTable[[2;;,2;;]] = SetAccuracy[reportTable[[2;;,2;;]],4];

		shapes = {"\[Times]","\[FilledCircle]", "\[EmptyCircle]",  "\[FilledUpTriangle]","\[FivePointedStar]", "\[FilledSquare]",  "\[FilledDiamond]"};
		lineStyles = {Black, {Black, Dotted}, {Black, Dashed}, {Black, DotDashed}, {Black, Thin}};

		Column[{
				Row[{
					Show[
						Graphics[Text[Style[shapes[[Position[svals,#[[4]]][[1,1]]]],ps, Bold,ColorData[cf][#[[3]]]],{#[[1]], #[[2]]} ]&/@g1, 
							Frame->True, FrameStyle->Black, AspectRatio->1,ImageSize->400,LabelStyle->20, ImageMargins->10, FrameLabel->{depAssoc[d1], depAssoc[d2]}, PlotRangePadding->Scaled[0.05]]
						,
						If[showLines, Plot[(lines[[#]])[xx], {xx, Min[g1[[;;,1]]], Max[g1[[;;,1]]]}, PlotStyle->lineStyles[[#]]]&/@Range[numVals], Graphics[]]
						,
						If[d1==33 && d2==53 || d1==53 && d2==33, Plot[xx, {xx, Min[g1[[;;,1]]], Max[g1[[;;,1]]]}, PlotStyle->{Black, Thickness[0.015]}], Graphics[]]
					]
					,
					Column[{
						BarLegend[{cf, gminmax}, LegendLabel->If[c==18, "Log ", ""]<>HEADER[[c]], LabelStyle->{20, Black}, LegendMarkerSize->200],
						If[showLines, 
							LineLegend[lineStyles[[1;;numVals]], Row/@Transpose[{Style[#, Bold, ps]&/@shapes[[1;;numVals]], ConstantArray["   ", numVals],Style[ToString[#]<>" "<>unitAssoc[s], 20]&/@svals}]]
							,
							Column[Row/@Transpose[{Style[#, Bold, ps]&/@shapes[[1;;numVals]], ConstantArray["   ", numVals],Style[ToString[#]<>" "<>unitAssoc[s], 20, FontFamily->"Arial"]&/@svals}]]
						]
					}]
				}]
				,
				Grid[reportTable, Alignment->Right, Dividers->{All, False}, ItemStyle->Directive[FontSize->12, FontFamily->"Arial"]]
		}, Center]
		
		,

		Column[{
			Row[{
				Show[
					Graphics[Text[Style["\[FilledCircle]",ps, Bold,ColorData[cf][#[[3]]]],{#[[1]], #[[2]]} ]&/@g1, 
						Frame->True, FrameStyle->Black, AspectRatio->1,ImageSize->400,LabelStyle->20, ImageMargins->10, FrameLabel->{depAssoc[d1], depAssoc[d2]}, PlotRangePadding->Scaled[0.05]]
					,
					If[showLines, Plot[line[xx], {xx, Min[g1[[;;,1]]], Max[g1[[;;,1]]]}, PlotStyle->Black], Graphics[]]
					,
					If[d1==33 && d2==53 || d1==53 && d2==33, Plot[xx, {xx, Min[g1[[;;,1]]], Max[g1[[;;,1]]]}, PlotStyle->{Black, Thickness[0.015]}], Graphics[]]
				]
				,
				Column[{BarLegend[{cf, gminmax}, LegendLabel->If[c==18, "Log ", ""]<>HEADER[[c]], LabelStyle->{20, Black}, LegendMarkerSize->200]}]
			}],
				Grid[reportTable, Alignment->Right, Dividers->{All, False}, ItemStyle->Directive[FontSize->12, FontFamily->"Arial"]]
		}, Center]
	]
]


covInterface:=Manipulate[Module[{thisGrid, xPlot, sio2, a2, v2, s2,runName, p2},
	{thisGrid, runName, sio2, a2, v2, s2, p2}=PIgetGrid[sio1, a1, v1, s1,p1,  2,2, onlyCont, maxArea, minArea, 0, 1000];
	plot = covPlot[thisGrid[[1]], var1, var2, color,shape, colorScheme, pointSize, showLines];
	fileName = StringJoin[BASEDIR, "\\plots\\corr plots\\", 
								runName, fileFriend[var1], "_",fileFriend[var2],"_", fileFriend[color],"_",
								fileFriend[maxArea],"_", fileFriend[minArea],"_",If[onlyCont, "cont",""]];
	Column[{plot, Button["Export", Export[StringJoin[fileName, ".pdf"], plot];Export[StringJoin[fileName, ".bmp"], plot];(* Print[fileName];*)]}]

]
, {{var1,33, "Var 1"}, (#->ToString[depAssoc[#]])&/@Keys[depAssoc]}
, {{var2,53, "Var 2"}, (#->ToString[depAssoc[#]])&/@Keys[depAssoc]}
, {{color,2, "Color"}, Join[(#->ToString[depAssoc[#]])&/@Keys[depAssoc], (#->ToString[nameAssoc[#]])&/@Keys[nameAssoc]]}
, {{shape,5, "Shape"}, Append[(#->nameAssoc[#])&/@{3, 4, 5, 7}, 0->"None"]}
, {{sio1, ALLSIOVALS}, Join[{#}&/@ALLSIOVALS, {ALLSIOVALS, {5., 6., 7., 8., 9., 10.}}]}
, {{a1, ALLAVALS}, DeleteCases[Subsets[ALLAVALS], {}]}
, {{v1, ALLVVALS}, DeleteCases[Subsets[ALLVVALS],{}]}
, {{s1, ALLSVALS}, DeleteCases[Subsets[ALLSVALS],{}]}
, {{p1, {120}}, DeleteCases[Subsets[ALLPVALS],{}]}
, {{onlyCont, True, "Only continuous lines"}, {True, False}}
, {{minArea, 0, "Minimum cross-sectional area"}, 0, 1}
, {{maxArea, 10, "Maximum cross-sectional area"}, 1, 10}
, {{colorScheme, "DarkRainbow"}, ColorData["Gradients"][[{8, 17, 20,26,32,36,37,40,41,45,51}]]}
, {{pointSize, 15}, 5, 30}
, {{showLines, False, "Show lines"}, {True, False}}
, ContinuousAction->False
]
