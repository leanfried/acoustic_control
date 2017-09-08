(* ::Package:: *)

(* ::Title:: *)
(*Keyence functions*)


(* ::Section:: *)
(*Initialization*)


(* ::Text:: *)
(*Last updated 06/01/16 by Leanne Friedrich*)


<<Functions_master.wl;
kFiles = FileNames[{"K*1.csv","K*2.csv", "K*3.csv"}, BASEDIR, Infinity];
	(*keyence files*)


(* ::Section:: *)
(*Analysis*)


(* ::Text:: *)
(*getWidthHeight gets the surface profile cross section and its width and height*)
(*extra comments are for processing each cross-section to analyze variations throughout the line - this takes a long time and I didn't end up using it, so it's commented out here*)


getWidthHeight[file_]:=Module[{fileArray, l, wh, p, fa, fawh, ks},
Catch[
(*x = Import[file];
y = ConstantArray[0, {801, 3}];
Do[
	l = clean[x[[;;, i]]];
	y[[i-399, ;;]] = getWH[l];
,{i, 400, 1200}]

(*Export[StringJoin[StringSplit[file, "."][[1]], "_wh.txt"], "csv"];*)
Mean[y]*)
ks = getKS[file];
	fileArray = Transpose[Import[file, "csv"][[10;;1190, 400;;1200]]];  (*import the file, crop the outer 10 pixels, and take the middle half of the line*)
(*	fa = (*clean/@fileArray;*)fileArray;
	fawh = getStats/@fa;
	fawh = DeleteCases[fawh, #[[6]]\[Equal]0&];*)
	l = Mean@fileArray; (*average the cross sections across the line*)
l = clean[l]; (*remove glare*)
(*Print[ListPlot[l]];*)
(*wh = getWH[l];*)
{p,wh} = getStats[l]; (*get stats about the profile*)
wh = wh/ks;
wh[[1]] = wh[[1]]*ks;
(*p = ListPlot[l, Epilog->{Red, Line[{#, wh[[1]]/2}&/@wh[[2;;3]]]}, PlotStyle->Black, PlotRange->All];*)
(*Print[Column[{p, Grid[{{"Height", "Left","Right", "FWHM"}, wh}]}]];*)
(*Export[StringReplace[file, ".csv"->"_profile.tiff"], p]; (*export an image of the profile*)
(*Export[StringReplace[file, ".csv"->"_allprofs.csv"], fawh];
Export[StringReplace[file, ".csv"->"_allprofsStats_f.csv"], 1.*{Mean[fawh], StandardDeviation[fawh]}];*)*)
Export[StringReplace[file, ".csv"-> "_stat_f.txt"], Transpose[{{"height (\[Micro]m)", "full width (\[Micro]m)", "full area (\[Micro]m^2)",
		 "full width/height", "full area/height (\[Micro]m)", "full width at half max (\[Micro]m)", 
		"full area at half max (\[Micro]m^2)", "full width at half max/height", "full area at half max/height (\[Micro]m)", "stdev (\[Micro]m)", "stdev/height"}, wh}], "csv"]; (*export stats of the profile*)
(*Export[StringReplace[file, ".csv"-> "_profile.txt"], l, "csv"]; (*export the profile*)*)
Return[0];]
]


getfileArray[file_]:=Mean@Transpose[Import[file, "csv"][[10;;1190, 400;;1200]]]


(* ::Text:: *)
(*sometimes if the anti-glare spray is not applied correctly, glare will appear on the sides of the line, showing up as deep grooves in the surface profile. clean figures out if there is large glare and replaces it with a straight line between the nearest non-glare regions*)


clean[l_]:=Module[{x, func, minn, leftCut, m, l3, rightCut,l2, minLoc},Catch[
	l2 = getderivatives[l, 20]; (*numerically determine the derivative of the profile, using an approximation radius of 20 pixels*)
	l3 = l;
(*Print[Min[l3]-Mean[Join[l3[[1;;5]], l3[[-5;;-1]]]]];
Print[Mean[Join[l[[1;;5]], l[[-5;;-1]]]]];*)
	If[Min[l2[[3]]]<-50 || Min[l]-Mean[Join[l[[1;;5]], l[[-5;;-1]]]]<-10,  
		(*if there is a steep decrease or if there is a point well below the background, fix glare*)
		minLoc = Ordering[l2[[3]],1]; (*location of the dip*)
		func = Interpolation[l]; (*get a function that interpolates the surface profile*)
		minn = Floor[x/.(Minimize[{func[x], Max[minLoc-40,2]<x<Min[minLoc+40, Length[l]]}, x][[2]])]; (*find the middle of the glare*)
			leftCut = Floor[x/.(Maximize[{func[x], minn-50<x<minn}, x][[2]])]; (*find the left edge of the glare*)
(*			Print[minn,"  ",  leftCut, "   ", Show[Plot[func[x], {x, 1, Length[l]}], ListPlot[l], ImageSize->Medium], Show[ListPlot[l2[[3]]],PlotRange->All, ImageSize->Medium]];
		*)	rightCut = 2*minn -leftCut; (*find the right edge of the glare*)
			m = (l[[rightCut]] - l[[leftCut]])/(rightCut - leftCut); (*find the slope between the left edge and right edge*)
			l3[[leftCut;;rightCut-1]] = ConstantArray[l3[[leftCut]], rightCut-leftCut] + m*Range[rightCut-leftCut]; (*replace the glare region with a straight line between the left and right edge*)
(*			Print[ListPlot[l3]];*)
	];
	l3 = l3 - EstimatedBackground[l3, 100]; (*subtract background*)
	Return[l3];
]]


(* ::Text:: *)
(*get the first derivative of a profile*)


getderivatives[l_, m_]:=Module[{k, l2},
l2 = ConstantArray[0, {4,Length[l]}];
k = 2*m-1;
l2[[1]] = Range[Length[l]]; (*row 1 is the pixel number (distance across the line) - Keyence analysis doesn't use this, but you can use it to help diagnose issues*)
l2[[2]] = l; (*row 2 is the profile*)
l2[[3, m;;-m]] = l[[k;;-1]] - l[[1;;-k]]; (*row 3 is the first derivative of the profile*)
(*l2[[4, m;;-m]] = l2[[3, k;;-1]] - l2[[3,1;;-k]];*)
l2
]


getStats[l_]:=Module[{maxx, baseline, minxleft, minxright, fw, fullArea, height, halfmax, left, right, fwhm, halfArea,p, fullSD},
	maxx = Ordering[l, -1][[1]]; (*x location of the peak*)
	baseline = Min[l]; (*height of substrate*)
	minxleft = Ordering[l[[1;;maxx]], 1][[1]]; (*location of left edge of line*)
	minxright = Ordering[l[[maxx;;-1]], 1][[1]]+maxx-1; (*location of right edge of line*)
	fw = minxright-minxleft; (*full width*)
	fullArea = Total[l[[minxleft;;minxright]]]-(minxright-minxleft)*baseline;
	height = Max[l]-Min[l];
	halfmax = height*0.5+Min[l];
	If[maxx>21 && maxx<(Length[l]-21), (*if your maximum height is at the edge of the image, something went wrong*)
		left = Max[Ordering[Abs[#-halfmax]&/@(l[[1;;maxx]]), {1,20}]]; (*left point of halfmax*)
		right = Min[Ordering[Abs[#-halfmax]&/@(l[[maxx;;-1]]), {1,20}]]+maxx-1; (*right point of half max*)
		fwhm = right-left;  (*full width at half max*)
		halfArea = Total[l[[left;;right]]]-(right-left)*baseline;
		,
		fwhm = 0;
		halfArea = 0;];
	p=ListPlot[l, Epilog->{Red, Line[{{left, halfmax}, {right, halfmax}}], Line[{{minxleft, baseline}, {minxright, baseline}}]}, 
					ImageSize->500,PlotStyle->Black, PlotRange->{{0, Length[l]}, {0, Max[l]}}]; (*graphical representation of stats*)
	fullSD = getSD[Transpose[{Range[minxright-minxleft+1], l[[minxleft;;minxright]]}]][[2]]; (*mean and standard deviation of the profile within the full width*)
	{p, {height, fw, fullArea, fw/height, fullArea/height, fwhm, halfArea, fwhm/height, halfArea/height, fullSD, fullSD/height} } (*return stats*)
	
]


fixstatf[f_]:=Module[{x, x2},
x = Import[f, "csv"][[;;,1]];
If[Length[x]==10, 
x2 = {{"height (micron)", "full width (px)", "full area (px*micron)", "full width/height (px/micron)", "full area/height (px)", "full width at half max (px)", 
		"full area at half max (px*micron)", "full width at half max/height (px/micron)", "full area at half max/height (px)", "stdev (px)", "stdev/height (px/micron)"}, 
		Append[x, x[[10]]/x[[1]]]};
Export[f, x2, "csv"]];]


fixstatf2[f_]:=Module[{x, ks},
x = Import[f, "csv"];
If[x[[1,1]] =="height (micron)",
ks = getKS[f];
x = x/ks;
x[[1,2]] = x[[1,2]]*ks;
x[[;;, 1]] = {"height (\[Micro]m)", "full width (\[Micro]m)", "full area (\[Micro]m^2)", "full width/height", "full area/height (\[Micro]m)", "full width at half max (\[Micro]m)", 
		"full area at half max (\[Micro]m^2)", "full width at half max/height", "full area at half max/height (\[Micro]m)", "stdev (\[Micro]m)", "stdev/height"};
Export[f, x, "csv"]]]


fixstatf2/@FileNames["*stat_f*", "*", Infinity];


(* ::Section:: *)
(*Troubleshooting*)


(* ::Text:: *)
(*compare finds folders where there is a discrepancy between the Keyence and ILM data and prints out their stats*)


compare[d_]:=Module[{ilmFiles, kFiles, ilm, k},
ilmFiles = FileNames["*ILM*", d];
kFiles = FileNames["*K*", d];
If[Length[ilmFiles]>0,
If[ Length[kFiles]>0, 
ilm= StringSplit[ilmFiles[[1]], "_"][[-2]];
k = StringSplit[kFiles[[1]], "_"][[-2]];
If[!StringMatchQ[ilm, k], Print[{d, ilm, k }]];
,
Print[d, "    ILM file but no k file"];
]
,
If[Length[kFiles]>0,
Print[d, "    K file but no ILM file"];]]
]


(*getWH[l_]:=Module[{halfmax, y, maxx, left, right},
	halfmax = (Max[l]-Min[l])*0.5+Min[l];
	maxx = Ordering[l, -1][[1]];
	If[maxx>21 && maxx<(Length[l]-21),
	left = Max[Ordering[Abs[#-halfmax]&/@(l[[1;;maxx]]), {1,20}]];
	right = Min[Ordering[Abs[#-halfmax]&/@(l[[maxx;;-1]]), {1,20}]]+maxx;
	{Max[l]-Min[l], left, right, right-left} 
		,
	ConstantArray[0, 4]]
]*)


(*KInterface:=Manipulate[Module[{fileArray, l, lmax, a, sigma, x, v2, v},
	fileArray = Import[kFiles[[i]], "csv"];
	l = Mean@Transpose@fileArray[[10;;1190, 400;;1200]];
	lmax = Max[l];
	v = NonlinearModelFit[{#, l[[#]]/Total[l]}&/@Range[Length[l]], 1/Sqrt[2*Pi*sigma^2]*Exp[-(x - a)^2/(2*sigma^2)], {{sigma, 300}, {a, 600}},  x];
	v2 = v[x]*Total[l];
	Column[{Show[{ListPlot[l, ImageSize->Medium],ListPlot[clean[l], PlotStyle->Blue], Plot[v2, {x, 0, 1200}, PlotStyle->Red, ImageSize->Medium, PlotRange->All]}]}]
], {i, 1, Length[kFiles]},{{aGuess, 600}, 1, 1200}, {{sigmaGuess, 300}, 1, 10000}]*)


(*redo[dir_]:=Module[{},
files = FileNames[{"K*stats.txt"}, dir];
Do[
If[Length[Import[files[[i]], "csv"]]<4,
f2 = FileNames[StringReplace[files[[i]], {"_stats.txt"\[Rule]".csv"}]];
getWidthHeight[f2[[1]]];
Print[f2];
]
,{i, Length[files]}];
]*)
