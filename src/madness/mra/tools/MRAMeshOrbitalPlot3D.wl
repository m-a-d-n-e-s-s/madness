(* ::Package:: *)

BeginPackage["MRAMeshOrbitalPlot3D`"]


(* ::Text:: *)
(*{*)
(* {Description: MRAMeshOrbitalPlot3D[fp,opts] returns a plot of the orbital stored via madness::plot_cubefile and madness::print_tree_jsonfile in files fp . cube and fp . tree . json, respectively . MRAMeshOrbitalPlot3D accepts all options recognized by Graphics3D and ListContourPlot3D functions, as well as the following additional options:, \[SpanFromLeft], \[SpanFromLeft]},*)
(* {Option, Default, Description},*)
(* {Zoom, 1, can be used to produce a plot in a zoomed-in section of the simulation cell . This does not need to match the zoom value given to plot_cubefile (that only affects the resolution/extent of the Gaussian Cube mesh)},*)
(* {MRAMeshCuboidDirectives, {EdgeForm[Thick]}, Specifies how the Cuboid objects comprising the MRA mesh are drawn . All Graphics3D directives that are applicable to Cuboid (except Opacity) can be specified . },*)
(* {MaxLevel, Infinity, Controls the highest refinement level of displayed mesh elements .},*)
(* {MinLevel, 0, Controls the lowest refinement level of displayed mesh elements .}*)
(*}*)


MRAMeshOrbitalPlot3D::usage="MRAMeshOrbitalPlot3D[fp,opts] returns a plot of the orbital stored via madness::plot_cubefile and madness::print_json_treefile in files fp.cube and fp.tree.json, respectively. MRAMeshOrbitalPlot3D accepts all options recognized by Graphics3D and ListContourPlot3D functions, as well as the following additional options: Zoom, MRAMeshCuboidDirectives, MinLevel, and MaxLevel."


Begin["`Private`"]


ReadTree[fileName_] :=
    Module[{jsonData, cellData, treeData, boxCoords, nodesData},
        jsonData = Import[fileName];
        cellData = "cell" /. jsonData;
        treeData = "tree" /. jsonData;
        boxCoords = {};
        Do[
            nodesData = "nodes" /. (ToString[n] /. treeData);
            Do[AppendTo[boxCoords, Prepend[Interpreter["Integer"] /@ 
                StringSplit[Part[nodesData[[i]], 1], {"[", "]", ","}], n]], {i, Length[
                nodesData]}];
            ,
            {n, 0, Length[treeData] - 1}
        ];
        Return[{cellData, boxCoords}]
    ];

BoxToXYZCoords[cell_, box_] :=
    Module[{S},
        S = Table[cell[[xyz, 2]] - cell[[xyz, 1]], {xyz, 3}];
        Return[{{cell[[1, 1]] + S[[1]] box[[2]] / (2 ^ box[[1]]), cell
            [[2, 1]] + S[[2]] box[[3]] / (2 ^ box[[1]]), cell[[3, 1]] + S[[3]] box
            [[4]] / (2 ^ box[[1]])}, {cell[[1, 1]] + S[[1]] (box[[2]] + 1) / (2 ^
             box[[1]]), cell[[2, 1]] + S[[2]] (box[[3]] + 1) / (2 ^ box[[1]]), cell
            [[3, 1]] + S[[3]] (box[[4]] + 1) / (2 ^ box[[1]])}}]
    ];

BoxToGraphics[cell_, box_, nmax_, omax_(* opacity value of the smallest
     boxes *), shadeexp_(* opacity of box at level n-1 is this times smaller
     than than of box at level n *),cuboidDirectives_List] :=
    Module[{},
        Return[Join[(*N.B. MUST BE FIRST to apply to Cuboid*)cuboidDirectives,{Opacity[omax / shadeexp ^ (nmax - box[[1]])], Cuboid 
            @@ BoxToXYZCoords[cell, box],EdgeForm[Thick]}]]
    ];

(*
MaxLevel,MinLevel: show boxes with resolution level [nmin,nmax]
Zoom: limit PlotRange to cell/Zoom
*)

Protect[MaxLevel,MinLevel,Zoom,MRAMeshCuboidDirectives];
Options[MRAMeshOrbitalPlot3D] = {MaxLevel -> Infinity, MinLevel -> 0, Zoom
     -> 1,MRAMeshCuboidDirectives->{EdgeForm[Thick]}};

(* plots orbital and its mesh superimposed *)

MRAMeshOrbitalPlot3D[filePrefix_, opt : OptionsPattern[{MRAMeshOrbitalPlot3D,
     Graphics3D, ListContourPlot3D, Show}]] :=
    Module[
        {simulationCell, boxTreeCoords, bohr2angstrom, selectedBoxes,
             boxGraphics, meshPlot, orbitalPlot, nmax, nmin, zoom, omax, shadeexp,
             actualNMax,cuboidDirectives, plotRange}
        ,
        (* process options *)
        nmax = OptionValue[MaxLevel];
        nmin = OptionValue[MinLevel];
        zoom = OptionValue[Zoom];
        cuboidDirectives=OptionValue[MRAMeshCuboidDirectives];
        
        (* these are hardwired since they are not needed for most users
             *)
        omax = 0;(* opacity value of the smallest boxes *)
        shadeexp = 1.9;(* opacity of box at level n-1 is this times smaller
             than than of box at level n *)
             
        {simulationCell, boxTreeCoords} = ReadTree[filePrefix <> ".tree.json"
            ]; (* the deduced hardwired value used by Wolfram when importing Gaussian
             Cube file *) bohr2angstrom = 0.529177249;
        simulationCell *= bohr2angstrom;
        (* override default PlotRange *)
        plotRange =
            If[OptionValue[PlotRange] === All,
                simulationCell / zoom
                ,
                OptionValue[PlotRange]
            ];
        selectedBoxes = Select[boxTreeCoords, #[[1]] <= nmax && #[[1]]
             >= nmin&];
        actualNMax = MaximalBy[selectedBoxes, #[[1]]&][[1, 1]];
        boxGraphics = Map[BoxToGraphics[simulationCell, #, actualNMax,
             omax, shadeexp, cuboidDirectives]&, selectedBoxes];
        meshPlot = Graphics3D[boxGraphics, Boxed -> False, Evaluate @
             FilterRules[{opt}, Options[Graphics3D]]];
        orbitalPlot = Import[filePrefix <> ".cube", "Graphics3D", Boxed
             -> False, Evaluate @ FilterRules[{opt}, Options[ListContourPlot3D]]]
            ;
        Return[Show[{orbitalPlot, meshPlot}, Evaluate @ FilterRules[{
            opt, PlotRange -> plotRange}, Options[Graphics3D]]]];
    ];



End[]


EndPackage[]
