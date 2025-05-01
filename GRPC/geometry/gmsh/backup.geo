// Glass RPC prototype

SetFactory("OpenCASCADE");

//Mesh.Algorithm = 8;
//Mesh.MeshSizeMax = 0.1;
//Mesh.MeshSizeFromPoints = 0;
//Mesh.MeshSizeExtendFromBoundary = 1;
//Mesh.MeshSizeFromCurvature = 0;
//Mesh.Smoothing = 100;

//Mesh.OptimizeNetgen = 1;
//Mesh.Algorithm = 7;
//Mesh.Algorithm3D = 4;
//Mesh.Algorithm3D = 10;
//Mesh.MeshSizeFromPoints = 1;
//Mesh.MeshSizeExtendFromBoundary = 1;
//Mesh.MeshSizeMax = 0.001;
Mesh.ElementOrder = 2;
//Mesh.SecondOrderLinear = 1;
Mesh.SurfaceFaces = 1;
Mesh.SurfaceEdges = 1;
Mesh.VolumeFaces = 1;
Mesh.VolumeEdges = 1;
Geometry.Points = 0;

Geometry.OCCTargetUnit = "M";

frame() = ShapeFromFile("frame.step");

bbox() = BoundingBox Volume{frame()};
xmin = bbox(0);
ymin = bbox(1);
zmin = bbox(2);
xmax = bbox(3);
ymax = bbox(4);
zmax = bbox(5);

Printf("Import frame with dimensions (meter) '%g','%g' '%g','%g' '%g','%g'", xmin, xmax, ymin, ymax, zmin, zmax);

nsp = 2; // number of spacers per axis
ngap = 4;
glelx = 0.20; // 20cm
glely = 0.20; // 20cm
glelw = 0.0011; // 1.1mm
glgapx = 0.15; // 15cm
glgapy = 0.15; // 15cm
glgapw = 0.0007; // 0.7mm
coatx = 0.15; // 15cm
coaty = 0.15; // 15cm
coatw = 0.00005; // 0.05mm
gap = 0.00024; // 0.24mm
sptapeh = 0.00007; // 0.07mm
spr = 0.002; // 2mm
sph = 0.00024; // 0.24mm

gwidth = ngap*glgapw+(ngap+1)*gap;
gwidthc = ngap*glgapw+(ngap+1)*sph;
Printf("Gap width is '%g' '%g'", gwidth, gwidthc);

Rotate {{1,0,0}, {0,0,0}, -Pi/2} { Volume{frame}; }
Translate {0.5*glelx, 0.5*glely, coatw+glelw+gwidth} { Volume{frame}; }

det = newv;
Box(det) = {0,0,glelw+coatw, glelx,glely,gwidth};
gasv = newv;
BooleanDifference(gasv) = {Volume{det}; Delete;}{Volume{frame};};

//coatings[0] = newv;
//Box(coatings[0]) = {(glelx-coatx)*0.5,(glely-coaty)*0.5,0,coatx,coaty,coatw};
//coatings += Translate {0, 0, glelw+coatw+ngap*(glgapw+gap)+gap+glelw+coatw} { Duplicata{ Volume{coatings[0]}; }};

plane = newv;
Box(plane) = {0,0,coatw, glgapx,glgapy,glgapw};

//electrodet = newv;
//Box(electrodet) = {0,0,coatw, glelx,glely,glelw};
//electrodeb = Translate {0, 0, glelw+coatw+ngap*(glgapw+gap)+gap} { Duplicata{ Volume{electrodet}; }};

planes = {};
For i In {0: ngap-1}
  planes += Translate {(glelx-glgapx)/2., (glely-glgapy)/2., i*(glgapw+gap)+glelw+gap+coatw} { Duplicata{ Volume{plane}; }};
EndFor

spacer = newv;
Cylinder(spacer) = {0, 0, 0, 0, 0, sph, spr};

spacers = {};

lspx = glelx/nsp;
lspy = glely/nsp;

For l In {1: ngap+1}
  For i In {1: nsp}
    For j In {1: nsp}
      spacers += Translate {lspx*i-lspx/2, lspy*j-lspy/2, glelw+coatw+(glgapw+gap)*(l-1)} { Duplicata{ Volume{spacer}; }};
    EndFor
  EndFor
EndFor

Recursive Delete{ Volume{spacer}; }
Recursive Delete{ Volume{plane}; }

For i In {0: #planes[]-1}
  gasv0 = gasv; gasv = newv;
  BooleanDifference(gasv) = {Volume{gasv0}; Delete;}{Volume{planes[i]};};
EndFor
For i In {0: #spacers[]-1}
  gasv0 = gasv; gasv = newv;  
  BooleanDifference(gasv) = {Volume{gasv0}; Delete;}{Volume{spacers[i]};};
EndFor

/*
bb = newv;
Box(bb) = {0,0,0, glelx,lspy/2,gwidth+2*glelw+2*coatw*(1.5)};

frame_cut = newv;
BooleanDifference(frame_cut) = {Volume{frame};}{Volume{bb};};
gasv_cut = newv;
BooleanDifference(gasv_cut) = {Volume{gasv};}{Volume{bb};};
electrodet_cut = newv;
BooleanDifference(electrodet_cut) = {Volume{electrodet};}{Volume{bb};};
electrodeb_cut = newv;
BooleanDifference(electrodeb_cut) = {Volume{electrodeb};}{Volume{bb};};

planes_cut = {};
For i In {0: #planes[]-1}
  cut = newv;
  BooleanDifference(cut) = {Volume{planes[i]};}{Volume{bb};};
  planes_cut += cut;
EndFor

coatings_cut = {};
For i In {0: #coatings[]-1}
  cut = newv;
  BooleanDifference(cut) = {Volume{coatings[i]};}{Volume{bb};};
  coatings_cut += cut;
EndFor

spacers_mylar_cut = {};
For i In {0: #spacers_mylar[]-1}
  cut() = BooleanDifference{Volume{spacers_mylar[i]};}{Volume{bb};};
  spacers_mylar_cut += cut();
EndFor

spacers_tapet_cut = {};
For i In {0: #spacers_tapet[]-1}
  cut() = BooleanDifference{Volume{spacers_tapet[i]};}{Volume{bb};};
  spacers_tapet_cut += cut();
EndFor

spacers_tapeb_cut = {};
For i In {0: #spacers_tapeb[]-1}
  cut() = BooleanDifference{Volume{spacers_tapeb[i]};}{Volume{bb};};
  spacers_tapeb_cut += cut();
EndFor

Recursive Delete{ Volume{bb}; }*/

//Coherence;
Mesh 3;

//Physical Volume("electrodet", 1) = {electrodet};
//Physical Volume("electrodeb", 2) = {electrodeb};
//Physical Volume("coatings", 3) = {coatings[]};
Physical Volume("planes", 4) = {planes[]};
//Physical Volume("spacers_mylar", 5) = {spacers_mylar[]};
Physical Volume("frame", 6) = {frame};
//Physical Volume("gasv", 7) = {gasv};
/*
Physical Volume("electrodet_cut", 11) = {electrodet_cut};
Physical Volume("electrodeb_cut", 12) = {electrodeb_cut};
Physical Volume("coatings_cut", 13) = {coatings_cut[]};
Physical Volume("planes_cut", 14) = {planes_cut[]};
Physical Volume("spacers_mylar_cut", 15) = {spacers_mylar_cut[]};
Physical Volume("spacers_tapet_cut", 16) = {spacers_tapet_cut[]};
Physical Volume("spacers_tapeb_cut", 17) = {spacers_tapeb_cut[]};
Physical Volume("frame_cut", 18) = {frame_cut};
Physical Volume("gasv_cut", 19) = {gasv_cut};
*/
//Color Snow { Surface{Boundary{Volume{electrodet, electrodeb};}}; }
//Color Snow { Volume{electrodet, electrodeb}; }
//Color Snow { Surface{Boundary{Volume{electrodet_cut, electrodeb_cut};}}; }
//Color Snow { Volume{electrodet_cut, electrodeb_cut}; }

//Color DarkGrey { Surface{Boundary{Volume{coatings[]};}}; }
//Color DarkGrey { Volume{coatings[]}; }
//Color DarkGrey { Surface{Boundary{Volume{coatings_cut[]};}}; }
//Color DarkGrey { Volume{coatings_cut[]}; }

//Color WhiteSmoke { Surface{Boundary{Volume{planes[]};}}; }
//Color WhiteSmoke { Volume{planes[]}; }
//Color WhiteSmoke { Surface{Boundary{Volume{planes_cut[]};}}; }
//Color WhiteSmoke { Volume{planes_cut[]}; }

//Color AntiqueWhite { Surface{Boundary{Volume{spacers_mylar[]};}}; }
//Color AntiqueWhite { Volume{spacers_mylar[]}; }
//Color AntiqueWhite { Surface{Boundary{Volume{spacers_mylar_cut[]};}}; }
//Color AntiqueWhite { Volume{spacers_mylar_cut[]}; }

//Color SandyBrown { Surface{Boundary{Volume{spacers_tapet[]};}}; }
//Color SandyBrown { Volume{spacers_tapet[]}; }
//Color SandyBrown { Surface{Boundary{Volume{spacers_tapet_cut[]};}}; }
//Color SandyBrown { Volume{spacers_tapet_cut[]}; }

//Color SandyBrown { Surface{Boundary{Volume{spacers_tapeb[]};}}; }
//Color SandyBrown { Volume{spacers_tapeb[]}; }
//Color SandyBrown { Surface{Boundary{Volume{spacers_tapeb_cut[]};}}; }
//Color SandyBrown { Volume{spacers_tapeb_cut[]}; }

//Color SandyBrown { Surface{Boundary{Volume{frame};}}; }
//Color SandyBrown { Volume{frame}; }
//Color SandyBrown { Surface{Boundary{Volume{frame_cut};}}; }
//Color SandyBrown { Volume{frame_cut}; }

//Color Green { Surface{Boundary{Physical Volume{9};}}; }
//Color Green { Physical Volume{9}; }
//Color Green { Surface{Boundary{Physical Volume{19};}}; }
//Color Green { Physical Volume{19}; }

//Color Green { Surface{Boundary{Volume{gasv};}}; }
//Color Green { Volume{gasv}; }
//Color Green { Surface{Boundary{Volume{gasv_cut};}}; }
//Color Green { Volume{gasv_cut}; }

//Recursive Hide{ Volume{electrodet, electrodeb}; }
//Recursive Hide{ Volume{coatings[]}; }
Recursive Hide{ Volume{frame}; }
//Recursive Hide{ Volume{planes[]}; }
//Recursive Hide{ Volume{spacers[]}; }
//Recursive Hide{ Volume{spacers_tapet[]}; }
//Recursive Hide{ Volume{spacers_tapeb[]}; }
//Recursive Hide{ Volume{gasv}; }

//Recursive Hide{ Volume{electrodet_cut, electrodeb_cut}; }
//Recursive Hide{ Volume{coatings_cut[]}; }
//Recursive Hide{ Volume{frame_cut}; }
//Recursive Hide{ Volume{planes_cut[]}; }
//Recursive Hide{ Volume{spacers_mylar_cut[]}; }
//Recursive Hide{ Volume{spacers_tapet_cut[]}; }
//Recursive Hide{ Volume{spacers_tapeb_cut[]}; }
////Recursive Hide{ Volume{gasv_cut}; }

//Recursive Delete{ Volume{electrodes[]}; }
//Recursive Delete{ Volume{coatings[]}; }
//Recursive Delete{ Volume{frame}; }
//Recursive Delete{ Volume{planes[]}; }
//Recursive Delete{ Volume{spacers_mylar[]}; }
//Recursive Delete{ Volume{spacers_tapet[]}; }
//Recursive Delete{ Volume{spacers_tapeb[]}; }

Mesh.Format = 1;
Save "multigap.msh";

//r() = BooleanIntersection{Volume{electrodes[1]};}{Volume{coatings[1]};};
//Printf("'%g'",#r()); // 1-overlap, 2-do not overlap
