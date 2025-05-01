// Glass RPC prototype

cut = 0;
break = 0;
span = 0.2;

SetFactory("OpenCASCADE");

Mesh.ElementOrder = 2;
Mesh.SurfaceFaces = 1;
Mesh.SurfaceEdges = 0;
Mesh.VolumeFaces = 1;
Mesh.VolumeEdges = 0;
Geometry.Points = 0;

Geometry.OCCTargetUnit = "M";

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

gwidth = ngap*glgapw+(ngap+1)*gap;
Printf("Gap width is '%g'", gwidth);

frame() = ShapeFromFile("frame.step");

bbox() = BoundingBox Volume{frame()};
xmin = bbox(0);
ymin = bbox(1);
zmin = bbox(2);
xmax = bbox(3);
ymax = bbox(4);
zmax = bbox(5);

Printf("Import frame with dimensions (meter) '%g','%g' '%g','%g' '%g','%g'", xmin, xmax, ymin, ymax, zmin, zmax);

Rotate {{1,0,0}, {0,0,0}, -Pi/2} { Volume{frame}; }
Translate {0.5*glelx, 0.5*glely, coatw+glelw+gwidth} { Volume{frame}; }

det = newv;
Box(det) = {0,0,glelw+coatw, glelx,glely,gwidth};
gasv = newv;
BooleanDifference(gasv) = {Volume{det}; Delete;}{Volume{frame};};

sph = gwidth;

coatings[0] = newv;
Box(coatings[0]) = {(glelx-coatx)*0.5,(glely-coaty)*0.5,0,coatx,coaty,coatw};
coatings += Translate {0, 0, glelw+coatw+ngap*(glgapw+gap)+gap+glelw} { Duplicata{ Volume{coatings[0]}; }};

electrodet = newv;
Box(electrodet) = {0,0,coatw, glelx,glely,glelw};
electrodeb = Translate {0, 0, glelw+ngap*(glgapw+gap)+gap} { Duplicata{ Volume{electrodet}; }};

plane = newv;
Box(plane) = {0,0,coatw, glgapx,glgapy,glgapw};

planes = {};
For i In {0: ngap-1}
  planes += Translate {(glelx-glgapx)/2., (glely-glgapy)/2., i*(glgapw+gap)+glelw+gap} { Duplicata{ Volume{plane}; }};  
EndFor

spacer = newv;
Cylinder(spacer) = {0, 0, 0, 0, 0, sph, spr};

spacers = {};

lspx = coatx/nsp;
lspy = coaty/nsp;

For i In {1: nsp}
  For j In {1: nsp}
    spacers += Translate {(glelx-coatx)/2+lspx*i-lspx/2, (glely-coaty)/2+lspy*j-lspy/2, glelw+coatw} { Duplicata{ Volume{spacer}; }};
  EndFor
EndFor

Recursive Delete{ Volume{spacer}; }
Recursive Delete{ Volume{plane}; }

sp = {};
For s In {0: #spacers[]-1}
  For i In {0: #planes[]-1}
    j = #planes[]-1-i;
    If (j == #planes[]-1)
      v() = BooleanDifference{Volume{spacers[s]};}{Volume{planes[j]};};
      sp += v[0];
    Else
      v() = BooleanDifference{Volume{v[1]};Delete;}{Volume{planes[j]};};
      sp += v;
    EndIf
  EndFor
  sp += v[1];
EndFor

Recursive Delete{ Volume{spacers[]}; }

For i In {0: #planes[]-1}
  gasv0 = gasv; gasv = newv;
  BooleanDifference(gasv) = {Volume{gasv0}; Delete;}{Volume{planes[i]};};
EndFor

For i In {0: #sp[]-1}
  gasv0 = gasv; gasv = newv;
  BooleanDifference(gasv) = {Volume{gasv0}; Delete;}{Volume{sp[i]};};
EndFor

If (cut)
  bb = newv;
  Box(bb) = {0,0,0, (glelx-coatx)/2+lspx/2,glely,gwidth+2*glelw+2*coatw};

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

//  sp_cut = {};
//  For i In {0: #sp[]-1}
//    cut() = BooleanDifference{Volume{sp[i]};}{Volume{bb};};
//    sp_cut += cut();
//  EndFor

  Recursive Delete{ Volume{bb}; }
  Recursive Delete{ Volume{electrodet}; }
  Recursive Delete{ Volume{electrodeb}; }
  Recursive Delete{ Volume{planes[]}; }
  Recursive Delete{ Volume{frame}; }
  Recursive Delete{ Volume{gasv}; }
//  Recursive Delete{ Volume{sp[]}; }
  Recursive Delete{ Volume{coatings[]}; }
  
  electrodet = electrodet_cut;
  electrodeb = electrodeb_cut;
  planes[] = planes_cut[];
  frame = frame_cut;
  gasv = gasv_cut;
  coatings[] = coatings_cut[];
EndIf

If (break)
  elem = {};
  psg = newv;
  Physical Volume(psg) = {coatings[0]};
  elem += psg;
  psg = newv;
  Physical Volume(psg) = {electrodet};
  elem += psg;
  For i In {0: #planes[]-1}
    psg = newv;
    Physical Volume(psg) = {};
    For j In {0: #planes[]-1}
      Physical Volume(psg) += sp[i+(#planes[]+1)*j];
    EndFor
    elem += psg;
    psg = newv;
    Physical Volume(psg) = {planes[i]};
    elem += psg;
    If (i == #planes[]/2-1)
      psg = newv;
      Physical Volume(psg) = {gasv};
      elem += psg;
      psg = newv;
      Physical Volume(psg) = {frame};
      elem += psg;
    EndIf
    If (i == #planes[]-1)
      psg = newv;
      Physical Volume(psg) = {};
      For j In {0: #planes[]-1}
        Physical Volume(psg) += sp[i+1+(#planes[]+1)*j];
      EndFor
      elem += psg;    
    EndIf
  EndFor
  psg = newv;
  Physical Volume(psg) = {electrodeb};
  elem += psg;
  psg = newv;
  Physical Volume(psg) = {coatings[1]};
  elem += psg;
  sep = span/#elem[];
  For i In {0: #elem[]-1}
    Translate {0, 0, -sep*i}{Physical Volume{elem[i]};}
  EndFor  
EndIf

Coherence;

Mesh 3;

If (break)
  Recursive Hide{ Volume{gasv}; }
EndIf
Recursive Delete{ Surface{:}; }

//Recursive Hide{ Volume{planes[]}; }
//Recursive Hide{ Volume{electrodet, electrodeb}; }
//Recursive Hide{ Volume{coatings[]}; }
//Recursive Hide{ Volume{frame}; }
//Recursive Hide{ Volume{sp[]}; }
//Recursive Hide{ Volume{gasv}; }

If (cut)
  electrodet = electrodet_cut;
  electrodeb = electrodeb_cut;
  coatings[] = coatings_cut[];
  planes[] = planes_cut[];
//  sp[] = sp_cut[];
  frame = frame_cut;
////  gasv = gasv_cut;
EndIf

Physical Volume("electrodet", 1) = {electrodet};
Physical Volume("electrodeb", 2) = {electrodeb};
Physical Volume("coatings", 3) = {coatings[]};
Physical Volume("planes", 4) = {planes[]};
Physical Volume("spacers", 5) = {sp[]};
Physical Volume("frame", 6) = {frame};
Physical Volume("gasv", 7) = {gasv};

Color Snow { Surface{Boundary{Volume{electrodet, electrodeb};}}; }
Color Snow { Volume{electrodet, electrodeb}; }
Color DarkGrey { Surface{Boundary{Volume{coatings[]};}}; }
Color DarkGrey { Volume{coatings[]}; }
Color WhiteSmoke { Surface{Boundary{Volume{planes[]};}}; }
Color WhiteSmoke { Volume{planes[]}; }
Color AntiqueWhite { Surface{Boundary{Volume{sp[]};}}; }
Color AntiqueWhite { Volume{sp[]}; }
Color SandyBrown { Surface{Boundary{Volume{frame};}}; }
Color SandyBrown { Volume{frame}; }
Color Green { Surface{Boundary{Volume{gasv};}}; }
Color Green { Volume{gasv}; }

Mesh.Format = 1;
Save "multigap.msh";
