#!/Users/kskovpen/.juliaup/bin/julia

import Gmsh: gmsh

gmsh.initialize()

#gmsh.option.setNumber("General.Verbosity", 2)

gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
gmsh.option.setNumber("Mesh.SurfaceEdges", 0)
gmsh.option.setNumber("Mesh.VolumeFaces", 1)
gmsh.option.setNumber("Mesh.VolumeEdges", 0)
gmsh.option.setNumber("Geometry.Points", 0)

gmsh.merge("multigap.msh")

phys = gmsh.model.getPhysicalGroups()
for gr in phys
  name = gmsh.model.getPhysicalName(gr[1], gr[2])
  if name == ""
    continue
  end
  ents = gmsh.model.getEntitiesForPhysicalName(name)
  if occursin("electrode", name)
    for ent in ents
      gmsh.model.setColor(ent, 163, 204, 231)
    end
  end
  if occursin("coating", name)
    for ent in ents
      gmsh.model.setColor(ent, 102, 111, 117)
    end
  end
  if occursin("frame", name)
    for ent in ents
      gmsh.model.setColor(ent, 244, 164,  96, 255)
    end
  end
  if occursin("gasv", name)
    for ent in ents
      gmsh.model.setColor(ent, 0, 255,   0, 255)
    end
  end
  if occursin("spacer", name)
    for ent in ents
      gmsh.model.setColor(ent, 250, 235, 215, 255)
    end
  end
  if occursin("plane", name)
    for ent in ents
      gmsh.model.setColor(ent, 245, 245, 245, 255)
    end
  end
end

#gmsh.model.removeEntities(ents)

gmsh.option.setNumber("General.Trackball", 0)
gmsh.option.setNumber("General.RotationX", -65)
gmsh.option.setNumber("General.RotationY", 0)
gmsh.option.setNumber("General.RotationZ", 125)
gmsh.option.setNumber("General.Light0", 1)
gmsh.option.setNumber("General.Light0X", 0)
gmsh.option.setNumber("General.Light0Y", 0)
gmsh.option.setNumber("General.Light0Z", 1)

gmsh.fltk.initialize()

gmsh.option.setNumber("Print.Width", 2000)
gmsh.write("img/cutview.png")
gmsh.finalize()
