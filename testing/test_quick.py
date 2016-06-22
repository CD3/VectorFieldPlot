import math
import numpy
import vectorfieldplot as vfp

def test_dipole():

  sources = vfp.SourceCollection()
  sources.add_source( vfp.PointCharge( [0, 1], 1 ) )
  sources.add_source( vfp.PointCharge( [0,-1],-1 ) )

  doc = vfp.FieldplotDocument( 'ElectricDipole', width=800,height=600,unit=100)
  doc.draw_sources(sources)

  N = 10
  for i in range(N):
      # compute the angle that the line will start off at
      angle = i * 2.*math.pi / (N-1)
      line = vfp.FieldLine( sources, [0+1e-1*math.cos(angle),1+1e-1*math.sin(angle) ] )
      doc.draw_fieldline(line)

  line = vfp.FieldLine( sources, [2,0] )
  # doc.draw_fieldline(line)
  # line = vfp.FieldLine( sources, [0.1,-2] )
  # doc.draw_fieldline(line)

  # line = vfp.EquipotentialLine( sources, [2,0] )
  # doc.draw_equipotentialline(line)
  line = vfp.EquipotentialLine( sources, [0.1,-2] )
  doc.draw_equipotentialline(line)

  doc.write()



