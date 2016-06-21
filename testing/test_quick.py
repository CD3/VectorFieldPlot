import math
import numpy
import vectorfieldplot as vfp

def test_dipole():
  # field = vfp.Field()
  # field.add_element('monopoles' , [ [ 1,0, 1]  # the positive charge
                                  # , [-1,0,-1]  # the negative charge
                                  # ] )
  # doc = vfp.FieldplotDocument( 'ElectricDipole', width=800,height=600,unit=100)
  # doc.draw_charges(field)

  # N = 10
  # for i in range(N):
      # # compute the angle that the line will start off at
      # angle = i * 2.*math.pi / (N-1)
      # line = vfp.FieldLine( field, [1,0], start_v=[math.cos( angle ) , math.sin( angle )],directions='forward' )
      # doc.draw_fieldline(line,arrows_style={'min_arrows':1,'max_arrows':1})

  # line = vfp.EquipotentialLine( field, [2,0] )
  # doc.draw_equipotentialline(line)
  # line = vfp.EquipotentialLine( field, [0.1,-2] )
  # doc.draw_equipotentialline(line)


  # doc.write()


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
      # doc.draw_fieldline(line)

  line = vfp.FieldLine( sources, [2,0] )
  # doc.draw_fieldline(line)
  # line = vfp.FieldLine( sources, [0.1,-2] )
  # doc.draw_fieldline(line)

  # line = vfp.EquipotentialLine( sources, [2,0] )
  # doc.draw_equipotentialline(line)
  line = vfp.EquipotentialLine( sources, [0.1,-2] )
  doc.draw_equipotentialline(line)

  doc.write()



