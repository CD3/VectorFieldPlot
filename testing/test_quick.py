import math, cmath
import numpy
import vectorfieldplot as vfp

class MySource(vfp.PointCharge):
  pass

def test_dipole():

  sources = vfp.SourceCollection()
  sources.add_source( vfp.Monopole( [ 1, 0], q = -1 ) )
  sources.add_source( vfp.Monopole( [-1, 0], q =  1 ) )
  # sources.add_source( vfp.Dipole( [-0.5, 0.5], [0,1] ) )
  # sources.add_source( vfp.Dipole( [0.5, 0], [0,1] ) )


  doc = vfp.FieldPlotDocument( 'ElectricDipole', width=800,height=600,unit=100)
  doc.draw_sources(sources)

  fieldlines = vfp.FieldLineCollection()
  # fieldlines.add_lines_around_point( sources, [1,0], numlines=10, halo=1e-1, delta_theta=math.pi/4, theta0=math.radians(45) )
  # fieldlines.add_lines_around_sources( sources, numlines=10, delta_theta=2*math.pi, filter = lambda s : s.q < 0 )
  # fieldlines.add_lines_around_sources( sources, numlines=4 , delta_theta=4*(2*math.pi/10), theta0 = math.pi, filter = lambda s : s.q > 0 )
  # fieldlines.add_lines_to_point_charges( sources, numlines=10 )
  # fieldlines.add_lines_to_dipole_charges( sources, 5, theta = math.pi/4, halo = 0.1 )

  fieldlines.add_lines_around_sources( sources.filtered(lambda s : s.q > 0), numlines=3 )
  fieldlines.add_lines_around_sources( sources.filtered(lambda s : s.q < 0), numlines=5 )
  fieldlines.add_lines_around_sources( sources, numlines=5, filterfunc = lambda s : s.q > 0 )

  doc.draw_fieldlines( fieldlines, config = {'halo' : 0.1} )
  # doc.draw_fieldlines( fieldlines, config = {'halo' : 0.1} )

  # equipotentiallines = vfp.EquipotentialLineCollection()
  # equipotentiallines.add_lines_along_line( sources, [ [1,0], [-1,0] ], halo=2e-1, numlines=10 )

  # doc.draw_equipotentiallines( equipotentiallines )



  doc.write()




