import math, cmath
import numpy
import vectorfieldplot as vfp

class MySource(vfp.PointCharge):
  pass

def test_dipole():

  sources = vfp.SourceCollection()
  # sources.add_source( vfp.Monopole( [ 1, 0], q = -1 ) )
  # sources.add_source( vfp.Monopole( [-1, 0], q =  1 ) )
  sources.add_source( vfp.Dipole( [0, 0], [0,1] ) )

  doc = vfp.FieldPlotDocument( 'ElectricDipole', width=800,height=600,unit=100)
  doc.draw_sources(sources)

  fieldlines = vfp.FieldLineCollection()
  fieldlines.add_lines_to_point_charges( sources, 10 )
  fieldlines.add_lines_to_dipole_charges( sources, 10 )

  doc.draw_fieldlines( fieldlines )



  doc.write()




