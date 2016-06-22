import math, cmath
import numpy
import vectorfieldplot as vfp

def test_utils():

  doc = vfp.svgwrite.Drawing()
  g = doc.g(id='group1')
  gg = doc.g(id='group2')
  g.add(gg)
  doc.add( g )
  print doc.elements

  assert id(vfp.get_elem_by_id(doc,'group1') ) == id(g)
  assert id(vfp.get_elem_by_id(doc,'group1.group2') ) == id(gg)

def test_dipole():
  sources = vfp.SourceCollection()
  sources.add_source( vfp.PointCharge( [0, 1], 1 ) )
  sources.add_source( vfp.PointCharge( [0,-1],-1 ) )

  doc = vfp.FieldplotDocument( 'ElectricDipole', width=800,height=600,unit=100)
  doc.draw_sources(sources)

  N = 20
  for i in range(N):
    # compute the angle that the line will start off at
    angle = i * 2.*math.pi / (N-1)
    p = [0,1]
    p[0] += 1e-1*math.cos(angle)
    p[1] += 1e-1*math.sin(angle)
    line = vfp.FieldLine( sources, p)
    doc.draw_fieldline(line)

  doc.write()



def test_linecharge():
  sources = vfp.SourceCollection()
  X = numpy.arange(-2,2,0.1)
  for x in X:
    sources.add_source( vfp.PointCharge( [x, 0], 1 ) )

  doc = vfp.FieldplotDocument( 'LineCharge', width=800,height=600,unit=100)
  doc.draw_sources(sources)

  for x in X:
    line = vfp.FieldLine( sources, [x,0.1] )
    doc.draw_fieldline(line)
    line = vfp.FieldLine( sources, [x,-0.1] )
    doc.draw_fieldline(line)

  doc.write()
