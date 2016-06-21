import math
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
  field = vfp.Field()
  field.add_element('monopoles' , [ [ 1,0, 1]  # the positive charge
                                  , [-1,0,-1]  # the negative charge
                                  ] )
  doc = vfp.FieldplotDocument( 'ElectricDipole', width=800,height=600,unit=100)
  doc.draw_charges(field)

  N = 20
  for i in range(N):
      # compute the angle that the line will start off at
      angle = i * 2.*math.pi / (N-1)
      line = vfp.FieldLine( field, [1,0], start_v=[math.cos( angle ) , math.sin( angle )],directions='forward' )
      doc.draw_fieldline(line,arrows_style={'min_arrows':1,'max_arrows':1})

  doc.write()


  field = vfp.Field()
  field.add_element('dipoles' , [ [ 1,0, 1, 0] ] )
  doc = vfp.FieldplotDocument( 'ElectricDipoleApprox', width=800,height=600,unit=100)

  N = 20
  for i in range(N):
      # compute the angle that the line will start off at
      angle = i * 2.*math.pi / (N-1)
      line = vfp.FieldLine( field, start_p=[1,0], start_d=[0.3*math.cos( angle ) , 0.3*math.sin( angle )], directions='both' )
      doc.draw_line(line,arrows_style={'min_arrows':1,'max_arrows':1})

  doc.write()



def test_linecharge():
  field = vfp.Field()
  X = numpy.arange(-2,2,1)
  field.add_element('monopoles' , [ [x, 0, 1] for x in X ] )

  doc = vfp.FieldplotDocument( 'LineCharge', width=800,height=600,unit=100)
  doc.draw_charges(field)

  for x in X:
    line = vfp.FieldLine( field, [x,0.1], directions='both' )
    doc.draw_line(line,arrows_style={'min_arrows':1,'max_arrows':1})
    line = vfp.FieldLine( field, [x,-0.1], directions='both' )
    doc.draw_line(line,arrows_style={'min_arrows':1,'max_arrows':1})

  doc.write()
