#!  /usr/bin/env python

'''An example script that generates a picture of the field lines for an electric dipole.'''

import math
import vectorfieldplot as vfp


doc = vfp.FieldplotDocument( 'ElectricFieldVectors', width=800,height=600,unit=100)
# create a field opbject
field = vfp.Field()
# add the point charges
field.add_element('monopoles' , [ [ 1,0, 1]  # the positive charge
                                , [-1,0,-1]  # the negative charge
                                ] )
doc.draw_charges(field)

line = vfp.FieldLine( field, [1,0], start_v=[0,1],directions='forward' )
doc.draw_line(line,arrows_style={'min_arrows':1,'max_arrows':1})
line = vfp.FieldLine( field, [1,0], start_v=[1,.1],directions='forward' )
doc.draw_line(line,arrows_style={'min_arrows':1,'max_arrows':1})

vectors = vfp.FieldVectors( field, maxscale=0.7, halo='auto' )
vectors.add_vectors_in_grid( -3.5, -2.5, 3.5, 3.5, 20, 10 )
doc.draw_vectors(vectors)

doc.write()
