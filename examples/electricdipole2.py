#!  /usr/bin/env python

'''An example script that generates a picture of the field lines for an electric dipole.'''

import math
import logging
import vectorfieldplot as vfp

logging.basicConfig(level=logging.DEBUG)


# create a document. we specify the file name and image size here
doc = vfp.FieldplotDocument( 'ElectricDipole2', width=800,height=600,unit=100)

# create a field opbject
field = vfp.Field()
# add the dipole
# note, parameters are should be [ r_x, r_y, p_x, p_y ], where r is the position vector and p is the dipole moment.
field.add_element('dipoles' , [ [ 0,0,1,0] ] )

# draw the charges for the field on the document
doc.draw_charges(field)

# start drawing the field lines
# we are going to draw 20 field lines comming off of the positive charge at uniformly spaced angles.
N = 100
for i in range(N):
    # compute the angle that the line will start off at
    angle =  -math.pi/2 + i * math.pi / (N-1)
    # generate the line
    # this takes the initial position, the x and y components of the initial direction, and whether or not should go forward
    # or backward.
    line = vfp.FieldLine( field, [0.1*math.cos(angle),0.1*math.sin(angle)], start_v=[math.cos( angle ) , math.sin( angle )],directions='forward' )
    # draw the lin on the document
    doc.draw_line(line,arrows_style={'min_arrows':1,'max_arrows':1})

# write the document
doc.write()
