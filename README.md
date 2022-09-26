VectorFieldPlot
===============

A python module for creating svg images of electric and magnetic field lines for user defined charge and current configurations.

This module is an adaptation of the code written by Geek3 to create high quality, physically correct images of electric and magnetic field lines.

[http://commons.wikimedia.org/wiki/User:Geek3/VectorFieldPlot](http://commons.wikimedia.org/wiki/User:Geek3/VectorFieldPlot)

The version posted by Geek3 was a single python script that generated the field lines. It contained a set of classes and utility functions for creating the lines,
and then a short program at the end that used these classes to create an image for a specific charge configuration. The images created by the script are fantastic, it does a great job.

This project just organized the classes that do the work into their own module so that is can be used by multiple programs at once.
This initial commit is a straight copy of version 1.3 posted by Geek3.

# Installing

```bash
$ pip install vectorfieldplot
```

# Using

Here is an example program that uses the module to generate a picture of the fields lines for an electric dipole
consisting of a positive and negative charge on the x-axis.

```python
#!  /usr/bin/env python

'''An example script that generates a picture of the field lines for an electric dipole.'''

import math
import vectorfieldplot as vfp


# create a document. we specify the file name and image size here
doc = vfp.FieldplotDocument( 'ElectricDipole', width=800,height=600,unit=100)

# create a field opbject
field = vfp.Field()
# add the point charges
field.add_element('monopoles' , [ [ 1,0, 1]  # the positive charge
                                , [-1,0,-1]  # the negative charge
                                ] )

# draw the charges for the field on the document
doc.draw_charges(field)

# start drawing the field lines
# we are going to draw 20 field lines comming off of the positive charge at uniformly spaced angles.
N = 20
for i in range(N):
    # compute the angle that the line will start off at
    angle = i * 2.*math.pi / (N-1)
    # generate the line
    # this takes the initial position, the x and y components of the initial direction, and whether or not should go forward
    # or backward.
    line = vfp.FieldLine( field, [1,0], start_v=[math.cos( angle ) , math.sin( angle )],directions='forward' )
    # draw the lin on the document
    doc.draw_line(line,arrows_style={'min_arrows':1,'max_arrows':1})

# write the document
doc.write()
```


This example will draw the field for a dipole directly using the dipole element
```python
#!  /usr/bin/env python

'''An example script that generates a picture of the field lines for an electric dipole.'''

import math
import vectorfieldplot as vfp


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
```

## Messages and logging

By default, the software only prints out messages that are significant warnings or errors.
To change this default, you can do the following:

```python
import logging
logging.basicConfig(level=logging.DEBUG) # can be DEBUG, INFO, WARNING, or ERROR
```
