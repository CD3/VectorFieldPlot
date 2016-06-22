# VectorFieldPlot

A python module for creating svg images of electric and magnetic field lines for user defined charge and current configurations.

This module started out as an adaptation of the
[code written by Geek3](http://commons.wikimedia.org/wiki/User:Geek3/VectorFieldPlot),
to create high quality, physically correct images of electric and magnetic field lines.
It has since been significantly rewritten, and is not as complete as the original version yet.

**Rewrite**
I started the rewrite when I wanted to add support for drawing equipotential lines on the images. I spent several days trying to figure
out how the code worked and started modifying small pieces until eventually I decided to do a rewrite so that I could learn about svg images
and how to calculate field lines.

## Examples

An image of a field is constructed in parts using several classes:

  - The `FieldPlotDocument` class represents the svg file. To build an image you draw items on the document,
    and then write it.

  - The `Source` class is used to represent sources of the field. This is a baseclass that will not be used directly,
    but all field source classes inherit from it. The following sources are currently implemented:

      - `PointCharge` represents a single electric point charge.

  - The 'SourceCollection` class is used to build a collection of sources into a system.

  - The `FieldLine` class is used to compute a field line for a collection of point sources.

  - The `EquipotentialLine` class is used to compute an equipotential line for a collection of point sources.


This example draws the electric field lines for an electric dipole.

```
import vectorfieldplot as vfp

sources = vfp.SourceCollection()
sources.add_source( vfp.PointCharge( [0, 1], 1 ) )
sources.add_source( vfp.PointCharge( [0,-1],-1 ) )

doc = vfp.FieldPlotDocument( 'ElectricDipole', width=800,height=600,unit=100)
doc.draw_sources(sources)

N = 10
for i in range(N):
  # compute the angle that the line will start off at
  angle = i * 2.*math.pi / (N-1)
  p = [1,0]
  p[0] = 1e-1*math.cos(angle)
  p[1] = 1e-1*math.sin(angle)
  line = vfp.FieldLine( field, p)
  doc.draw_fieldline(line,arrows_style={'min_arrows':1,'max_arrows':1})

doc.write()
```

This example does the same thing, but uses the FieldLineCollection class to automatically genertae the lines
that leave (or enter) each point charge.

```
import vectorfieldplot as vfp

sources = vfp.SourceCollection()
sources.add_source( vfp.PointCharge( [0, 1], q= 1 ) )
sources.add_source( vfp.PointCharge( [0,-1], q=-1 ) )

doc = vfp.FieldPlotDocument( 'ElectricDipole', width=800,height=600,unit=100)
doc.draw_sources(sources)

fieldlines = vfp.FieldLineCollection()
fieldlines.add_lines_to_point_charges( sources, 10 )

doc.draw_fieldlines( fieldlines )

doc.write()
```

More examples can be found in the [testing](testing) directory.

## Sources

To add support for new field sources (such as an infinite line of charge) you just need to create a new class that derives the Source class and implements
the `E` or `B` (or both) methods. The FieldLine class will use these methods with a single argument that is the coordinate to evaluate the field at,
and that is all that is required for the library to draw the field lines. You can also
add support for drawing the source itself by implementing a method named `_make_<ClassName>_drawing` and adding it to the FieldPlotDocument class. This method will be passed
the source instance and a scale number. It should return an SVG element that will draw an image of the source in the correct spot.
