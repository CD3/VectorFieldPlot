VectorFieldPlot
===============

A python module for creating svg images of electric and magnetic field lines for user defined charge and current configurations.

This module is an adaptation of the code written by Geek3 to create high quality, physically correct images of electric and magnetic field lines.

[http://commons.wikimedia.org/wiki/User:Geek3/VectorFieldPlot](http://commons.wikimedia.org/wiki/User:Geek3/VectorFieldPlot)

The version posted by Geek3 was a single python script that generated the field lines. It contained a set of classes and utility functions for creating the lines,
and then a short program at the end that used these classes to create an image for a specific charge configuration. The images created by the script are fantastic, it does a great job.

This project just organized the classes that do the work into their own module so that is can be used by multiple programs at once.
This initial commit is a straight copy of version 1.3 posted by Geek3.


