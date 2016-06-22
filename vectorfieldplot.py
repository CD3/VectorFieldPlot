#! /bin/env python
# -*- coding: utf8 -*-
 
'''
VectorFieldPlot - plots electric and magnetic fieldlines in svg
http://commons.wikimedia.org/wiki/User:Geek3/VectorFieldPlot
 
Copyright (C) 2010 Geek3
Copyright (C) 2015 CD3
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation;
either version 3 of the License, or (at your option) any later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this program; if not, see http://www.gnu.org/licenses/

Todo:

  Handle point source termination (more elegantly)
  Field line "cleaning". Need a way to make sure line spacing actually reflects field strength.

'''
 
import sys
from math import *
from cmath import phase

import svgwrite
from svg.path import Path, Line, parse_path
import scipy as sc
import scipy.integrate as ig
 
 
version = '2.0-alpha'
 
# some helper functions
def vabs(x):
    '''
    euclidian vector norm for any kind of vector
    '''
    return sqrt(sum([i**2 for i in x]))
 
def vnorm(x):
    '''
    vector normalisation
    '''
    d = vabs(x)
    if d != 0.: return sc.array(x) / d
    return sc.array(x)
 
def rot(xy, phi):
    '''
    2D vector rotation
    '''
    s = sin(phi); c = cos(phi)
    return sc.array([c * xy[0] - s * xy[1], c * xy[1] + s * xy[0]])
 

def get_elem_by_id(parent,key):
  keys = key.split('.')
  key = keys[0]

  for e in parent.elements:
    if 'id' in e.attribs and e['id'] == key:
      if len(keys) > 1:
        return get_elem_by_id(e,'.'.join(keys[1:]))
      else:
        return e

  return None

def in_bounds(bounds,p):
  '''Check if a point is within a given bounding box.'''
  if bounds is None:
    return True
  return p[0] >= bounds['x0'] and p[0] <= bounds['x1'] and p[1] >= bounds['y0'] and p[1] <= bounds['y1']
 
 







class FieldPlotDocument:
    '''
    Class representing an image.
    '''
    def __init__ (self, name, width=800, height=600, unit=100, center=None, license='CC-BY', author="UNKNOWN"):
        self.name = name
        self.width = float(width)
        self.height = float(height)
        self.unit = float(unit)
        self.license = license
        if center == None: self.center = [width / 2., height / 2.]
        else: self.center = [float(i) for i in center]




        # the SVG document
        self.dwg = svgwrite.Drawing(filename=name+'.svg', height=self.height, width=self.width)
        dwg = self.dwg
        # defs
        grad = dwg.radialGradient(id='glare',center=(0.65,0.7), r=0.75)
        grad.add_stop_color( color='#ffffff', offset='0',   opacity='0.7')
        grad.add_stop_color( color='#ffffff', offset='0.1', opacity='0.5')
        grad.add_stop_color( color='#ffffff', offset='0.6', opacity='0')
        grad.add_stop_color( color='#000000', offset='0.6', opacity='0')
        grad.add_stop_color( color='#000000', offset='0.75',opacity='0.05')
        grad.add_stop_color( color='#000000', offset='0.85',opacity='0.15')
        grad.add_stop_color( color='#000000', offset='1',   opacity='0.5')
        dwg.defs.add( grad )

        arrows = dwg.defs.add(dwg.g(id='arrows'))

        # background
        r = dwg.rect( id='background', insert=(0,0), size=(self.width,self.height), fill='white' )
        dwg.add(r)

        # image group
        self.img = dwg.g(id='image')
        dwg.add(self.img)
        self.img.translate(self.center[0], self.center[1])
        self.img.scale( self.unit, -self.unit )

        # add containers (groups for different elements)
        self.img.add(dwg.g(id='fieldlines'))
        self.img.add(dwg.g(id='equipotentiallines'))
        # make sure sources are on top
        self.img.add(dwg.g(id='sources'))

 
        # title and description
        desc = ''
        desc += 'created by '+author+' using VectorFieldPlot ' + version + ' (https://github.com/CD3/VectorFieldPlot)\n'
        desc += 'license: '+self.license+' (https://creativecommons.org/licenses/)\n'
        dwg.set_desc( title=name, desc=desc ) 



    def _add_arrow(self,color):
      '''Add an arrow for the given color to the document defs.'''
      id = color+'_arrow'
      if get_elem_by_id( self.dwg.defs, 'arrows.'+id ) is None:
        arrow_geo = {'x_nock':0.3,'x_head':3.8,'x_tail':-2.2,'width':4.5}
        path = parse_path( 'M 0.3,0 L -2.2,2.2 L 3.8,0 L -2.2,-2.2 Z' )
        arrow = self.dwg.path(path.d(), id=id, stroke='none', fill=color)
        arrow.scale(1./self.unit)

        get_elem_by_id( self.dwg.defs, 'arrows' ).add(arrow)

    def _get_bounds(self):
      '''Get the bounds for the image.'''
      bounds = {}
      bounds['x0'] = -self.center[0] / self.unit
      bounds['y0'] = -(self.height - self.center[1]) / self.unit
      bounds['x1'] = (self.width - self.center[0]) / self.unit
      bounds['y1'] =  self.center[1] / self.unit

      return bounds

    def _make_PointCharge_drawing(self, source, scale):
      '''Create an SVG element for a point charge.'''
      dwg = self.dwg

      g = dwg.g()
      c = dwg.circle(r=14, fill=('blue' if source.q < 0 else 'red') )
      g.add( c )
      c = dwg.circle(r=14,fill='url(#glare)',stroke='black',stroke_width=2)
      g.add( c )
      p = dwg.path()
      if source.q < 0:
        p.push('M 8,2 H -8 V -2 H 8 V 2 Z')
      else:
        p.push('M 2,2 V 8 H -2 V 2 H -8 V -2')
        p.push(' H -2 V -8 H 2 V -2 H 8 V 2 H 2 Z')
      g.add(p)
      g.translate(*source.pos())
      g.scale(float(scale)/self.unit)

      return g

    def draw_sources(self, sources, scale=1., stypes='All'):
      '''Draw the sources on the image.'''
      dwg = self.dwg
      img = self.img


      container = get_elem_by_id(self.img,'sources')
      for i,source in enumerate(sources.sources):
        if stypes == 'All' or isinstance( source, stypes ):
          g = dwg.g(id="None")
          try:
            g = getattr(self, '_make_'+source.__class__.__name__+'_drawing')( source, scale )
          except:
            pass
          container.add( g )



    def _make_line_drawing(self,line,linewidth,linecolor):
      bounds = self._get_bounds()
      line = line.get_line(bounds,)
      nodes = line['nodes']
      if len(nodes) < 2:
        nodes.append(nodes[0])
      path = Path()
      start = end = complex(*nodes[0]['p'])
      for i in range(1,len(nodes)):
        end = complex(*nodes[i]['p'])
        path.append( Line( start,end ) )
        start = end
      if line['closed']:
        path.append( Line( path[-1].end, path[0].start ) )
        path.closed = True

      line = self.dwg.path( path.d(), stroke=linecolor,stroke_width=linewidth/self.unit,fill='none' )

      return line

    def draw_fieldline(self, line, linewidth=2, linecolor='black', arrowstyle = {'num':1} ):
      '''Draw a field line on the drawing.'''
      container = get_elem_by_id(self.img,'fieldlines')
      line = self._make_line_drawing( line, linewidth, linecolor )
      group = self.dwg.g( id='fieldline{0}'.format( len(container.elements) ) )

      # now draw arrows (if needed)
      if arrowstyle['num'] > 0:
        self._add_arrow(linecolor)
        arrowscale = linewidth
        if 'scale' in arrowstyle:
          arrowscale *= arrowstyle['scale']
        path = parse_path( ' '.join(line.commands) )
        l = path.length()
        dl = l / (arrowstyle['num'] + 1)
        l = dl
        s = 0
        for segment in path:
          if s < l and s + segment.length() >= l:
            # put arrow at the center of the segment
            pos = (segment.start + segment.end)/2
            # get the displacement vector between the segments end points
            dir = segment.end - segment.start
            arrow = self.dwg.use( get_elem_by_id( self.dwg.defs, 'arrows.'+linecolor+'_arrow') )
            arrow.translate( pos.real, pos.imag )
            arrow.rotate(degrees(phase(dir)))
            arrow.scale( arrowscale )

            group.add(arrow)

          s += segment.length()


      group.add(line)
      container.add(group)

    def draw_fieldlines(self, lc, *args, **kwargs):
      for l in lc.get_lines():
        self.draw_fieldline(l,*args,**kwargs)

    def draw_equipotentialline(self, line, linewidth=2, linecolor='red'):
      '''Draw a equipotential line on the drawing.'''
      container = get_elem_by_id(self.img,'equipotentiallines')
      line = self._make_line_drawing( line, linewidth, linecolor )
      group = self.dwg.g( id='equipotentialline{0}'.format( len(container.elements) ) )
      group.add(line)
      container.add(group)


    def write(self, filename=None):
      '''Write the image.'''
      self.dwg.save()
      print 'image written to', self.dwg.filename






class FieldLine:
  '''Class that calculates field lines.'''

  def __init__(self, sources, p_start, backward=False):

    self.sources = sources
    self.p0 = sc.array(p_start)
    self.backward = backward

  def get_line(self, bounds=None, ds=1e-1, maxn=1000, ftype='E'):
    '''Get a set of points that lie on a field line passing through the point p.'''

    line = { 'closed' : False, 'nodes' : [] }
    nodes = line['nodes']

    # f is the function that returns the tangent vector for the line at all points in space
    try:
      field = getattr(self.sources,ftype)
    except:
      raise BaseException("ERROR: Unknown ftype "+ftype)

    dir = -1 if self.backward else 1
    f = lambda s, p: dir*vnorm(field(p))

    p = self.p0
    s = 0

    integrator = ig.ode( f )
    integrator.set_integrator('vode',method='bdf')
    integrator.set_initial_value(p,s)

    n = 0
    out = False
    while n <= maxn and in_bounds(bounds,p) and not line['closed'] and integrator.successful():
      s += ds
      plast = p
      p = integrator.integrate(s)
      nodes.append( {'p' : p } )
      n += 1
      if vabs(p-plast) < .9*ds:
        break
      if vabs(p - self.p0) > ds:
        out = True
      if out and vabs(p - self.p0) < ds:
        line['closed'] = True

    if not line['closed'] and n < maxn:
      # we probably ran out of bounds
      # try to go backwards
      p = self.p0
      integrator.set_initial_value(p,s)
      p_last = nodes[-1]['p']
      while n <= maxn and in_bounds(bounds,p) and not line['closed'] and integrator.successful():
        s -= ds
        plast = p
        p = integrator.integrate(s)
        nodes.insert(0, {'p' : p } )
        n += 1
        if vabs(p-plast) < .9*ds:
          break
        if out and vabs(p - p_last) < ds:
          line['closed'] = True

      if self.backward:
        line['nodes'].reverse()

    return line

class EquipotentialLine:
  '''Class that calculates equipotential lines.'''

  def __init__(self, sources, p_start):

    self.sources = sources
    self.p0 = sc.array(p_start)

  def get_line(self, bounds=None, ds=1e-2, maxn=1000):
    '''Get a set of points that lie on an equipotential line passing through the point p.'''

    line = { 'closed' : False, 'nodes' : [] }
    nodes = line['nodes']

    # equipotential lines are perpindicular to field lines (and therefore the field vectors)
    # so the "direction field" for equipotential lines is just the vector field rotated by 90 degrees.

    # f is the function that returns the tangent vector for the line at all points in space
    f = lambda s, p: rot( (vnorm(self.sources.E(p))), pi/2 )
    p = self.p0
    s = 0

    integrator = ig.ode( f )
    integrator.set_integrator('lsoda')
    integrator.set_initial_value(p,s)

    n = 0
    out = False
    while n <= maxn and in_bounds(bounds,p) and not line['closed']:
      s += ds
      p = integrator.integrate(s)
      nodes.append( {'p' : p } )
      n += 1
      if vabs(p - self.p0) > ds:
        out = True
      if out and vabs(p - self.p0) < ds:
        line['closed'] = True

    if not line['closed'] and n < maxn:
      # we probably ran out of bounds
      # try to go backwards
      p = self.p0
      integrator.set_initial_value(p,s)
      p_last = nodes[-1]['p']
      while n <= maxn and in_bounds(bounds,p) and not line['closed']:
        s -= ds
        p = integrator.integrate(s)
        nodes.insert(0, {'p' : p } )
        n += 1
        if vabs(p - p_last) < ds:
          line['closed'] = True


    return line



PoleSources = []

class Source(object):
  '''Class represnting a field source.'''
  def E(self,r):
    '''Returns the electric field due to the source at a given point.'''
    return sc.array([0,0])
  def B(self,r):
    '''Returns the magnetic field due to the source at a given point.'''
    return sc.array([0,0])
  def V(self,r):
    '''Returns the electric potential (with respect to ground at infinity) due to the source at a given point.'''
    return 0

class PointCharge(Source):
  '''A point charge source.'''
  def __init__(self, r, q=1 ):
    self.r = sc.array(r)
    self.q = q

  def E(self,r):
    Exy = sc.zeros(2)
    r = sc.array(r)
    dr = r - self.r
    d = vabs(dr)
    if d != 0.:
      Exy += self.q * dr / d**3

    return Exy

  def pos(self):
    return self.r

PoleSources.append(PointCharge)

class SourceCollection:
  '''A collection of field sources.'''
  def __init__(self):
    self.sources = []

  def add_source(self,s):
    self.sources.append(s)

  def get_nearest_pole(self,r):
    ns = None
    nd = None
    for s in self.sources:
      if isinstance( s, PoleSources ):
        d = vabs(s.pos() - r)
        if nd is None:
          nd = 2*d

        if d < nd:
          ns = s
          nd = d

    return ns

  def E(self,r):
    Exy = sc.zeros(2)
    for s in self.sources:
      Exy += s.E(r)
    return Exy

  def B(self,r):
    Bxy = sc.zeros(2)
    for s in self.sources:
      Bxy += s.B(r)
    return Bxy

  def V(self,r):
    Vxy = 0
    for s in self.sources:
      Vxy += s.V(r)
    return Vxy






class LineCollection(object):
  def __init__(self,*args,**kwargs):
    self.lines = []

  def get_lines(self):
    return self.lines



class FieldLineCollection(LineCollection):
  '''A collection of field lines with methods to help create them.'''
  def __init__(self,*args,**kwargs):
    super(FieldLineCollection,self).__init__(*args,**kwargs)


  def add_lines_to_point_charges( self, sources, numlines, halo=1e-1, theta0 = 0 ):
    '''Creates lines that will leave (or enter) all point sources at uniformly spaced angles.
       For example, to make sure that at least 10 lines leaving (or enter) each point charge in
       the drawing.
    '''

    lines = []
    for s in sources.sources:
      if isinstance(s, PointCharge):
        rs = s.pos()
        dtheta = 2*pi/numlines
        for i in range(numlines):
          theta = theta0 + i*dtheta
          dr = halo*sc.array( [cos(theta), sin(theta)] )
          lines.append( FieldLine( sources, rs + dr, s.q < 0 ) )

    self.lines += lines

