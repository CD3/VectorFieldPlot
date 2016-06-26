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

  Handle source termination (more elegantly)
  Field line "cleaning". Need a way to make sure line spacing actually reflects field strength.
  add minn to line classes

'''
 
import sys, inspect
from math import *
from cmath import phase

import dpath.util
import svgwrite
import svg.path
import scipy as sc
import scipy.integrate as ig
import scipy.optimize as op
 
 
version = '2.0-alpha'
 
# some helper functions
def vmag(x):
    '''
    euclidian vector norm for any kind of vector
    '''
    return sqrt(sum([i**2 for i in x]))

vabs = vmag

def vdir(x):
    '''
    direction of 2D vector
    '''
    return atan2( x[1], x[0] )
 
def vnorm(x):
  '''
  vector normalization
  '''
  d = vabs(x)
  if d != 0.: return sc.array(x) / d
  return sc.array(x)
 
def vrot(xy, phi):
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
 
 

arrowhead = svg.path.parse_path( 'M 0.3,0 L -2.2,2.2 L 3.8,0 L -2.2,-2.2 Z' )






class FieldPlotDocument(object):
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

        arrowhead_marker = dwg.marker( id='black_arrow_marker', orient='auto', insert=(1,5), markerWidth=3, markerHeight=3, viewBox=(0,0,10,10) )
        arrowhead_marker.add( dwg.path('M 0 0 L 10 5 L 0 10 z') )
        dwg.defs.add( arrowhead_marker )


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
        path = arrowhead
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

    def draw_sources(self, sources, config={}):
      '''Draw the sources on the image.'''
      scale = config.get('scale',1.)
      stypes = config.get('stypes','All')

      dwg = self.dwg
      img = self.img


      container = get_elem_by_id(self.img,'sources')
      for i,source in enumerate(sources.sources):
        if stypes == 'All' or isinstance( source, stypes ):
          g = dwg.g(id="None")
          # look for a drawing method for the source or any of its base classes
          for cls in inspect.getmro(source.__class__):
            try:
              mn = '_make_'+cls.__name__+'_drawing'
              g = getattr(self, mn)( source, scale )
              break
            except:
              pass
          if not g is None:
            container.add( g )



    def _make_line_drawing(self, line, config={}):
      linewidth = config.get('linewidth', 2)
      linecolor = config.get('linecolor','black')

      nodes = line['nodes']
      if len(nodes) < 2:
        nodes.append(nodes[0])
      path = svg.path.Path()
      start = end = complex(*nodes[0]['p'])
      for i in range(1,len(nodes)):
        end = complex(*nodes[i]['p'])
        path.append( svg.path.Line( start,end ) )
        start = end
      if line['closed']:
        path.append( svg.path.Line( path[-1].end, path[0].start ) )
        path.closed = True

      line = self.dwg.path( path.d(), stroke=linecolor,stroke_width=linewidth/self.unit,fill='none' )

      return line

    def draw_fieldline(self, line, config={}):
      '''Draw a field line on the drawing.'''
      container = get_elem_by_id(self.img,'fieldlines')

      linewidth = config.get('linewidth', 2)
      linecolor = config.get('linecolor', 'black')

      line = line.get_line(self._get_bounds(), config)
      line = self._make_line_drawing( line, config )
      group = self.dwg.g( id='fieldline{0}'.format( len(container.elements) ) )

      arrowstyle = config.get('arrowstyle', {'spacing':1})
      # now draw arrows (if needed)
      if not arrowstyle is None:
        self._add_arrow(linecolor)
        arrowscale = linewidth
        if 'scale' in arrowstyle:
          arrowscale *= arrowstyle['scale']

        if not 'num' in arrowstyle:
          arrowstyle['num'] = 1


        path       = svg.path.parse_path( ' '.join(line.commands) )
        pathlength = path.length()
        if not 'spacing' in arrowstyle:
          arrowstyle['spacing'] = pathlength / (arrowstyle['num'] + 1)

        spacing = arrowstyle['spacing']
        if spacing < pathlength/2:
          spacing = pathlength / 2.

        arrowpos = spacing
        s = 0
        for segment in path:
          if s < arrowpos and s + segment.length() >= arrowpos:
            # put arrow at the center of the segment
            pos = (segment.start + segment.end)/2
            # get the displacement vector between the segments end points
            dir = segment.end - segment.start
            arrow = self.dwg.use( get_elem_by_id( self.dwg.defs, 'arrows.'+linecolor+'_arrow') )
            arrow.translate( pos.real, pos.imag )
            arrow.rotate(degrees(phase(dir)))
            arrow.scale( arrowscale )
            group.add(arrow)

            arrowpos += spacing

          s += segment.length()


      group.add(line)
      container.add(group)

    def draw_fieldlines(self, lc, config={}):
      for l in lc.get_lines():
        self.draw_fieldline(l, config)

    def draw_equipotentialline(self, line, config={}):
      '''Draw a equipotential line on the drawing.'''
      container = get_elem_by_id(self.img,'equipotentiallines')

      config['linewidth'] = config.get('linewidth', 2)
      config['linecolor'] = config.get('linecolor', 'red')
      line = line.get_line(self._get_bounds(),config)
      line = self._make_line_drawing( line, config )
      group = self.dwg.g( id='equipotentialline{0}'.format( len(container.elements) ) )
      group.add(line)
      container.add(group)

    def draw_equipotentiallines(self, lc, config={}):
      for l in lc.get_lines():
        self.draw_equipotentialline(l, config)



    def write(self, filename=None):
      '''Write the image.'''
      self.dwg.save()
      print 'image written to', self.dwg.filename



class Line(object):
  def __init__(self,*args,**kwargs):
    pass

  def get_line(self, bounds, config={}):
    return None

class FieldLine(Line):
  '''Class that calculates field lines.'''

  def __init__(self, sources, p_start, backward=False):
    super(FieldLine,self).__init__()
    self.sources = sources
    self.p0 = sc.array(p_start)
    self.backward = backward

  def get_line(self, bounds=None, config={}):
    '''Get a set of points that lie on a field line passing through the point p.'''

    ftype = config.get('ftype', 'E')
    ds    = config.get('ds',1e-1)
    maxn  = config.get('maxn',1000)
    minn  = config.get('minn',10)
    halo  = config.get('halo',ds)

    line = { 'closed' : False, 'nodes' : [] }
    nodes = line['nodes']

    # f is the function that returns the tangent vector for the line at all points in space
    try:
      field = getattr(self.sources,ftype)
    except:
      raise BaseException("ERROR: Unknown ftype "+ftype)

    dir = -1 if self.backward else 1
    f = lambda s, p: dir*vnorm(field(p))

    visited_sources = {}
    for s in self.sources.sources:
      visited_sources[id(s)] = 0


    def trace(ds,n):

      p = integrator.y
      add = (lambda d : nodes.append(d)) if ds > 0 else (lambda d : nodes.insert(0,d))
      last_source = None
      while n <= maxn and in_bounds(bounds,p) and not line['closed'] and integrator.successful():
        p = integrator.integrate( integrator.t + ds )
        add( {'p' : p } )
        n += 1

        if vabs(p - pend) < ds/2: # if we get within 1/2 of a step, consider the path closed.
          line['closed'] = True
          return n


        # check to see if we are in a souce halo
        current_source = None
        for s in self.sources.sources:
          if s.dist(p) < halo:
            current_source = s
        if not current_source is None: # Yes, we are
          if not 'sources' in line:
            line['sources'] = []
          line['sources'].append(current_source)

          if isinstance( current_source, tuple(TerminalSources) ):
            p = current_source.r
            add( {'p' : p } )
            return n
            

          elif isinstance( current_source, tuple(PassThroughSources) ):
            # check that we aren't just in the same one from last time
            if id(current_source) != id(last_source): # No, we aren't
              # check to see if we have been here before
              if visited_sources[id(current_source)] > 0: # Yes, we have
                return n
              else: # No, we haven't
                visited_sources[id(current_source)] += 1

                # special cases

                # if the source is a dipole, and we are going into it
                # we want to jump to the other side.
                if isinstance(s, Dipole):
                  rp = p - s.r
                  if sc.dot( rp, f(0,p) ) < 0: # traveling into the source
                    # jump to the other side.
                    #    ^  .  if we are here...
                    #    |
                    # ---*------
                    #    |
                    #       .  ...then we want to jump to here
                    #
                    rp -= 2*sc.dot( rp, s.p )*s.p/vabs(s.p)
                    p = s.r + rp
          else:
            return n

        last_source = current_source


      return n


    integrator = ig.ode( f )
    integrator.set_integrator('lsoda')
    integrator.set_initial_value(self.p0,0)

    n = 0
    pend = self.p0
    n += trace(ds,n)

    if not line['closed'] and n < maxn:
      # we probably ran out of bounds
      # try to go backwards
      pend = nodes[-1]['p']
      integrator.set_initial_value(self.p0,0)
      # if the start point was inside of a source halo, then it
      # will be counted in the visited sources list and the trace
      # function will exit immediatly instead of trying to trace
      # through the source backward. we need to decrease the visited count
      # so we can to through it if needed.
      s = self.sources.get_nearest_source(self.p0)
      if s.dist(self.p0) < halo and visited_sources[id(s)] > 0:
        visited_sources[id(s)] -= 1

      n += trace(-ds,n)



    if self.backward:
      line['nodes'].reverse()

    if len(line['nodes']) < minn:
      c = config.copy()
      c['ds'] = ds/2
      line = self.get_line(bounds,c)

    return line

class EquipotentialLine(Line):
  '''Class that calculates equipotential lines.'''

  def __init__(self, sources, p_start):
    super(EquipotentialLine,self).__init__()

    self.sources = sources
    self.p0 = sc.array(p_start)

  def get_line(self, bounds=None, config={}):
    '''Get a set of points that lie on an equipotential line passing through the point p.'''
    ds   = config.get('ds',1e-1)
    maxn = config.get('maxn',1000)

    line = { 'closed' : False, 'nodes' : [] }
    nodes = line['nodes']

    # equipotential lines are perpindicular to field lines (and therefore the field vectors)
    # so the "direction field" for equipotential lines is just the vector field rotated by 90 degrees.

    # f is the function that returns the tangent vector for the line at all points in space
    f = lambda s, p: vrot( (vnorm(self.sources.E(p))), pi/2 )
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



TerminalSources = []
PassThroughSources = []

class Source(object):
  '''Class represnting a field source.'''
  def __init__(self,**kwargs):
    pass
  def E(self,r):
    '''Returns the electric field due to the source at a given point.'''
    return sc.array([0,0])
  def B(self,r):
    '''Returns the magnetic field due to the source at a given point.'''
    return sc.array([0,0])
  def V(self,r):
    '''Returns the electric potential (with respect to ground at infinity) due to the source at a given point.'''
    return 0
  def pos(self):
    return self.r
  def dist(self,r):
    return float('inf')

class Monopole(Source):
  '''A point charge source.'''
  def __init__(self, r, q=1, **kwargs):
    super(Monopole,self).__init__(**kwargs)
    self.r = sc.array(r)
    self.q = q

  def dist(self,r):
    return vmag(r-self.r)

  def E(self,r):
    Exy = sc.zeros(2)
    r = sc.array(r)
    dr = r - self.r
    d = vabs(dr)
    Exy += self.q * dr / d**3

    return Exy

  def V(self,r):
    Vxy = 0
    r = sc.array(r)
    dr = r - self.r
    d = vabs(dr)
    Vxy += self.q / d

    return Vxy

TerminalSources.append(Monopole)

# implement method to draw point charges
def _make_Monopole_drawing(self, source, scale):
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

FieldPlotDocument._make_Monopole_drawing = _make_Monopole_drawing

class PointCharge(Monopole):
  pass

class Dipole(Source):
  '''A dipole charge source.'''
  def __init__(self, r, p, **kwargs ):
    if 'halo' not in kwargs:
      kwargs['halo'] = 1e-1
    super(Dipole,self).__init__(**kwargs)
    self.r = sc.array(r,dtype=float)
    self.p = sc.array(p,dtype=float)

  def dist(self,r):
    return vmag(r-self.r)

  def E(self,r):
    Exy = sc.zeros(2)
    r = sc.array(r,dtype=float)
    dr = r - self.r
    d = vabs(dr)
    p = self.p
    Exy += (3.*sc.dot(p,dr)*dr - p*d**2) / (4.*pi*d**5)

    return Exy

  B = E

PassThroughSources.append(Dipole)

def _make_Dipole_drawing(self, source, scale):
  '''Create an SVG element for a point charge.'''
  dwg = self.dwg

  g = dwg.g()
  c = dwg.circle(r=14, fill='gray')
  g.add( c )
  c = dwg.circle(r=14,fill='url(#glare)',stroke='black',stroke_width=2)
  g.add( c )
  p = dwg.path(fill='black')
  p.push('M 2,2 H -8 V -2 H 6 V 2 Z')
  p.push('M 3 0 L 2 6 L 10 0 L 2 -6 L 3 0 z')
  g.add(p)
  g.translate(*source.r)
  g.rotate(degrees(vdir(source.p)))
  g.scale(float(scale)/self.unit)

  return g

FieldPlotDocument._make_Dipole_drawing = _make_Dipole_drawing



class SourceCollection(object):
  '''A collection of field sources.'''
  def __init__(self):
    self.sources = []

  def add_source(self,s):
    self.sources.append(s)

  def get_nearest_source(self, r, filter = lambda s : True ):
    distance = float('inf')
    nearest_source = None
    for source in self.sources:
      if filter(source) and source.dist(r) < distance:
        nearest_source = source
        distance = source.dist(r)
    return nearest_source
      
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

  def add_lines_around_point( self, sources, p, **config ):
    '''Creates lines that will pass through a set of points along a halo around
       a given point.
    '''
    numlines    = config['numlines']
    halo        = config.get('halo', 1e-1)
    delta_theta = config.get('delta_theta', 2*pi)
    theta0      = config.get('theta0', 0)
    backward    = config.get('backward', False )

    p = sc.array(p)

    lines = []
    if (delta_theta - 2*pi) < 1e-5:
      delta_theta -= delta_theta/(numlines)

    if numlines > 1:
      dtheta = delta_theta/(numlines-1)
    else:
      dtheta = 0
      delta_theta = 0

    for i in range(numlines):
      theta = (theta0 - delta_theta/2) + i*dtheta
      dr = halo*sc.array( [cos(theta), sin(theta)] )
      lines.append( FieldLine( sources, p + dr, backward ) )

    self.lines += lines

  def add_lines_around_sources( self, sources, **config ):
    filter = config.get('filter', lambda s : True)
    for s in sources.sources:
      if filter(s):
        self.add_lines_around_point( sources, s.pos(), **config )


  def add_lines_to_point_charges( self, sources, **config ):
    '''Creates lines that will leave (or enter) all point sources at uniformly spaced angles.
       For example, to make sure that at least 10 lines leaving (or enter) each point charge in
       the drawing.
    '''
    config['filter'] = lambda s : isinstance(s,Monopole)
    for s in sources.sources:
      self.add_lines_around_point( sources, s.pos(), **config )



  def add_lines_to_dipole_charges( self, sources, numlines, theta = pi/2, halo=1e-1 ):
    config['filter'] = lambda s : isinstance(s,Dipole)
    for s in sources.sources:
      self.add_lines_around_point( sources, s.pos(), **config )

class EquipotentialLineCollection(LineCollection):
  '''A collection of equipotential lines with methods to help create them.'''
  def __init__(self,*args,**kwargs):
    super(EquipotentialLineCollection,self).__init__(*args,**kwargs)

  def add_lines_along_line(self,sources,l,**config):
    numlines = config['numlines']
    halo     = config.get('halo', 1e-1)
    equal_spacing = config.get('equal_spacing', False)

    pstart = sc.array(l[0],dtype='float')
    pend   = sc.array(l[1],dtype='float')

    rhat = vnorm(pend-pstart)

    pstart += halo*rhat
    pend   -= halo*rhat

    L = vabs(pend-pstart)

    if equal_spacing:
      dl = L/(numlines-1)
      points = [ pstart + i*dl*rhat for i in range(numlines) ]
    else:
      points = []
      Vstart = sources.V(pstart)
      Vend   = sources.V(pend)
      dV = (Vend-Vstart)/(numlines-1)
      points.append( pstart )
      points.append( pend )

      for i in range(1,numlines-1):
        # V = 
        V = Vstart + i*dV
        f = lambda s : sources.V( pstart + rhat*s ) - V
        s = op.brentq( f, 0, L )
        points.append( pstart + s*rhat )


    lines = []
    for p in points:
      lines.append( EquipotentialLine( sources, p) )

    self.lines += lines



    pass
