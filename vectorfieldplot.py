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
'''
 
import sys
from math import *
import bisect

import svgwrite
from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier
import scipy as sc
import scipy.optimize as op
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
    if d != 0.: return sc.array(x) / vabs(x)
    return sc.array(x)
 
def rot(xy, phi):
    '''
    2D vector rotation
    '''
    s = sin(phi); c = cos(phi)
    return sc.array([c * xy[0] - s * xy[1], c * xy[1] + s * xy[0]])
 
def cosv(v1, v2):
    '''
    find the cosine of the angle between two vectors
    '''
    d1 = sum(v1**2); d2 = sum(v2**2)
    if d1 != 1.: d1 = sqrt(d1)
    if d2 != 1.: d2 = sqrt(d2)
    if d1 * d2 == 0.: return 1.
    return sc.dot(v1, v2) / (d1 * d2)
 
def sinv(v1, v2):
    '''
    find the sine of the angle between two vectors
    '''
    d1 = sum(v1**2); d2 = sum(v2**2)
    if d1 != 1.: d1 = sqrt(d1)
    if d2 != 1.: d2 = sqrt(d2)
    if d1 * d2 == 0.: return 0.
    return (v1[0] * v2[1] - v1[1] * v2[0]) / (d1 * d2)
 
def angle_dif(a1, a2):
    return ((a2 - a1 + pi) % (2. * pi)) - pi
 
def list_interpolate(l, t):
    n = max(0, bisect.bisect_right(l, t) - 1)
    s = None
    if t < l[0]:
        if l[1] - l[0] == 0.:
            s = 0.
        else:
            s = (t - l[0]) / float(l[1] - l[0])
    elif t >= l[-1]:
        n = len(l) - 2
        if l[-1] - l[-2] == 0.:
            s = 1.
        else:
            s = (t - l[-2]) / float(l[-1] - l[-2])
    else:
        s = (t - l[n]) / (l[n+1] - l[n])
    return n, s
 
def pretty_vec(p):
    return '{0:> 9.5f},{1:> 9.5f}'.format(p[0], p[1])
 
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

def grad( f, p ):
  '''Compute the gradient of a function.'''

  dx = 1e-5
  dy = 1e-5

  f0 = f(p)
  p[0] += dx
  fx = (f(p) - f0) / dx
  p[0] -= dx
  p[1] += dy
  fy = (f(p) - f0) / dy
  p[1] -= dy

  return sc.array([fx,fy])

def in_bounds(bounds,p):
  '''Check if a point is within a given bounding box.'''
  if bounds is None:
    return True
  return p[0] >= bounds['x0'] and p[0] <= bounds['x1'] and p[1] >= bounds['y0'] and p[1] <= bounds['y1']
 
 







class FieldplotDocument:
    '''
    Class representing an image.
    '''
    def __init__ (self, name, width=800, height=600, digits=3.5, unit=100, center=None, license='CC-BY', author="UNKNOWN"):
        self.name = name
        self.width = float(width)
        self.height = float(height)
        self.digits = float(digits)
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


        # misc info
        self.arrow_geo = {'x_nock':0.3,'x_head':3.8,'x_tail':-2.2,'width':4.5}





    def __add_arrow(self,color):
      arrow = self.dwg.path(id=color+'_arrow', stroke='none', fill=color)
      arrow.push( 'M {0},0 L {1},{3} L {2},0 L {1},-{3} L {0},0 Z'.format(
          self.arrow_geo['x_nock'], self.arrow_geo['x_tail'],
          self.arrow_geo['x_head'], self.arrow_geo['width'] / 2.))
      arrow.scale(1./self.unit)

      get_elem_by_id( self.dwg.defs, 'arrows' ).add(arrow)

    def __get_arrows_on_polylines(self, polylines, arrows_style, linewidth, linecolor='#000000', closed=False):
        '''
        draws arrows on polylines.
        options in "arrows_style":
        min_arrows: minimum number of arrows per segment
        max_arrows: maximum number of arrows per segment (None: no limit)
        dist: optimum distance between arrows
        scale: relative size of arrows to linewidth
        offsets [start_offset, mid_end, mid_start, end_offset]
        fixed_ends [True, False, False, True]:
        	make first/last arrow distance invariable
        '''
        min_arrows = 1; max_arrows = None; arrows_dist = 1.; scale = linewidth
        offsets = 4 * [0.5]; fixed_ends = 4 * [False]
        if 'min_arrows' in arrows_style:
            min_arrows = arrows_style['min_arrows']
        if 'max_arrows' in arrows_style:
            max_arrows = arrows_style['max_arrows']
        if 'dist' in arrows_style:
            arrows_dist = arrows_style['dist']
        if 'scale' in arrows_style:
            scale *= arrows_style['scale']
        if 'offsets' in arrows_style:
            offsets = arrows_style['offsets']
        if 'fixed_ends' in arrows_style:
            fixed_ends = arrows_style['fixed_ends']
 
        arrows = []
        arrows_def = get_elem_by_id( self.dwg.defs, 'arrows' )
        for j, polyline in enumerate(polylines):
            line_points = polyline['path']
            mina = min_arrows
            maxa = max_arrows
            # measure drawn path length
            lines_dist = [0.]
            for i in range(1, len(line_points)):
                lines_dist.append(lines_dist[-1]
                    + vabs(line_points[i] - line_points[i-1]))
 
            offs = [offsets[2], offsets[1]]
            fixed = [fixed_ends[2], fixed_ends[1]]
            if polyline['start']:
                offs[0] = offsets[0]
                fixed[0] = fixed_ends[0]
            if polyline['end']:
                offs[1] = offsets[3]
                fixed[1] = fixed_ends[3]
 
            d01 = [0., lines_dist[-1]]
            for i in [0, 1]:
                if fixed[i]:
                    d01[i] += offs[i] * arrows_dist * [1., -1.][i]
                    mina -= 1
                    if not maxa is None: maxa -= 1
            if d01[1] - d01[0] < 0.: break
            elif d01[1] - d01[0] == 0.: d_list = [d01[0]]
            else:
                d_list = []
                if fixed[0]: d_list.append(d01[0])
                if maxa > 0 or maxa == None:
                    number_intervals = (d01[1] - d01[0]) / arrows_dist
                    number_offsets = 0.
                    for i in [0, 1]:
                        if fixed[i]: number_offsets += .5
                        else: number_offsets += offs[i] - .5
                    n = int(number_intervals - number_offsets + 0.5)
                    n = max(n, mina)
                    if not maxa is None: n = min(n, maxa)
                    if n > 0:
                        d = (d01[1] - d01[0]) / float(n + number_offsets)
                        if fixed[0]: d_start = d01[0] + d
                        else: d_start = offs[0] * d
                        for i in range(n):
                            d_list.append(d_start + i * d)
                if fixed[1]: d_list.append(d01[1])
 
            geo = self.arrow_geo # shortcut
            #### arrow drawing ####
            for d1 in d_list:
                # calculate arrow position and direction
                if d1 < 0. or d1 > lines_dist[-1]: continue
                d0 = d1 + (geo['x_nock'] * scale + 2.5*linewidth *
                    (geo['x_tail'] - geo['x_nock']) / geo['width']) / self.unit
                if closed and d0 < 0.: d0 = lines_dist[-1] + d0
                d2 = d1 + (geo['x_head'] * scale + linewidth *
                    (geo['x_tail'] - geo['x_head']) / geo['width']) / self.unit
                if closed and d2 > lines_dist[-1]: d1 -= lines_dist[-1]
                i0, s0 = list_interpolate(lines_dist, d0)
                i1, s1 = list_interpolate(lines_dist, d1)
                i2, s2 = list_interpolate(lines_dist, d2)
                p0 = line_points[i0] + s0 * (line_points[i0+1]-line_points[i0])
                p1 = line_points[i1] + s1 * (line_points[i1+1]-line_points[i1])
                p2 = line_points[i2] + s2 * (line_points[i2+1]-line_points[i2])
                p = None; angle = None
                if vabs(p2-p1) <= .1**self.digits or (d2 <= d0 and not closed):
                    v = line_points[i1+1] - line_points[i1]
                    p = p1
                    angle = atan2(v[1], v[0])
                else:
                    v = p2 - p0
                    p = p0 + sc.dot(p1 - p0, v) * v / vabs(v)**2
                    angle = atan2(v[1], v[0])
 
                arrow = self.dwg.use( get_elem_by_id( arrows_def, linecolor+"_arrow" ) )
                arrow.translate( p[0], p[1] )
                arrow.rotate(degrees(angle))
                arrow.scale(scale)

                arrows.append(arrow)

        return arrows

    def __get_bounds(self):
        bounds = {}
        bounds['x0'] = -self.center[0] / self.unit
        bounds['y0'] = -(self.height - self.center[1]) / self.unit
        bounds['x1'] = (self.width - self.center[0]) / self.unit
        bounds['y1'] =  self.center[1] / self.unit

        return bounds


    def __make_pointcharge(self, source, id, scale):
      dwg = self.dwg

      g = dwg.g(id=id)
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

      dwg = self.dwg
      img = self.img


      container = get_elem_by_id(self.img,'sources')
      for i,source in enumerate(sources.sources):
        if stypes == 'All' or isinstance( source, stypes ):
          g = dwg.g(id="None")
          if isinstance( source, PointCharge ):
            g = self.__make_pointcharge( source, 'charge{0}'.format(i), scale )
          container.add( g )


    def __draw_line(self,line,linewidth,linecolor):
      bounds = self.__get_bounds()
      line = line.get_line(bounds)
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


      container = get_elem_by_id(self.img,'equipotentiallines')
      line = self.dwg.path( path.d(), stroke=linecolor,stroke_width=linewidth/self.unit,fill='none' )
      container.add( line )

    def draw_fieldline(self, line, linewidth=2, linecolor='black'):
      self.__draw_line( line, linewidth, linecolor )

    def draw_equipotentialline(self, line, linewidth=2, linecolor='red'):
      self.__draw_line( line, linewidth, linecolor )

    def write(self, filename=None):
      self.dwg.save()
      print 'image written to', self.dwg.filename
























class FieldVectors:
  '''
  calculates field vectors
  '''
  def __init__(self, field, minscale=0, maxscale=1, unit=None, halo=None):
    self.minscale = minscale
    self.maxscale = maxscale
    self.unit = unit
    self.halo = halo
    self.field = field
    self.vectors = []
    self.scaled_vectors = []

  def add_vectors_in_grid(self, xmin, ymin, xmax, ymax, xN, yN):
    xmin = 1.*xmin
    xmax = 1.*xmax
    dx = (xmax - xmin) / (xN - 1)

    ymin = 1.*ymin
    ymax = 1.*ymax
    dy = (ymax - ymin) / (yN - 1)

    if self.unit is None:
      self.unit = min(dx,dy)

    XYs = []
    for x in sc.linspace( xmin, xmax, num=xN ):
      for y in sc.linspace( ymin, ymax, num=yN ):
        XYs.append( (x,y) )

    self.__add_vectors(XYs)
    self.__create_scaled_vectors( )

  def add_vectors_on_line(self, line, ds):
    '''Adds vectors'''

  def __add_vectors(self, XYs):
    for x,y in XYs:
      f = self.field.F([x,y])
      self.vectors.append( (x,y,f) )

  def __create_scaled_vectors(self, trans = None ):
    if trans is None:
      trans = self.__log_scale
    def inhalo(p):
      if self.halo is None:
        return False

      pd  = self.__get_distance_to_nearest_pole( sc.array(p) )
      halo = self.halo
      if self.halo == 'auto':
        halo = self.maxscale*self.unit

      return  pd and pd < halo

    self.scaled_vectors = []
    minlen = maxlen = None
    for i in range(len(self.vectors)):
      x,y,f = self.vectors[i]
      if inhalo([x,y]):
        continue

      if minlen is None or vabs(f) < minlen:
        minlen = vabs(f)
      if maxlen is None or vabs(f) > maxlen:
        maxlen = vabs(f)

    self.scaled_vectors = []
    for i in range(len(self.vectors)):
      x,y,f = self.vectors[i]
      if inhalo([x,y]):
        continue
      self.scaled_vectors.append( (x,y,trans( f, self.minscale*self.unit, self.maxscale*self.unit, minlen, maxlen )) )

  def __log_scale(self, f, newmin, newmax, oldmin, oldmax):
    oldlen = vabs(f)
    ratio = log10( oldlen / oldmin ) / log10( oldmax / oldmin )
    newlen = newmin + (newmax - newmin) * ratio

    fx = newlen * f[0] / oldlen
    fy = newlen * f[1] / oldlen

    return [fx,fy]

  def __get_distance_to_nearest_pole(self, p, v=None):
      '''
      returns distance to nearest pole
      '''
      d_near = None
      for ptype, poles in self.field.elements.iteritems():
          if ptype not in ['monopoles', 'dipoles'] or len(poles) == 0:
              continue
          for pole in poles:
              d = vabs(pole[:2] - p)
              if d_near is None or d < d_near:
                d_near = d

      return d_near



class FieldLine:
  '''calculates field lines.'''

  def __init__(self, sources, p_start):

    self.sources = sources
    self.p0 = sc.array(p_start)

  def get_line(self, bounds=None, ds=1e-1, maxn=1000, ftype='E',backward=False):
    '''Get a set of points that lie on a field line passing through the point p.'''

    line = { 'closed' : False, 'nodes' : [] }
    nodes = line['nodes']

    # f is the function that returns the tangent vector for the line at all points in space
    if ftype == 'E':
      if backward:
        f = lambda s, p: vnorm(self.sources.E(p))
      else:
        f = lambda s, p: -vnorm(self.sources.E(p))
    elif ftype == 'B':
      if backward:
        f = lambda s, p: vnorm(self.sources.B(p))
      else:
        f = lambda s, p: -vnorm(self.sources.B(p))
    else:
      raise BaseException("ERROR: Unknown ftype "+ftype)

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


    return line


class EquipotentialLine:
  '''calculates equipotential lines.'''

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
      print out,p,self.p0,vabs(p-self.p0)
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


    




class Source(object):
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

  
class SourceCollection:
  def __init__(self):
    self.sources = []

    self.pole_sources = ( PointCharge)

  def add_source(self,s):
    self.sources.append(s)

  def get_nearest_pole(self,r):
    ns = None
    nd = None
    for s in self.sources:
      if isinstance( s, self.pole_sources ):
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







 
