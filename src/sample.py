from vectorfieldplot import *

doc = FieldplotDocument('VFPt_cylindrical_magnet_thumb', width=320, height=165, unit=30, commons=True) 
field = Field({'monopoles':[ [-1, 0, 1] ]})
line = FieldLine(field, [-1, 0], start_v=[0, 1], directions='forward')
doc.draw_line(line)
doc.draw_charges(field)
doc.draw_object('circle', {'cx':0, 'cy':0, 'r':1})
doc.write()
