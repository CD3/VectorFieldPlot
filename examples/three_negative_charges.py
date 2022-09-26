from vectorfieldplot import *

"""From the example https://commons.wikimedia.org/wiki/File:VFPt_charges_plus_minus_thumb.svg"""

doc = FieldplotDocument('three_negative_charges', width=220, height=165, unit=30., commons=True)
field = Field({'monopoles':[[0.5,-2,-1], [2,0,-1], [0.5,2,-1]]})
doc.draw_charges(field, scale=11./14.)
n = 60
dy = (2. - -2)/(n-1)
for i in range(n):
    line = FieldLine(field, [-2,-2 + i*dy], start_v=[1, 0], directions='forward')
    doc.draw_line(line, linewidth=1.0, arrows_style={
        'min_arrows':1, 'offsets':[.65,.5,.5,.65],
        'scale':1.4})
doc.write()
