from vectorfieldplot import *

"""From the example https://commons.wikimedia.org/wiki/File:VFPt_charges_plus_minus_thumb.svg"""

doc = FieldplotDocument('VFPt_charges_plus_minus_thumb', width=220, height=165, unit=30., commons=True)
field = Field({'monopoles':[[-1,0,1], [1,0,-1]]})
doc.draw_charges(field, scale=11./14.)
n = 16
for i in range(n):
    a = 2.0 * pi * (0.5 + i) / n
    if abs(i - (n-1)/2.) > 6: fix = 4 * [False]
    else: fix = [True, False, False, True]
    line = FieldLine(field, [-1,0], start_v=[cos(a), sin(a)],
        directions='forward')
    doc.draw_line(line, linewidth=1.0, arrows_style={
        'dist':1.6, 'min_arrows':1, 'offsets':[.65,.5,.5,.65],
        'fixed_ends':fix, 'scale':1.4})
doc.write()
