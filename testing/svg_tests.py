#! /usr/bin/python

import svgwrite

doc = svgwrite.Drawing(filename='test.svg', size = ("800px", "600px"))

doc.add(doc.rect(insert = (0, 0), size = ("200px", "100px"),
                                  stroke_width = "1",
                                  stroke = "black",
                                  fill = "rgb(255,255,0)"))

doc.add(doc.text("Hello World", insert = (210, 110)))

print(doc.tostring())

doc.save()
