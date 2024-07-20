#!/usr/bin/env/python3

from PIL import Image
import os

# List of image paths which you want to stack
images_1_12 = ['stacked_images_title.png',
               'CM056809.1_gwen_genes.png', 
               'CM056810.1_gwen_genes.png', 
               'CM056811.1_gwen_genes.png', 
               'CM056812.1_gwen_genes.png',
               'CM056813.1_gwen_genes.png', 
               'CM056814.1_gwen_genes.png',
               'CM056815.1_gwen_genes.png',            
               'CM056816.1_gwen_genes.png',
               'CM056817.1_gwen_genes.png',                 
               'CM056818.1_gwen_genes.png', 
               'CM056819.1_gwen_genes.png', 
               'CM056820.1_gwen_genes.png'
               ]



images = [Image.open(x) for x in images_1_12]

widths, heights = zip(*(i.size for i in images))

newWidth = max(widths)
newHeight = sum(heights)

new_im = Image.new('RGB', (newWidth, newHeight))

y_offset = 0
for im in images:
    w, h = im.size
    new_im.paste(im, (0, y_offset))
    y_offset += h

new_im.save('stacked_images.pdf')
