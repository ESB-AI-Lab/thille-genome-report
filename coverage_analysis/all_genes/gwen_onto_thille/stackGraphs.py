#!/usr/bin/env/python3

from PIL import Image
import os

# List of image paths which you want to stack
images_1_15 = ['1_15_title.png',
               'ptg000010l_gwen_genes.png', 
               'ptg000001l_gwen_genes.png', 
               'ptg000011l_gwen_genes.png', 
               'ptg000027l_gwen_genes.png',
               'ptg000033l_gwen_genes.png', 
               'ptg000043l_gwen_genes.png',
               'ptg000032l_gwen_genes.png',            
               'ptg000015l_gwen_genes.png',
               'ptg000018l_gwen_genes.png',                 
               'ptg000012l_gwen_genes.png', 
               'ptg000023l_gwen_genes.png', 
               'ptg000016l_gwen_genes.png',
               'ptg000017l_gwen_genes.png',
               'ptg000006l_gwen_genes.png',
               'ptg000054l_gwen_genes.png'
               ]

images_16_30  = ['16_30_title.png',
               'ptg000005l_gwen_genes.png', 
               'ptg000028l_gwen_genes.png', 
               'ptg000025l_gwen_genes.png', 
               'ptg000013l_gwen_genes.png',
               'ptg000007l_gwen_genes.png', 
               'ptg000021l_gwen_genes.png',
               'ptg000008l_gwen_genes.png',            
               'ptg000004l_gwen_genes.png',
               'ptg000026l_gwen_genes.png',                 
               'ptg000019l_gwen_genes.png', 
               'ptg000045l_gwen_genes.png', 
               'ptg000036l_gwen_genes.png',
               'ptg000002l_gwen_genes.png', 
               'ptg000096l_gwen_genes.png',
               'ptg000041c_gwen_genes.png'
               ]

imageList = [images_1_15, images_16_30]

for i in range(2):
    images = [Image.open(x) for x in imageList[i]]

    widths, heights = zip(*(i.size for i in images))

    newWidth = max(widths)
    newHeight = sum(heights)

    new_im = Image.new('RGB', (newWidth, newHeight))

    y_offset = 0
    for im in images:
        w, h = im.size
        new_im.paste(im, (0, y_offset))
        y_offset += h

    if i == 0: new_im.save('stacked_1-15.pdf')
    if i == 1: new_im.save('stacked_16-30.pdf')
