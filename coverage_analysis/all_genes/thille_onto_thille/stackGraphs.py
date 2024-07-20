#!/usr/bin/env/python3

from PIL import Image
import os

# List of image paths which you want to stack
images_1_15 = ['1_15_title.png',
               'thille_coverage_output.cov_ptg000010l_smoothing30k.png', 
               'thille_coverage_output.cov_ptg000001l_smoothing30k.png', 
               'thille_coverage_output.cov_ptg000011l_smoothing30k.png', 
               'thille_coverage_output.cov_ptg000027l_smoothing30k.png',
               'thille_coverage_output.cov_ptg000033l_smoothing30k.png', 
               'thille_coverage_output.cov_ptg000043l_smoothing30k.png',
               'thille_coverage_output.cov_ptg000032l_smoothing30k.png',            
               'thille_coverage_output.cov_ptg000015l_smoothing30k.png',
               'thille_coverage_output.cov_ptg000018l_smoothing30k.png',                 
               'thille_coverage_output.cov_ptg000012l_smoothing30k.png', 
               'thille_coverage_output.cov_ptg000023l_smoothing30k.png', 
               'thille_coverage_output.cov_ptg000016l_smoothing30k.png',
               'thille_coverage_output.cov_ptg000017l_smoothing30k.png',
               'thille_coverage_output.cov_ptg000006l_smoothing30k.png',
               'thille_coverage_output.cov_ptg000054l_smoothing30k.png'
               ]

images_16_30  = ['16_30_title.png',
               'thille_coverage_output.cov_ptg000005l_smoothing30k.png', 
               'thille_coverage_output.cov_ptg000028l_smoothing30k.png', 
               'thille_coverage_output.cov_ptg000025l_smoothing30k.png', 
               'thille_coverage_output.cov_ptg000013l_smoothing30k.png',
               'thille_coverage_output.cov_ptg000007l_smoothing30k.png', 
               'thille_coverage_output.cov_ptg000021l_smoothing30k.png',
               'thille_coverage_output.cov_ptg000008l_smoothing30k.png',            
               'thille_coverage_output.cov_ptg000004l_smoothing30k.png',
               'thille_coverage_output.cov_ptg000026l_smoothing30k.png',                 
               'thille_coverage_output.cov_ptg000019l_smoothing30k.png', 
               'thille_coverage_output.cov_ptg000045l_smoothing30k.png', 
               'thille_coverage_output.cov_ptg000036l_smoothing30k.png',
               'thille_coverage_output.cov_ptg000002l_smoothing30k.png', 
               'thille_coverage_output.cov_ptg000096l_smoothing30k.png',
               'thille_coverage_output.cov_ptg000041c_smoothing30k.png'
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

    if i == 0: new_im.save('stacked_1-15.png')
    if i == 1: new_im.save('stacked_16-30.png')
