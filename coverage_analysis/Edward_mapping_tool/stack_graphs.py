#!/usr/bin/env/python3

import matplotlib.pyplot as plt
from typing import List
import sys
from PIL import Image

title = sys.argv[1]
graph_title = sys.argv[2]

def stack_graph(imageList: List[str], reduceFactor: int):
    images = [Image.open(x) for x in imageList]
    widths, heights = zip(*(i.size for i in images))
    newWidth = max(widths) // reduceFactor
    newHeight = sum(heights) // reduceFactor
    new_im = Image.new('RGB', (newWidth, newHeight))
    y_offset = 0
    for im in images:
        w, h = im.size
        w //= reduceFactor
        h //= reduceFactor
        new_im.paste(im.resize((w, h)), (0, y_offset))
        y_offset += h  
    new_im.save(title + '/stacked_graphs.png')

#generate a title
fig = plt.figure(figsize=(20,1))
fig.suptitle('\n' + graph_title, fontsize=20)
plt.savefig(title + '/graph_title.png', dpi=600)
    
imageList = [title + '/graph_title.png']

with open('jobfile_' + title + '.txt', 'r') as file:
    for line in file:
        imageList.append(title + "/" + line[:-1] + ".png")

reduceFactor = 4

stack_graph(imageList, reduceFactor)
