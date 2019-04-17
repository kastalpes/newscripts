#!/usr/bin/env python

help_string="""
 usage: % pager.py [options] files
 Makes web pages to easily view multiple (timesteps, fields, simulations) quickly.
 Run this in a directory with a number of png files, 
 it will produce oot.html.
 If no *files* arguments are passed, all png files in the current working directory are used.

 Expects all files to be formatted
 simname_0004_quantity.png
 or
 simname_n0004_quantity.png
 and will sort timestep (n0004, etc) into rows, "quantity" into columns, and 
 "simname" into small squares at each time/quantity position.
 
 For yt, make your save look like
 pw.save("%s_%04d"%(sim,frame))
 with no type suffix, and the "quantity" will be produced in a way that makes sense.

 Files in the list or current working directory that don't match are stored in 
 "skipped.txt"


 """

row_order_help = \
"""Column/row order. 
 1: (row, column, inner) = frame, field, sim
 2: (row, column, inner) = sim, field, frame
 """
# >>> python pager.py -h 
# for more options


import sys
import glob
import re
import pdb
import numpy as np

from optparse import OptionParser
parser = OptionParser(help_string)
parser.add_option("-w", "--width", dest="width",action="store",help = "width", default="300")
parser.add_option("-t", "--title", dest="title",action="store",help = "title", default=None)
parser.add_option("-n", "--name", dest="name",action="store",help = "fame", default="oot.html")
parser.add_option("-c", "--number_of_columns", dest="number_of_columns", help = "number of inner columns", default=2)
parser.add_option("-o", "--column_order", dest = "column_order", help = row_order_help, default=1)
parser.add_option("-k", "--caption_file", dest="caption_file",help="space separated list of <run><caption>", default='captions.txt')
parser.add_option("-z", "--zoom_sequence", action='store_true', dest="zoom_sequence",help="If this is a sequence of zooms, the names are parsed differently", default=False)
parser.add_option("-d", "--directory_trim", action='store_true',dest='directory_trim', help="If sims are in sub directories, and there's a one-to-one directory-sim relation, remove direcory names", default=False)
parser.add_option("-x", "--xtra", action='store_true',dest='xtra_skipped', help="Put links to extra files that dont match the ax19_0019_projectin.png format", default=True)
options, args = parser.parse_args()
#title=options.title
width = options.width

#this_fname_temp = 'p33_%s_%04d_2d-Profile_%s_%s_cell_mass.png'
#filename_template = r'([^_]*)_(\d\d\d\d)_([^_]*)_%s_(.*)_cell_mass.png'
#filename_template = r'([^_]*)_(\d\d\d\d)_(.*).png'
#filename_template = r'(.*)_(\d\d\d\d)_(.*).png' #pretty good version
if options.zoom_sequence:
    filename_template = r'(.*)_n{0,1}(\d\d\d\d)_(zoom\d+){0,1}_(.*).png' #with the zoom.
else:
    filename_template = r'(.*)_[TD]{0,1}D{0,1}n{0,1}(\d\d\d\d)_(.*).png' #pretty good version
this_fname_temp = '%s_%04d_%s.png'
framelist = []
fieldlist = []
simlist = []
framere = re.compile(filename_template)
files_skipped = []

max_zoom = -1
name_dict={}
if len(args) > 0:
    fnames = args
else:
    fnames = glob.glob("*png")
for fname in fnames:
    mmm= framere.match(fname)
    if mmm:
        sim = mmm.group(1)
        if options.directory_trim:
            sim = sim.split("/")[-1]
        if sim not in simlist:
            simlist.append(sim)
            name_dict[sim]={}
        frame = int(mmm.group(2))
        field_or_zoom = mmm.group(3)
        if field_or_zoom.startswith('zoom'):
            field = mmm.group(4)
            zoom = int(mmm.group(3)[5:])
            max_zoom = max([max_zoom,zoom])
        else:
            field = mmm.group(3)
            zoom = -1
        if frame not in framelist:
            framelist.append(frame)
        if field not in fieldlist:
            fieldlist.append(field)

        if not frame in name_dict[sim]:
            name_dict[sim][frame] = {}

        if not field in name_dict[sim][frame]:
            name_dict[sim][frame][field] = {}

        #if not name_dict[sim][frame].has_key(zoom):
        #    name_dict[sim][frame][field][zoom] = {}

        if zoom in name_dict[sim][frame][field]:
            print( "ERROR name collision (keeping first)"%( name_dict[sim][frame][field][zoom], fname))
        else:
            name_dict[sim][frame][field][zoom] = fname

    else:
        files_skipped.append(fname)

fptr_skipped = open("skipped.txt","w")
for fname in files_skipped:
    fptr_skipped.write("%s\n"%fname)
fptr_skipped.close()

framelist = sorted(framelist)
simlist = sorted(simlist)
fieldlist = sorted(fieldlist)

caption = {}
if len(glob.glob( options.caption_file ) ):
    fptr = open(options.caption_file,'r')
    for line in fptr:
        spl = line.split(" ")
        if len(spl) > 1:
            caption[spl[0]] = line[ line.index(" "):].strip()
    fptr.close()
else:
    fptr = open(options.caption_file,'w')
    for sim in simlist:
        fptr.write("%s ---\n"%sim)
    print( "wrote a caption file", options.caption_file)
    fptr.close()




print( simlist)
print(framelist)
print(fieldlist)
print("zoom", max_zoom)
title = "%s"*len(simlist)%tuple(simlist)
if options.title is not None:
    title = options.title
fptr = open(options.name,'w')
fptr.write('<html>\n')
fptr.write('<title>%s</title>'%title)

if options.xtra_skipped:
    still_skipped=[]
    fptr.write('<table boerder="1">\n')
    fptr.write('<tr>')
    for nf, fname in enumerate(files_skipped):
        if fname.split('.')[-1] in ['png','jpg']:
            img_tag = '<td><a href="%s"><img src="%s" width='+width+'></a></td>'
            fptr.write(img_tag%(fname,fname))
        else:
            still_skipped.append(fname)
    fptr.write('</tr>')
    if len(still_skipped):
        #haha, this does nothing because we are only reading png in the first place.
        fptr.write('<tr>')
        for nf, fname in enumerate(still_skipped):
            fptr.write("<td><a href='%s'>%s</a></td>"%(fname,fname))

        fptr.write('</tr>')

    fptr.write('</table>\n')
        



fptr.write('<table border="1">\n')
class label_tool():
    def __init__(self,order):
        self.order = int(order)
        if self.order == 1:
            self.outer_list = [-1]+framelist
            self.second_list = fieldlist
            self.inner_list = simlist
            self.id={'fr':0,'fi':1,'r':2} 
            self.inner_label_template = "n%04d %s<br>" #take outer, second value
            self.inner_caption_template = "%s (%s)"
        if self.order == 2:
            self.outer_list = ['']+simlist
            self.second_list = fieldlist
            self.inner_list = framelist
            self.id={'fr':2,'fi':1,'r':0} 
            self.inner_label_template = "%s %s<br>" #take outer, second value
            self.inner_caption_template = "%s (%s)"
        if self.order == 3:
            self.outer_list = ['']+simlist
            self.second_list = framelist
            self.inner_list = fieldlist
            self.id={'fr':1,'fi':2,'r':0} 
            self.inner_label_template = "%s %d<br>" #take outer, second value
            self.inner_caption_template = "%s (%s)"
        if self.order == 4:
            self.outer_list = [-1]+framelist
            self.second_list = simlist
            self.inner_list = fieldlist
            self.id={'fr':0,'fi':2,'r':1} 
            self.inner_label_template = "n%04d %s<br>" #take outer, second value
            self.inner_caption_template = "%s (%s)"
        if self.order == 5:
            self.outer_list = ['']+fieldlist
            self.second_list = simlist
            self.inner_list = framelist
            self.id={'fr':2,'fi':0,'r':1} 
            self.inner_label_template = "%s %s<br>" #take outer, second value
            self.inner_caption_template = "%d (%s)"
    def set_values(self,*args):
        self.frame=args[self.id['fr']]
        self.field=args[self.id['fi']]
        self.run=args[self.id['r']]
    def inner_label(self,outer_value,second_value):
        return self.inner_label_template%(outer_value, second_value)
    def inner_caption(self,inner_value):
        return self.inner_caption_template%(inner_value, caption.get(inner_value,"---"))
LT = label_tool(options.column_order)


for nouter,outer_value in enumerate(LT.outer_list):
    fptr.write('<tr>')
    fptr.write('<td class="td_frame_number"> %s </td>'%str(outer_value))
    for second_value in LT.second_list:
        fptr.write('<td class="td_figure_table">')
        if nouter ==0:
            fptr.write("%s"%(second_value))
        else:
            for zoom in range(-1,max_zoom+1):
                img_tag = '<a h<figure><a href="%s"><img src="%s" width='+width+'></a><figcaption>%s</figcaption></figure>'
                fptr.write(LT.inner_label(outer_value,second_value))

                fptr.write('<table border="2"><tr>\n')
                for n,inner_value in enumerate(LT.inner_list):
                    LT.set_values(outer_value,second_value,inner_value)
                    run = LT.run; frame = LT.frame; field=LT.field
                    
                    fptr.write('<td class="td_image">')
                    if frame in name_dict[run] and field in name_dict[run][frame] and zoom in  name_dict[run][frame][field]:
                        this_fname = name_dict[run][frame][field].pop(zoom)
                        fptr.write(img_tag%(this_fname,this_fname,LT.inner_caption(inner_value)))
                    else:
                        this_fname = this_fname_temp%(run,frame,field)
                        fptr.write("%s<br>"%this_fname)
                    fptr.write("</td>")
                    if (n+1)%int(options.number_of_columns) == 0:
                        fptr.write("</tr><tr>\n")
                
                fptr.write('</tr></table><br>\n')

        fptr.write('</td>\n')
    fptr.write('</tr>\n')

fptr.write('</table>\n')
fptr.write('</html>\n')
fptr.close()

#Make sure we got everyting.
for sim in name_dict.keys():
    for frame in name_dict[sim].keys():
        for field in name_dict[sim][frame].keys():
            if len(name_dict[sim][frame][field].keys()) > 0:
                print( "PARSE ERROR: did not properly treat", name_dict[sim][frame])


#p33_ai01_0025_2d-Profile_density_HeI_Density_cell_mass.png 
#end
