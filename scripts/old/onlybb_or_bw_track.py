######
#!/usr/bin/python3
#####
#!/home/rcf-40/amalthom/panases_soft/anaconda3/bin/python3.4
import sys, os, re,glob

if len(sys.argv)<4:
	print("python3 <script> <configfile> <chiporatacorbigbed> <bborbw> <count(default=0)>\nExiting!!!!\n")
	exit(0)
header = '''\
##################################
# Tracks                  #
##################################
'''

supertrack_template = '''\
track {group_name}
superTrack on
group {group_name}
shortLabel {group_name}
longLabel {group_name}
visibility full

'''

pileup_template = '''\
  track {track_name}-pileup 
  parent {parent}
  priority {priority}
  longLabel {label}-pileup
  shortLabel {label}-pileup
  type bigWig
  autoScale off
  viewLimits 0:3
  gridDefault on
  color {color}
  visibility full
  bigDataUrl {filename}

'''

peak_template = '''\
  track {track_name}
  parent {parent}
  visibility dense
  priority {priority}
  longLabel {label}
  shortLabel {label}
  type bigBed
  color {color}
  bigDataUrl {filename}

'''
count=0
if len(sys.argv)==5:
	count=int(sys.argv[4].strip())
priorities = {}
priority = count 

#chip_tf_supertrack = supertrack_template.format(group_name = 'chip-tf')
#chip_chrom_supertrack = supertrack_template.format(group_name = 'chip-chrom')
#chip_input_supertrack = supertrack_template.format(group_name = 'chip-input')
#rna_supertrack = supertrack_template.format(group_name = 'rna')
atac_supertrack = supertrack_template.format(group_name = sys.argv[2])

infile=open(sys.argv[1])
line_colors = {}
bw_names =[]
if sys.argv[3]=="bw":
	rename=".bw"
if sys.argv[3]=="bb":
	rename=".bb"
for line in infile:
	if line !="" and line !="\n":
		elements=line.rstrip().split("\t")
		line_colors[elements[0][:-3]]=elements[1]
		bw_names.append(elements[0])
infile.close()


bw_base=sys.argv[2]+"-pileups/"
bb_base=sys.argv[2]+"-peaks/"

for filename in bw_names:
	if filename!="" and filename!='\n':
		name=filename[:-3]
		#print(name)
		#line = name[6:]
		#print(line)
		color = line_colors[name]
		#print(color)
		priority += 1
		
		if name not in priorities:
			priorities[name] = priority
		if sys.argv[3]=="bb":
			atac_supertrack += peak_template.format(track_name = name, parent = sys.argv[2], priority = priorities[name],\
												label = name, color = color, filename = bb_base+name+".bb")
		if sys.argv[3]== "bw":
			atac_supertrack += pileup_template.format(track_name = name, parent = sys.argv[2], priority = priorities[name]+0.5, \
													label = name, color = color, filename = bw_base+name+".bw")

outfile=open(sys.argv[3]+"_trackDb.txt",'w')
outfile.write(header+"\n"+atac_supertrack)
outfile.close()

#print(chip_tf_supertrack)
#print(chip_chrom_supertrack)
#print(chip_input_supertrack)
#print(rna_supertrack)


