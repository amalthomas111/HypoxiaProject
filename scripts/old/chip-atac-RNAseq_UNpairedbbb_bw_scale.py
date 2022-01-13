######
#/usr/bin/python3
#####
#!/home/rcf-40/amalthom/panases_soft/anaconda3/bin/python3.4
import sys, os, re,glob

if len(sys.argv)!=3:
	print("python3 <script> <configfile> <chip/atac/RNA-seq>\nExiting!!!!\n")
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
  shortLabel {label}
  type bigWig
  autoScale off
  maxHeightPixels 24:24:24
  viewLimits 0:20
  yLineOnOff on
  gridDefault on
  color {color}
  visibility full
  bigDataUrl {filename}
  
'''

peak_template = '''\
  track {track_name}-peak
  parent {parent}
  visibility dense
  priority {priority}
  longLabel {label}-peak
  shortLabel {label}
  type bigBed
  color {color} 
  bigDataUrl {filename}

'''

priorities = {}
priority = 0 

#chip_tf_supertrack = supertrack_template.format(group_name = 'chip-tf')
#chip_chrom_supertrack = supertrack_template.format(group_name = 'chip-chrom')
#chip_input_supertrack = supertrack_template.format(group_name = 'chip-input')
#rna_supertrack = supertrack_template.format(group_name = 'rna')
atac_supertrack = supertrack_template.format(group_name = sys.argv[2])

#line_colors = {}
bw_names =[]
strip=".bw"
#print(sys.argv[2])
if sys.argv[2]=="atac" or sys.argv[2] == "chip":
	color="197,172,24"
if sys.argv[2]=="RNAseq":
	color="128,128,128"
infile=open(sys.argv[1])
for line in infile:
	if line !="" and line !="\n":
		elements=line.rstrip()
		#line_colors[elements.strip(strip)]=color
		bw_names.append(elements)
infile.close()

#bw_base=sys.argv[2]+"-pileups/"
#bb_base=sys.argv[2]+"-peaks/"

bw_base=""
bb_base=""

count=0
for filename in bw_names:
	if filename!="" and filename!='\n':
		if filename.find(".bw")!=-1:
			track="bigwig"
		elif filename.find(".bb")!=-1:
			track="bb"
		name=filename[:-3]
		#print(name)
		#line = name[6:]
		#print(line)
		color = color
		#print(color)
		priority += 1
		
		#if name not in priorities:
		#	priorities[name] = priority
		if track=="bb":
			atac_supertrack += peak_template.format(track_name = name, parent = sys.argv[2], priority = priority,\
											label = name, color = color, filename = bb_base+name+".bb")
		if track=="bigwig":
			atac_supertrack += pileup_template.format(track_name = name, parent = sys.argv[2], priority = priority, \
												label = name, color = color, filename = bw_base+name+".bw")

outfile=open("trackDb_scale.txt",'w')
outfile.write(header+"\n"+atac_supertrack)
outfile.close()

#print(chip_tf_supertrack)
#print(chip_chrom_supertrack)
#print(chip_input_supertrack)
#print(rna_supertrack)

