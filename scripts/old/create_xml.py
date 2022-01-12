import sys, os, re,glob

if len(sys.argv)!=3:
	print("python3 <script> <xmllist(file)> <replace string:e.g .sorted.chr.nodup.filt.bam\nExiting!!!!\n")
	exit(0)

xmlformat = '''\
<ReplicateCorrelation>
	<BamFiles>
{bamsection}
	</BamFiles>
	<BamNames>
{namesection}
	</BamNames>
</ReplicateCorrelation>

'''
inputfile=open(sys.argv[1])
replace_string=sys.argv[2]

for line in inputfile:
	if line!="\n" or line!="\n":
		elements=line.rstrip().split(',')
		number=len(elements)
		bamfile=bamnames=""
		for i in range(0,number):
			bamfile=bamfile+"\t\t<File.R"+str(i+1)+">"+elements[i]+"</File.R"+str(i+1)+">\n"
			newname=elements[i].replace(replace_string,'')
			bamnames=bamnames+"\t\t<File.R"+str(i+1)+">"+newname+"</File.R"+str(i+1)+">\n"
			if i==0:
				header=newname
		#print(bamfile,bamnames,sep="\n")
		outfile=open("correlationxml_"+header+".txt",'w')
		#print(xmlformat.format(bamsection=bamfile,namesection=bamnames))
		outfile.write(xmlformat.format(bamsection=bamfile,namesection=bamnames))
		outfile.close()
