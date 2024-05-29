# convert chiapet "bed" files downloaded from ucsc database to bed then index by tabix
# remember to use ucsc bedSort before zipping !!!

import sys, os

if len(sys.argv) != 3:
	print('Usage: {0} <input> <output full name>'.format(sys.argv[0]))
	sys.exit()

def topos(string):
	t = string.split(':')
	t2 = t[1].split('-')
	return (t[0], int(t2[0]), int(t2[1]))


infile,outfile = sys.argv[1:]

fout = open(outfile, 'w')
n = 1

with open(infile) as fin:

	#### skip header
	next(fin)

	for line in fin:

		lst=line.rstrip().split('\t')
	
		#### one anchor
		v1 = lst[0:3]
		v11 = lst[3] +  ':' + lst[4] +  '-' + lst[5]
		
		#### second anchor
		v2 = lst[3:6]
		v22 = lst[0] + ':' + lst[1] +  '-' + lst[2]

		#### x should be number of contacts. 
		x = lst[9]

		#### formatted output
		fout.write('{0[0]}\t{0[1]}\t{0[2]}\t{1},{2}\t{3}\t{4}\n'.format(v1, v11, x, n, '+'))
		n += 1
		fout.write('{0[0]}\t{0[1]}\t{0[2]}\t{1},{2}\t{3}\t{4}\n'.format(v2, v22, x, n, '-'))
		n += 1

fout.close()

#os.system('bedSort {0} {0}'.format(outfile))
#os.system('bgzip {0}'.format(outfile))
#os.system('tabix -p bed {0}.gz'.format(outfile))

