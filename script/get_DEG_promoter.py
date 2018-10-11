import sys
import re

promoter = {}

with open(sys.argv[1]) as f:
	for line in f:
		lst = line.strip().split('\t')
		promoter[lst[0]] = lst[1]

with open(sys.argv[2]) as f:
	for line in f:
		lst = line.strip().split('\t')
		if lst[0] != 'Gene':
			lst[0] = re.sub('\..+', '', lst[0])
			if lst[0] in promoter:
				#pass
				print ">%s\n%s" % (lst[0], promoter[lst[0]])
			else:
				pass
				#print lst[0]
