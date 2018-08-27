import sys

counter = 0
with open(sys.argv[1]) as f:
	for line in f:
		s=line[:-1]
		id=s.split("\t")[0]
		if id == "id":
			continue
		counter+=1
		if counter % 10000 == 0:
			print(id)
		o = open(sys.argv[2]+'/'+str(id), 'a+');
		o.write(s + '\n');
		o.close()
