import os

nonequi = []
for sd in os.listdir("./workspace/"):
	if os.path.isfile("./workspace/" + sd + "/equi.gsd"):
		continue
	try:
		os.system("python project.py run -o equilibrate -n 1 -j " + sd)
	except:
		nonequi.append(sd)
print(nonequi)	
