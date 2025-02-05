import os

direc = "../workspace/"
for j in os.listdir(direc):
    polymerize_file = direc + j + "/polymerize.gsd"
    if not os.path.isfile(polymerize_file):
        continue
    if os.path.isfile("./" + j + ".txt"):
        os.system("rm " + polymerize_file)
        continue
    polymerize1_file = direc + j + "/polymerize1.gsd"
    print(j)
    os.system("./bin/sample_graph -i " + polymerize_file + " -o " + j + ".txt")
    os.system("python trim_gsd.py -i " + polymerize_file + " -o polymerize1.gsd" + " -p 10")
    print("cut")
    if os.path.isfile("./" + j + ".txt"):
    	os.system("rm " + polymerize_file)
    	print("removed")
    #os.system("rm " + polymerize_file)
