import os

direc = "../workspace/"
for j in os.listdir(direc):
    polymerize_file = direc + j + "/polymerize.gsd"
    if os.path.isfile(polymerize_file):
        print(j)
