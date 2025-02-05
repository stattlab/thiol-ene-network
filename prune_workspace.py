import os

for i in os.listdir("./workspace"):
	if os.path.isdir("./workspace/" + i) and os.path.isfile("./workspace/" + i + "/alternating_q.json"):
		os.remove("./workspace/" + i + "/polymerize.gsd")
