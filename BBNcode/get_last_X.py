elist = [
    ["n"],
    ["H_1", "p", "^1H"],
    ["H_2", "d", "D"],
    ["He_3", "^3He"],
    ["H_3", "t", "T"],
    ["He_4", "^4He"],
    ["Be_7", "^7Be"],
    ["Li_7", "^7Li"],
    ["Li_6", "^6Li"],
    ["He_6", "^6He"]
]

import pickle
import os

if not os.path.isfile("smart_cache.pickle"):
    exit()


with open("smart_cache.pickle", "rb") as f:
    i, X_ans, Tres, ode_params = pickle.load(f)

for j in range(len(X_ans[-1])):
    print("{}:\t{}".format(elist[j][0], X_ans[-1][j]))