import pickle
import os

if not os.path.isfile("smart_cache.pickle"):
    exit()


with open("smart_cache.pickle", "rb") as f:
    i, X_ans, Tres, ode_params, elements = pickle.load(f)

for j in range(len(X_ans[-1])):
    print("{}:\t{}".format(elements[j][0], X_ans[-1][j]))