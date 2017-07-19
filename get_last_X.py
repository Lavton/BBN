import pickle
import os

if not os.path.isfile("smart_cache.pickle"):
    exit()


with open("smart_cache.pickle", "rb") as f:
    i, X_ans, Tres, ode_params, element_stuct = pickle.load(f)

for j in range(len(X_ans[-1])):
    print("{}:\t{}".format(element_stuct[j][0][0], X_ans[-1][j]))