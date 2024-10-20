import numpy as np
import scipy
from itertools import combinations, chain
import time
import random

# Predefined cost vector (Paper: Polynomial calculation of the Shapley value based on sampling)
# a = np.array(list(range(1, 11)))
# s = np.array([8, 12, 6, 14, 8, 9, 13, 10, 10, 10])
# seed random number generator
random.seed(1)
a = [random.randint(1, 10) for _ in range(10)]
s = np.array([2, 2, 1, 3, 2, 2, 3, 2, 2, 2])
cost = np.repeat(a, s)
print(cost)

def airport_game(coalition):
    coalition = list(coalition)
    return np.max(cost[coalition])

def powerset(iterable):
    "powerset([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))

start = time.time()

NP = 20
player=2
num = np.setdiff1d(np.array(list(range(1,NP+1))), [player])
u = list()
for ind_S_no_i in powerset(num):
    u.append(airport_game(np.append(ind_S_no_i, player)) - airport_game(ind_S_no_i))

end = time.time()
print(end - start)

print(np.max(u))