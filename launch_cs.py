import argparse
import random
import numpy as np
from itertools import combinations
import cvxpy as cp
import ast


def bernoulli_matrix(m, n):
    arr = np.random.random([m, n])
    sp = np.random.binomial(1, arr)
    sp = np.float64(sp)
    sp[sp == 1] = 1./np.sqrt(m)
    sp[sp == 0] = -1./np.sqrt(m)
    return sp

def random_combination(iterable, r):
    "Random selection from itertools.combinations(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(range(n), r))
    return tuple(pool[i] for i in indices)

def run_cs(v, args):
    filename = f"nodes/indices_{args.g}.txt"
    with open(filename) as f:
        lines = f.readlines()
        ind_EV_loc = np.sort(ast.literal_eval(lines[0]))

    m = np.int32(args.m) #100
    t = np.int32(args.t) #50
    n = len(ind_EV_loc)
    eps = 0.001

    np.random.seed(19365)

    requests = []

    combs = list()
    for i in range(len(ind_EV_loc)):
        combs.append([p for p in combinations(ind_EV_loc, i+1)])

    m_combs = list()
    for c in combs:
        for p in c:
            m_combs.append(list(p))

    def utility(coalition, m_combs):
        if len(coalition) == 0:
            return 0.
        else:
            coalition = np.sort(coalition)
            requests.append(np.array(coalition))
            idx = np.where(list(map(lambda x: np.array_equal(coalition, x), m_combs)))[0][0]
            return v[idx]
    
    A = bernoulli_matrix(m, n)
    y = {}

    for cur_t in range(t):
        row_idx = np.array(random_combination(ind_EV_loc, np.random.randint(1, n + 1)))
        
        phi_arr = []
        
        for i in ind_EV_loc:
            permutation = row_idx.copy()

            if i in permutation:
                permutation_w_i = permutation.copy()
                permutation_wo_i = permutation[permutation != i]
            elif i not in permutation:
                permutation_wo_i = permutation.copy()
                permutation_w_i = np.append(permutation, i)

            u1 = utility(permutation_w_i, m_combs)
            u2 = utility(permutation_wo_i, m_combs)
            phi = u1 - u2

            phi_arr = np.append(phi_arr, phi)
        
        y_m = []
        for j in range(m):
            y_sum = 0.0
            for i in range(n):
                y_sum += A[j, i]*phi_arr[i]
            y_m = np.append(y_m, y_sum)
        
        y[cur_t] = y_m

    y_bar = np.sum(np.array(list(y.values())).T, axis=1) * (1/t)
    s_bar = utility(ind_EV_loc, m_combs) * (1/n)

    # Create variable.
    x_l1 = cp.Variable(shape=(n,1))

    # Create constraint.
    constraints = [cp.norm(A @ (s_bar + x_l1) - y_bar[:,np.newaxis]) <= eps]

    # Form objective.
    obj = cp.Minimize(cp.norm(x_l1, 1))

    # Form and solve problem.
    prob = cp.Problem(obj, constraints)
    prob.solve()

    CS_res = s_bar + x_l1.value

    result = []
    for arr in requests:
        flag = False
        for c in result:
            if np.array_equal(arr, c):
                flag = True
            
        if not flag:
            result.append(arr)

    print(f"Used size:{len(result)}")
    print(f"Total size:{len(v)}")

    return CS_res


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-i', '--input',
        help='Path to the file (.txt is accepted)',
        required=True,
    )
    parser.add_argument(
        '--g',
        help='Grid id',
        required=True,
    )
    parser.add_argument(
        '--c',
        help='Case id',
        required=True,
    )
    parser.add_argument(
        '--m',
        help='Number of CS measurements',
        default=100,
        required=False,
    )
    parser.add_argument(
        '--t',
        help='Number of CS iterations',
        default=50,
        required=False,
    )

    args = parser.parse_args()

    with open(args.input) as f:
        lines = f.readlines()
        v = ast.literal_eval(lines[0])
        CS_res = run_cs(v, args)
    
    res = np.array(CS_res.T[0])
    np.set_printoptions(linewidth=np.inf)
    np.set_printoptions(precision=4)
    print(str(f"CS_res:{res}"))