import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from math import ceil

plt.rcParams["figure.figsize"] = (5,3)

def get_bm(t):
    w = np.random.default_rng().normal(0.0, size=t)
    return np.cumsum(w)

def get_ou(t, alpha, theta, beta):
    w = np.random.default_rng().normal(0.0, size=t)
    x = [np.random.default_rng().normal(0.0, size=1)]
    for i in range(t):
        dx = alpha * (theta - x[i]) + beta*w[i]
        x.append(x[i]+dx)
    return x

t = 100
n = 23

# bm = np.asarray([get_bm(t) for i in range(n)])
# plt.figure()
# for i in range(9, n):
#     plt.plot(bm[i], c="black")
# for i in range(4, 9):
#     plt.plot(bm[i], c="blue")
# for i in range(4):
#     plt.plot(bm[i], c="red")
# plt.ylim(-20, 20)
# plt.title("Neutral")
# plt.savefig("figures/models/bm.png")
# # plt.show()s

ou = np.asarray([get_ou(t, 0.05, 10, .95) for i in range(n)])
plt.figure()
for i in range(9, n):
    plt.plot(ou[i], c="black")
for i in range(4, 9):
    plt.plot(ou[i], c="blue")
for i in range(4):
    plt.plot(ou[i], c="red")
plt.ylim(-20, 20)
plt.title("Constrained")
plt.savefig("figures/models/constrained.png")
# plt.show()

ou1 = np.asarray([get_ou(t, 0.05, -10, .95) for i in range(9)])
ou2 = np.asarray([get_ou(t, 0.05, 10, .95) for i in range(14)])
plt.figure()
for i in range(4):
    plt.plot(ou1[i], c="red")
for i in range(4, 9):
    plt.plot(ou1[i], c="blue")
for i in range(0, 14):
    plt.plot(ou2[i], c="black")
plt.ylim(-20, 20)
plt.title("Adaptive - Clade")
plt.savefig("figures/models/clade.png")
# plt.show()

ou1 = np.asarray([get_ou(t, 0.05, -10, .95) for i in range(4)])
ou2 = np.asarray([get_ou(t, 0.05, 10, .95) for i in range(19)])
plt.figure()
for i in range(4):
    plt.plot(ou1[i], c="red")
for i in range(5, 19):
    plt.plot(ou2[i], c="black")
for i in range(0, 5):
    plt.plot(ou2[i], c="blue")
plt.ylim(-20, 20)
plt.title("Adaptive - Aggressive")
plt.savefig("figures/models/agg.png")