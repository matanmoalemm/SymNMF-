import numpy as np
import pandas as pd
import symnmfmodule as C
import sys
np.random.seed(1234)

def read_data(file_path):
    file = open(file_path)
    data = file.readlines()
    file.close()
    for i in range(len(data)):
        data[i] = data[i].split(",")
    if data[0] == '\n':
        return 0
    for x in data:
        for i in range(len(x)):
            x[i] = float(x[i])

    return data

def sym(data):
    n = len(data)
    A = np.zeros((n,n))
    for i in range(n):
        x_i = np.array(data[i])
        for j in range(n):
            if i != j:
                x_j = np.array(data[j])
                A[i][j] = np.exp(-(np.linalg.norm(x_i - x_j))**2/2)
    return A

def norm(A):
    n = A.shape[0]
    D = np.zeros((n,n))
    sums = np.sum(A, axis = 0)
    for i in range(n):
        D[i][i] = 1/np.sqrt(sums[i])
    W = np.dot(D,A)
    W = np.dot(W,D)
    return W


def initialize_H(W,k):
    n = W.shape[0]
    m = np.average(W)
    H = np.zeros((n,k))
    val = 2*np.sqrt(m/k)
    for i in range(n):
        for j in range(k):
            H[i][j] = np.random.uniform(low = 0 , high = val)
    return H

def symnmf(data,k):
        W = norm(sym(data))
        H = initialize_H(W,k)
        H = H.tolist()
        W = W.tolist()
        finalH = C.symnmf_c(H,W,k)
        return finalH

def printPoints(list):
    for x in list:
        print(",".join(f"{coord:.4f}" for coord in x))

def main(args):

    data = read_data(args[2])
    k = int(args[0])
    if k > len(data):
        print("An Error Has Occurred")
        return None
    if args[1] == "symnmf":
        printPoints(symnmf(data,k))
    elif args[1] == 'sym':
        A = C.sym_c(data)
        printPoints(A)
    elif args[1] == 'ddg':
        D = C.ddg_c(data)
        printPoints(D)
    elif args[1] == 'norm':
        W = C.norm_c(data)
        printPoints(W)
    else:
        print("An Error Has Occurred")
        return None

if __name__ == "__main__":

    main(sys.argv[1:])
