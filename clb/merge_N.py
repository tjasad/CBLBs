import numpy as np

def merge_N(x,y):
    x1 = np.append(x, np.zeros([x.shape[0], y.shape[1]]), axis=1)
    y1 = np.append(np.zeros([y.shape[0], x.shape[1]]), y, axis=1)
    return np.append(x1,y1,axis=0)


x = np.array([[0, 1],[0, 1]])
y = np.array([[1,2,3],[4,5,6],[7,8,9]])
print(merge_N(y,x))
