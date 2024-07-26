import numpy as np

arr = np.array([[1, 2, 3], [4, np.nan, 6]], dtype=np.float32)

cdef float[:,:] arr_cy = arr

print("numpy")
print(arr)
print("cython")
for i in range(arr.shape[0]):
    for j in range(arr.shape[1]):
        v = arr_cy[i, j]
        print(v)
