import numpy as np
import time

# Create a large array with zeros and a single non-zero element
size = 10**7
arr = np.zeros(size, dtype=int)
arr[size // 2] = 1  # Place the first non-zero element in the middle

# Using np.flatnonzero
start = time.time()
result_np = np.flatnonzero(arr)[0]  # Fetch the first non-zero element
print(f"np.flatnonzero took {time.time() - start:.5f} seconds, result: {result_np}")

# Using manual search
start = time.time()
result_manual = -1
for i, val in enumerate(arr):
    if val != 0:
        result_manual = i
        break
print(f"Manual search took {time.time() - start:.5f} seconds, result: {result_manual}")

# Validate correctness
assert result_np == result_manual
