import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv("./swap_time t1.csv")

# Plot
plt.figure(figsize=(10, 6))
plt.plot(data['MatrixSize'], np.log10(data['RowSwapTime']), label='Row Swap ($log_{10}(t)$)', marker='o')
plt.plot(data['MatrixSize'], np.log10(data['ColumnSwapTime']), label='Column Swap ($log_{10}(t)$)', marker='s')

plt.xscale("log")
plt.xticks(data["MatrixSize"], labels=data["MatrixSize"])
plt.xlabel("Matrix Size", fontsize=12)
plt.ylabel("Time $log_{10}(t)$", fontsize=12)
plt.title("Matrix Swap Performance", fontsize=12)
plt.legend(fontsize=12)
plt.grid(True, which="both")

# Save plot
plt.savefig("swapPerformance .png", dpi=300)
plt.show()