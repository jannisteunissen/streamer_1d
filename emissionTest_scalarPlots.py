import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


f = pd.read_csv("output/emissionLEATest/emissionTest_fluid_scalars.txt", delim_whitespace=True)


plt.plot(f["time"], f["sigmapos_l"], "r.")
plt.plot(f["time"], f["sigmaneg_l"], "b.")
plt.legend(("Positive", "Negative"))
plt.show()