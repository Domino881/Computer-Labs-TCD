import matplotlib.pyplot as plt
import numpy as np

f = open("example.txt", "r")
x = [float(x) for x in f.readline().split(' ') if not x in ['\n','']]
y = [[float(x) for x in f.readline().split(' ') if not x in ['\n','']] for i in range(8)]
# mags = [float(x) for x in f.readline().split(' ') if not x in ['\n','']]

y = [[a[i] for a in y] for i in range(30)]
y_vals = [np.mean(a) for a in y]
y_errs = [np.std(a) for a in y]
print(y)

# plt.plot(x,y_vals)
plt.errorbar(x,y_vals,y_errs)
plt.xscale('log', base=2.0)
# F, ax = plt.subplots(1,1, figsize=(10,4))
# ax[0].plot(temps, ens)
# ax[1].plot(temps, mags)

# data = [float(x) for x in f.readline().split(' ') if not x in ['\n','']]
# plt.plot(data)