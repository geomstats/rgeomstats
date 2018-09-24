
repl_python()
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
os.environ['QT_QPA_PLATFORM_PLUGIN_PATH'] = 'C:/Users/lundj/Documents/Anaconda3/Library/plugins/platforms'

t = np.arange(0.0, 2.0, 0.01)
s = 1 + np.sin(2 * np.pi * t)

fig, ax = plt.subplots()
ax.plot(t, s)

ax.set(xlabel='time (s)', ylabel='voltage (mV)',
       title='About as simple as it gets, folks')
ax.grid()

plt.show()
