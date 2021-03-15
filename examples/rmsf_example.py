import matplotlib.pyplot as plt
import seaborn as sns
from cygnus import XvgFile

sns.set_style('ticks')
sns.set_palette('colorblind')

rmsf_data = XvgFile('./a-holo-rmsf_p885_avg.xvg')
rmsf_data.process_data()
print(rmsf_data.x_label)
print(rmsf_data.y_label)

fig, ax = plt.subplots()
rmsf_data.plot()
ax.set_ylabel('RMSF (nm)')
fig.tight_layout()
plt.show()
