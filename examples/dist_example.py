import matplotlib.pyplot as plt
import seaborn as sns
from cygnus import XvgFile

sns.set_style("ticks")
sns.set_palette("colorblind")

dist = XvgFile('./gmx_dist_output.xvg')
dist.process_data()

fig, ax = plt.subplots()
sns.histplot(x=dist.y_data*10)
fig.tight_layout()
plt.show()
