import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from cygnus.process_xvg import PMF, Histo
sns.set_style("ticks")
sns.set_palette("colorblind")

wham_profile = PMF('profile.xvg')
print(wham_profile.xvg_file)
wham_profile.process_data()
wham_profile.get_dG()

fig, ax = plt.subplots()
wham_profile.plot()
plt.show()
