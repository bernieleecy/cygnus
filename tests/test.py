import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import cygnus as cyg

trim33a_rmsf = cyg.RMSFData('a-holo-rmsf_p885_avg.xvg',
                            label='TRIM33A (peptide bound)')
trim33b_rmsf = cyg.RMSFData('b-rmsf-avg.xvg',
                            label='TRIM33B (peptide bound)')
old_trim33a_rmsf = cyg.RMSFData('a-holo-rmsf_avg.xvg',
                                label='TRIM33A (peptide bound, old)')

fig, ax = plt.subplots(figsize=(6,4))
print('running compare rmsf')
cyg.compare_rmsf([trim33a_rmsf, trim33b_rmsf, old_trim33a_rmsf], ax=ax)
plt.show()
