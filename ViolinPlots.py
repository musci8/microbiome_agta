import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import seaborn as sns
import matplotlib.cm as cm

sns.set_style('white')
sns.set_context('talk')
mpl.rcParams.update({'text.usetex': True})
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()


colors = ['#1b9e77','#d95f02','#7570b3','#e7298a']

df_violin = pd.read_csv('RatioDyadsSocialASVs.csv')
motifs = df_violin.motif.unique()
motifs1 = ['different\nhouseholds','same\nhousehold']
motifs2 = ['Friends\nsame camp','Friends\ndifferent camps']
motifs3 = ['different camps','same camp']
motifs4 = list(set(motifs).difference(motifs1).difference(motifs2).difference(motifs3))
Ms = [motifs1,motifs2,motifs3,motifs4]


trimmed = []

for m in motifs:
    a = df_violin.query('motif==@m')
    t = np.percentile(a.score,99)
    trimmed.append(a.query('score<@t'))
    
df_violin_trimmed = pd.concat(trimmed,ignore_index=True)


fig3 = plt.figure(constrained_layout=True,figsize=(10,8))
gs = fig3.add_gridspec(3, 2)
ax1 = fig3.add_subplot(gs[:, 1])
A = df_violin_trimmed.query('motif in @motifs4')

sns.violinplot(A.score,y=A.motif,scale='width',cut=0,ax=ax1,width=.8,color=colors[0])
ax1.axvline(x=1,color='red',ls='--',lw=3)

ax1.set_xlim([1e-1,40])
ax1.set_ylabel('')
ax2 = fig3.add_subplot(gs[0, 0])
A = df_violin_trimmed.query('motif in @motifs1')

sns.violinplot(A.score,y=A.motif,scale='width',cut=0,ax=ax2,width=.8,color=colors[1])
ax2.axvline(x=1,color='red',ls='--',lw=3)

ax2.set_xlim([1e-1,40])
ax2.set_ylabel('')
ax2.set_xticks([])
ax2.set_xlabel('')
ax3 = fig3.add_subplot(gs[1, 0])
A = df_violin_trimmed.query('motif in @motifs2')

sns.violinplot(A.score,y=A.motif,scale='width',cut=0,ax=ax3,width=.8,color=colors[2])
ax3.axvline(x=1,color='red',ls='--',lw=3)

ax3.set_xlim([1e-1,40])
ax3.set_ylabel('')
ax3.set_xticks([])
ax3.set_xlabel('')
ax4 = fig3.add_subplot(gs[2, 0])
A = df_violin_trimmed.query('motif in @motifs3')

sns.violinplot(A.score,y=A.motif,scale='width',cut=0,ax=ax4,width=.8,color=colors[3])
ax4.axvline(x=1,color='red',ls='--',lw=3)

ax4.set_xlim([1e-1,40])
ax4.set_ylabel('')

plt.tight_layout()
plt.savefig('ViolinPlotsSocialOtus.png')

df_violin_ns = pd.read_csv('RatioDyadsOtherASVs.csv')


fig3 = plt.figure(constrained_layout=True,figsize=(10,8))
gs = fig3.add_gridspec(3, 2)
ax1 = fig3.add_subplot(gs[:, 1])
A = df_violin_ns.query('motif in @motifs4')

sns.violinplot(A.score,y=A.motif,scale='width',cut=0,ax=ax1,width=.8,color=colors[0])
ax1.axvline(x=1,color='red',ls='--',lw=3)

ax1.set_xlim([1e-1,40])
ax1.set_ylabel('')
ax2 = fig3.add_subplot(gs[0, 0])
A = df_violin_ns.query('motif in @motifs1')

sns.violinplot(A.score,y=A.motif,scale='width',cut=0,ax=ax2,width=.8,color=colors[1])
ax2.axvline(x=1,color='red',ls='--',lw=3)

ax2.set_xlim([1e-1,40])
ax2.set_ylabel('')
ax2.set_xticks([])
ax2.set_xlabel('')
ax3 = fig3.add_subplot(gs[1, 0])
A = df_violin_ns.query('motif in @motifs2')

sns.violinplot(A.score,y=A.motif,scale='width',cut=0,ax=ax3,width=.8,color=colors[2])
ax3.axvline(x=1,color='red',ls='--',lw=3)

ax3.set_xlim([1e-1,40])
ax3.set_ylabel('')
ax3.set_xticks([])
ax3.set_xlabel('')
ax4 = fig3.add_subplot(gs[2, 0])
A = df_violin_ns.query('motif in @motifs3')

sns.violinplot(A.score,y=A.motif,scale='width',cut=0,ax=ax4,width=.8,color=colors[3])
ax4.axvline(x=1,color='red',ls='--',lw=3)

ax4.set_xlim([1e-1,40])
ax4.set_ylabel('')

plt.tight_layout()
plt.savefig('ViolinPlotsOtherOtus.png')
