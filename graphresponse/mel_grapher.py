""" Graphing module """
import graphresponse.subsetdata as sub
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = "monospace"
from matplotlib.lines import Line2D

specie = 'Melanogaster'


repA, repB = sub.create_dictionary(species=specie)


plt.rcParams.update({
    "lines.color": "white",
    "patch.edgecolor": "white",
    "text.color": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": "dimgray",
    "axes.labelcolor": "white",
    "legend.fontsize": 18,
    "xtick.color": "white",
    "ytick.color": "white",
    "grid.color": "black",
    "figure.facecolor": "black",
    "figure.edgecolor": "black",
    "savefig.facecolor": "black",
    "savefig.edgecolor": "black"})


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))
ax.set_title('Drosophila melanogaster', fontsize='24', style='italic')
for key in repA.keys():
    if key != 'control':
        if 'up' in key:
            ax.scatter(sub.get_gens(repA[key]), sub.get_resp(repA[key]), marker=',', color='darkorchid')
            ax.plot(sub.get_gens(repA[key]), sub.get_resp(repA[key]), color='darkorchid', linewidth=4.0)
            ax.errorbar(sub.get_gens(repA[key]), sub.get_resp(repA[key]), yerr=sub.make_error(repA[key], specie),
                        color='darkorchid', linewidth=4.0)
            ax.scatter(sub.get_gens(repB[key]), sub.get_resp(repB[key]), alpha=0.4, marker='o', color='darkorchid')
            ax.plot(sub.get_gens(repB[key]), sub.get_resp(repB[key]), alpha=0.4, linestyle=(0, (5, 1)),
                    color='darkorchid', linewidth=4.0)
            ax.errorbar(sub.get_gens(repB[key]), sub.get_resp(repB[key]), yerr=sub.make_error(repB[key], specie),
                        color='darkorchid', alpha=0.4, linewidth=4.0)
        else:
            ax.scatter(sub.get_gens(repA[key]), sub.get_resp(repA[key]), marker=',', color='darkgreen')
            ax.plot(sub.get_gens(repA[key]), sub.get_resp(repA[key]), marker=',', color='darkgreen', linewidth=4.0)
            ax.errorbar(sub.get_gens(repA[key]), sub.get_resp(repA[key]), yerr=sub.make_error(repA[key], specie),
                        color='darkgreen', linewidth=4.0)
            ax.scatter(sub.get_gens(repB[key]), sub.get_resp(repB[key]), alpha=0.4, marker='o', color='darkgreen')
            ax.plot(sub.get_gens(repB[key]), sub.get_resp(repB[key]), alpha=0.4, linestyle=(0, (5, 1)),
                    color='darkgreen', linewidth=4.0)
            ax.errorbar(sub.get_gens(repB[key]), sub.get_resp(repB[key]), yerr=sub.make_error(repB[key], specie),
                        color='darkgreen', alpha=0.4, linewidth=4.0)

custom_lines = [Line2D([0], [0], color='darkorchid', lw=4, label='Up A'),
                Line2D([0], [0], color='darkorchid', lw=4, alpha=0.4, label='Up B'),
                Line2D([0], [0], color='darkgreen', lw=4, label='Down A'),
                Line2D([0], [0], color='darkgreen', lw=4, alpha=0.4, label='Down B')]
legen = ax.legend(handles=custom_lines, loc='upper left')
plt.setp(legen.get_texts(), color='k')
xticks = list(range(0, 20, 2))
xtickslabels = [str(x) for x in xticks]
ax.set_xticks(xticks)
ax.set_xticklabels(xtickslabels)
plt.setp(ax.get_xticklabels(), fontsize=18)
plt.setp(ax.get_yticklabels(), fontsize=18)
# horitzontal line at Generation 15
ax.axvline(x=15, color='black', linestyle=(0, (3, 10, 1, 10, 1, 10)))
ax.set_ylabel('Phenotypic Change in Standard Deviation', fontsize=18)
ax.set_xlabel('Generation', fontsize=18)
# fig.show()
# plt.tight_layout()
fig.savefig('Melanogaster2018_PhenotypicResponse.png', bbox_inches='tight')
