import os,sys
import pandas as pd 
import matplotlib.pyplot as plt

# Import CSVs
antigen_df = pd.read_csv('./metadata/chip_atlas_antigen_list.csv')
celltype_df = pd.read_csv('./metadata/chip_atlas_celltype_list.csv')
file_df = pd.read_csv('./metadata/chip_atlas_file_list.csv')

# Graphics
colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99','#ffb3e6','#c2c2f0','ivory']
gen_dict = {'hg19': 'H. sapiens', 'mm9': 'M. musculus', 'sacCer3': 'S. cerevisiae', 'dm3':'D. melanogaster','ce10': 'C. elegans', 'rn6':'R. norvegicus'}

# Count Series
gen_counts = pd.Series.value_counts(antigen_df['Genome assembly']).rename(gen_dict)
antigen_counts = antigen_df.groupby(['Genome assembly','Antigen class']).sum().rename(gen_dict).unstack().rename_axis(None)

# Plot Antigen Class Breakdown by Genome Assembly
antigen_counts['Number of experiments'].plot(kind='barh', color=colors, 
     align='center',figsize=(6,6),stacked=True,width=0.75, edgecolor='k')
plt.title('Antigen Class by Genome Assembly in ChIP Atlas\n',loc='left')
plt.axes().legend(loc='upper right')
plt.xlabel('Number of Experiments')
plt.savefig('antigen_by_organism.pdf', format='pdf', dpi=1000, bbox_inches = 'tight')

# Genome Assembly Breakdown
plt.figure()
gen_counts.plot(kind='pie', colors = colors, autopct='%1.1f%%',
	startangle=250, wedgeprops={'linewidth':1, 'edgecolor':'black'}, figsize=(6,6))
plt.axis('equal')
plt.title('ChIP Atlas Breakdown by Host Genome Assembly\n\n',loc='left')
plt.tight_layout()
plt.ylabel('')
plt.savefig('gen_breakdown.pdf', format='pdf', dpi=1000, bbox_inches = 'tight')