import os,sys
import pandas as pd 
import matplotlib.pyplot as plt

# Initializes dataframes necessary for metadata analysis: antigen_df, celltype_df, 
#	file_df (type pandas dataframe)
# Initializes color scheme for plots: colors (type list)
# Initializes dictionary mapping genome assemblies to organism: gen_dict (type dict)
#
# Returns: antigen_df, celltype_df, file_df, colors, gen_dict
def init_dfs():
	colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99','#ffb3e6','#c2c2f0','ivory']
	gen_dict = {'hg19': 'H. sapiens', 'mm9': 'M. musculus', 'sacCer3': 'S. cerevisiae', 'dm3':'D. melanogaster','ce10': 'C. elegans', 'rn6':'R. norvegicus'}
	return pd.read_csv('./metadata/chip_atlas_antigen_list.csv'), pd.read_csv('./metadata/chip_atlas_celltype_list.csv'), pd.read_csv('./metadata/chip_atlas_file_list.csv'), colors, gen_dict

# Plots available data for cell class breakdown and assay antigens for given target
# genome assembly
#
# Returns: None (saves plots as .pdfs)
def plot_gen_metastats(target='hg19'):
	antigen_df, celltype_df, file_df, colors, gen_dict = init_dfs()

	# Class Breakdown
	plt.figure()
	celltype_df[celltype_df['Genome assembly'] == target].groupby('Cell type class').sum().rename_axis(None).plot(kind='barh', 
	  	color=colors, align='center',figsize=(6,6),stacked=True, width=0.75, edgecolor='k')
	plt.title('Cell Class Type Breakdown for '+gen_dict[target]+'\n',loc='left')
	plt.xlabel('Number of Experiments')
	plt.ylabel('')
	plt.legend().set_visible(False)
	plt.savefig(target+'_class_breakdown.pdf', format='pdf', dpi=1000, bbox_inches = 'tight')

	# Antigen Breakdown
	plt.figure()
	antigen_counts = antigen_df.groupby(['Genome assembly','Antigen class']).sum().rename(gen_dict).unstack().rename_axis(None)
	antigen_counts.loc[gen_dict[target]]['Number of experiments'].plot(kind='pie', colors = colors, autopct='%1.1f%%',
	    startangle=110, wedgeprops={'linewidth':1, 'edgecolor':'black'}, figsize=(6,6))
	plt.axis('equal')
	plt.title('Antigen Class Breakdown for '+gen_dict[target]+'\n',loc='left')
	plt.ylabel('')
	plt.savefig(target+'_ant_breakdown.pdf', format='pdf', dpi=1000, bbox_inches = 'tight')

# Returns list of experiment IDs for input params
def get_experiment_list(target='hg19',cell_class='Blood',cell_type='Lymphoblastoid cell line'):
	antigen_df, celltype_df, file_df, colors, gen_dict = init_dfs()
	ctype_df = celltype_df[(celltype_df['Genome assembly'] == target) & (celltype_df['Cell type class'] == cell_class) & (celltype_df['Cell type'] == cell_type)].drop(['Genome assembly','Cell type class','Cell type'],axis=1)
	print(int(ctype_df['Number of experiments'])," experiments found for: \n",gen_dict[target],' --> ',cell_class,' --> ',cell_type)
	return pd.Series.tolist(ctype_df['Experimental IDs included'])[0].split(',')

#plot_gen_metastats()
#experiments = get_experiment_list()
#print(experiments)


