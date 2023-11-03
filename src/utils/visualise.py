import numpy as np
import scanpy as sc
import scanpy.external as sce
sc.settings.set_figure_params(dpi=80)

from ..configs import config

#Plotting
import matplotlib.pyplot as plt
import seaborn as sns


def densityQCs(adataObj, hue=None):   
    #Plot them in line so they take up less space
    fig, ax = plt.subplots(1, 4, figsize=(20,5))
    fig.tight_layout(pad=2)   #space between plots
    
    if hue != None:
        hue_s = adata.obs[hue].astype('string')
    else:
        hue_s = None

    ### Genes ---------------
    d1 = sns.kdeplot(np.log10(adataObj.obs['n_genes_by_counts']), fill=True, color='cornflowerblue', hue=hue_s, ax=ax[0])
    min_x, max_x = d1.get_xlim() 

    #Threshold lines and fill
    if config.PROTO.QC.MIN_GENES != None:
        d1.axvline(np.log10(config.PROTO.QC.MIN_GENES), 0, 1, c='red')  #set manually for chosen threshold
        d1.axvspan(min_x, np.log10(config.PROTO.QC.MIN_GENES), alpha=0.2, color='red')
    if config.PROTO.QC.MAX_GENES != None:
        d1.axvline(np.log10(config.PROTO.QC.MAX_GENES), c='red')
        d1.axvspan(np.log10(config.PROTO.QC.MAX_GENES), max_x, alpha=0.2, color='red')

    ### UMI ---------------
    d2 = sns.kdeplot(np.log10(adataObj.obs['total_counts']), fill=True, color='forestgreen', hue=hue_s, ax=ax[1])
    min_x, max_x = d2.get_xlim() 
        
    if config.PROTO.QC.MIN_COUNTS != None:
        d2.axvline(np.log10(config.PROTO.QC.MIN_COUNTS), 0, 1, c='red')  #set manually for chosen threshold
        d2.axvspan(min_x, np.log10(config.PROTO.QC.MIN_COUNTS), alpha=0.2, color='red')
    if config.PROTO.QC.MAX_COUNTS != None:
        d2.axvline(np.log10(config.PROTO.QC.MAX_COUNTS), c='red')
        d2.axvspan(np.log10(config.PROTO.QC.MAX_COUNTS), max_x, alpha=0.2, color='red')

    ### Mito % ---------------
    d3 = sns.kdeplot(adataObj.obs['pct_counts_mito'], fill=True, color='coral', hue=hue_s, ax=ax[2])
    min_x, max_x = d3.get_xlim() 

    #Threshold lines and fill
    if config.PROTO.QC.MT_PERCENTAGE != None:
        d3.axvline(config.PROTO.QC.MT_PERCENTAGE, 0, 1, c='red')  #set manually for chosen threshold
        d3.axvspan(config.PROTO.QC.MT_PERCENTAGE, max_x, alpha=0.2, color='red')


    ### Ribo % ---------------
    d4 = sns.kdeplot(adataObj.obs['pct_counts_ribo'], fill=True, color='orchid', hue=hue_s, ax=ax[3])
    min_x, max_x = d4.get_xlim() 
    #ax[3].legend(loc='center left', bbox_to_anchor=(1.0, 1.0)) #upper right
    
    #Threshold lines and fill
    if config.PROTO.QC.MT_PERCENTAGE != None:
        d4.axvline(config.PROTO.QC.RB_PERCENTAGE, 0, 1, c='red')  #set manually for chosen threshold
        d4.axvspan(config.PROTO.QC.RB_PERCENTAGE, max_x, alpha=0.2, color='red')
    
    #Remove additional legends at need
    if hue != None:
        ax[0].get_legend().remove()
        ax[1].get_legend().remove()
        ax[2].get_legend().remove()
        
    # Remove all borders
    sns.despine(bottom = False, left = True)


def plot_bar(ax, x):
    ax.bar(x.index.astype(str), x.values)
    ax.set_xticks(
        ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha='right')


def volcano_plot(
    adata, group_key, group_name="clusters", groupby="time",
    title=None):
    cell_type = "_".join(group_key.split("_")[1:])
    result = sc.get.rank_genes_groups_df(adata, group=cell_type, key=group_key).copy()
    result["-logQ"] = -np.log(result["pvals"].astype("float"))
    lowqval_de = result.loc[abs(result["logfoldchanges"]) > config.PROTO.GLM.DEG_THRESHOLD.LFC]
    other_de = result.loc[abs(result["logfoldchanges"]) <= config.PROTO.GLM.DEG_THRESHOLD.LFC]

    fig, ax = plt.subplots(figsize=(7, 5))
    sns.regplot(
        x=other_de["logfoldchanges"],
        y=other_de["-logQ"],
        fit_reg=False,
        scatter_kws={"s": 6},
    )
    sns.regplot(
        x=lowqval_de["logfoldchanges"],
        y=lowqval_de["-logQ"],
        fit_reg=False,
        scatter_kws={"s": 6},
    )
    ax.set_xlabel("log2 FC")
    ax.set_ylabel("-log Q-value")

    if title is None:
        title = group_key.replace("_", " ")
    plt.title(title)
    plt.show()


def plot_heatmap(
    adata, group_key, group_name="cluster", groupby="isHuman"):
    cell_type = "_".join(group_key.split("_")[1:])
    res = sc.get.rank_genes_groups_df(adata, group=cell_type, key=group_key)
    res.index = res["names"].values
    res = res[
        (res["pvals_adj"] < config.PROTO.GLM.DEG_THRESHOLD.FDR) & (abs(res["logfoldchanges"]) > config.PROTO.GLM.DEG_THRESHOLD.LFC)
    ].sort_values(by=["logfoldchanges"])
    print(f"Plotting {len(res)} genes...")
    markers = list(res.index)
    sc.pl.heatmap(
        adata,#[adata.obs[group_name] == cell_type].copy(),
        markers,
        groupby=[groupby, 'line'],
        swap_axes=True,
    )
    return markers