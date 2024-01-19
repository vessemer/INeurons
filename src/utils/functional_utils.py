import matplotlib.pyplot as plt
import seaborn.objects as so
import subprocess

import pandas as pd
import numpy as np
import scanpy as sc
import decoupler as dc
from pathlib import Path
import gseapy as gp

from ..configs import config
from ..utils import utils as us


def load_GO():
    if config.PROTO.ENRICHMENT.REACTOME == 'C2':
        print('Initialising C2 reactome...')
        if not Path(
            config.PATHS.CSV/"c2.reactome.v7.5.1.symbols.gmt").is_file():
            subprocess.run([
                "wget", "-O", "symbols.gmt",
                "https://figshare.com/ndownloader/files/35233771" ])
            shutil.move(
                'symbols.gmt', config.PATHS.CSV/'c2.reactome.v7.5.1.symbols.gmt')

        reactome = us.gmt_to_decoupler(
            config.PATHS.CSV/"c2.reactome.v7.5.1.symbols.gmt")

    elif config.PROTO.ENRICHMENT.REACTOME == 'MSigDB':
        print('Initialising MSigDB reactome...')
        msigdb = dc.get_resource("MSigDB")
        # Get reactome pathways
        collection =    config.PROTO.ENRICHMENT.COLLECTION
        reactome = msigdb.query("collection == @collection")
        # Filter duplicates
        reactome = reactome[~reactome.duplicated(("geneset", "genesymbol"))]

    # Filtering genesets to match behaviour of fgsea
    geneset_size = reactome.groupby("geneset").size()
    gsea_genesets = geneset_size.index[ (geneset_size > 5) ]
    #& (geneset_size < 500)]
    reactome = reactome[reactome["geneset"].isin(gsea_genesets)]

    return reactome


def load_collectri():
    # Retrieve CollecTRI gene regulatory network
    collectri = dc.get_collectri(
        organism='human', split_complexes=False)
    collectri.head()
    return collectri


def run_GSEA(selected, reactome):
    print('Launched GSEA')
    scores = selected[['logfoldchanges']]
    scores, norm, pvals = dc.run_gsea(
        scores.T,
        reactome,
        source="geneset",
        target="genesymbol",
    )

    gsea_results = (
        pd.concat({"score": scores.T, "norm": norm.T, "pval": pvals.T}, axis=1)
        .droplevel(level=1, axis=1)
        .sort_values("pval")
    )

    so.Plot(
        data=(
            gsea_results.head(20).assign(
                **{"-log10(pval)": lambda x: -np.log10(x["pval"])}
            )
        ),
        x="-log10(pval)",
        y="source",
    ).add(so.Bar()).plot()
    plt.show()

    return gsea_results


def run_ULM(selected, collectri, volcano=None):
    print('Launched ULM for TF slope t-stat')
    mat = selected[['logfoldchanges']].T.rename(
        index={'logfoldchanges': 'isHuman'})
    tf_acts, tf_pvals = dc.run_ulm(mat=mat, net=collectri)
    ' '.join(tf_pvals.T[tf_pvals.T.isHuman < .01].index)

    dc.plot_barplot(
        tf_acts, 'isHuman', top=25, vertical=True,
        save=config.PATHS.LOGS/'appendix'/'ulm{}.png'.format(selected.pid.values[0]))

    # Extract logFCs and pvals
    logFCs = selected[['logfoldchanges']].T.rename(index={'logfoldchanges': 'isHuman'})
    pvals = selected[['pvals_adj']].T.rename(index={'pvals_adj': 'isHuman'})

    # Plot
    if volcano is not None:
        dc.plot_volcano(
            logFCs, pvals, 'isHuman',
            name=volcano, net=collectri,
            top=10, sign_thr=0.05, lFCs_thr=0.5)

    return tf_acts, tf_pvals


def run_ORA(selected, reactome, threshold=.05):
    print('Launched ORA')
    # Run ora
    enr_pvals = dc.get_ora_df(
        df=selected[selected.pvals_adj < threshold],
        net=reactome,
        source='geneset',
        target='genesymbol',
    )
    enr_pvals = enr_pvals.sort_values('FDR p-value')

    dc.plot_dotplot(
        enr_pvals.head(30), x='Combined score', y='Term', 
        s='Odds ratio', c='FDR p-value', scale=0.5, 
        figsize=(7,10),
        save=config.PATHS.LOGS/'appendix'/'ora{}.png'.format(selected.pid.values[0]))

    for set_ in config.PROTO.ENRICHMENT.GENE_SETS:
        if set_ not in reactome.collection: continue
        dc.plot_running_score(
            selected,
            stat='scores',
            net=reactome,
            source='geneset',
            target='genesymbol',
            set_name='KEGG_NEUROTROPHIN_SIGNALING_PATHWAY'
        )
    return enr_pvals


def up_down_ORA(selected, sets):
    # subset up or down regulated genes
    threshold = config.PROTO.GLM.DEG_THRESHOLD.FDR
    lfc = 0.#config.PROTO.GLM.DEG_THRESHOLD.LFC
    degs_sig = selected[selected.pvals_adj < threshold]
    degs_up = degs_sig[degs_sig.logfoldchanges > lfc]
    degs_dw = degs_sig[degs_sig.logfoldchanges < -lfc]

    enr_up = gp.enrichr(list(degs_up.index), gene_sets=sets)
    gp.dotplot(
        enr_up.res2d, figsize=(3,5), 
        title="Up", cmap = plt.cm.autumn_r, cutoff=.1)
    plt.show()

    enr_dw = gp.enrichr(list(degs_dw.index), gene_sets=sets)
    gp.dotplot(
        enr_dw.res2d, figsize=(3,5), title="Down", 
        cmap = plt.cm.winter_r, size=5)
    plt.show()

    # concat results
    enr_up.res2d['UP_DW'] = "UP"
    enr_dw.res2d['UP_DW'] = "DOWN"
    enr_res = pd.concat([
        enr_up.res2d.sort_values('Adjusted P-value').head(20),
        enr_dw.res2d.sort_values('Adjusted P-value').head(20)])

    # display multi-datasets
    ax = gp.dotplot(enr_res,figsize=(5,7),
                    x='UP_DW',
                    x_order = ["UP","DOWN"],
                    title="GO_BP",
                    # cmap = NbDr.reversed(),
                    size=7,
                    show_ring=True, cutoff=.1)
    ax.set_xlabel("")
    plt.show()
    plt.savefig(config.PATHS.LOGS/'appendix'/'up_down_ora{}.png'.format(selected.pid.values[0]))