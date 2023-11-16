import matplotlib.pyplot as plt
import seaborn.objects as so
import subprocess

import pandas as pd
import numpy as np
import scanpy as sc
import decoupler as dc
from pathlib import Path

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

    gsea_genesets = geneset_size.index[
        (geneset_size > 5) & (geneset_size < 500)]
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


def run_ULM(selected, collectri):
    print('Launched ULM for TF slope t-stat')
    mat = selected[['logfoldchanges']].T.rename(
        index={'logfoldchanges': 'isHuman'})
    tf_acts, tf_pvals = dc.run_ulm(mat=mat, net=collectri)
    ' '.join(tf_pvals.T[tf_pvals.T.isHuman < .01].index)

    dc.plot_barplot(tf_acts, 'isHuman', top=25, vertical=True)

    # Extract logFCs and pvals
    logFCs = selected[['logfoldchanges']].T.rename(index={'logfoldchanges': 'isHuman'})
    pvals = selected[['pvals_adj']].T.rename(index={'pvals_adj': 'isHuman'})

    # Plot
    dc.plot_volcano(
        logFCs, pvals, 'isHuman', name='CTNNB1', net=collectri,
        top=10, sign_thr=0.05, lFCs_thr=0.5)

    return tf_acts, tf_pvals


def run_ORA(selected, reactome, threshold=.005):
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
        s='Odds ratio', c='FDR p-value', scale=0.5, figsize=(7,10))

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