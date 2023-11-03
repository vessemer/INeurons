import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path

from ..configs import config


NUM_OF_CELL_PER_DONOR = 75


def aggregate_and_filter(
    adata,
    cell_identity,
    donor_key="sample",
    cell_identity_key="clusters",
    obs_to_keep=[],# which additional metadata to keep, e.g. gender, age, etc.
    replicates_per_patient=3
):

    # subset adata to the given cell identity
    adata_cell_pop = adata[adata.obs[cell_identity_key] == cell_identity].copy()
    # check which donors to keep according to the number of cells specified with NUM_OF_CELL_PER_DONOR
    size_by_donor = adata_cell_pop.obs.groupby([donor_key]).size()
    donors_to_drop = [
        donor
        for donor in size_by_donor.index
        if size_by_donor[donor] <= NUM_OF_CELL_PER_DONOR
    ]
    if len(donors_to_drop) > 0:
        print("Dropping the following samples:")
        print(donors_to_drop)
    df = list()

    adata_cell_pop.obs[donor_key] = adata_cell_pop.obs[donor_key].astype("category")
    for i, donor in enumerate(donors := adata_cell_pop.obs[donor_key].cat.categories):
        print(f"\tProcessing donor {i+1} out of {len(donors)}...", end="\r")
        if donor not in donors_to_drop:
            adata_donor = adata_cell_pop[adata_cell_pop.obs[donor_key] == donor]
            # create replicates for each donor
            indices = list(adata_donor.obs_names)
            random.shuffle(indices)
            indices = np.array_split(np.array(indices), replicates_per_patient)
            for i, rep_idx in enumerate(indices):
                adata_replicate = adata_donor[rep_idx]
                # specify how to aggregate: sum gene expression for each gene for each donor and also keep the condition information
                agg_dict = {gene: "sum" for gene in adata_replicate.var_names}
                for obs in obs_to_keep:
                    agg_dict[obs] = "first"
                # create a df with all genes, donor and condition info
                df_donor = pd.DataFrame(adata_replicate.X.A)
                df_donor.index = adata_replicate.obs_names
                df_donor.columns = adata_replicate.var_names
                #df_donor = df_donor.join(adata_replicate.obs[obs_to_keep])
                # aggregate
                df_donor = df_donor.sum()#.agg(agg_dict)
                df_donor = pd.concat([df_donor, adata_replicate.obs[obs_to_keep].iloc[0]])
                df_donor[donor_key] = donor
                df_donor.name = f"donor_{donor}_{i}"
                df.append(df_donor)

    df = pd.DataFrame(df, columns=[*adata_cell_pop.var_names, *obs_to_keep])
    print("\n")
    # create AnnData object from the df
    adata_cell_pop = sc.AnnData(
        df[adata_cell_pop.var_names],
        obs=df.drop(columns=adata_cell_pop.var_names),
    )
    adata_cell_pop.var['holder'] = np.arange(
        len(adata_cell_pop.var.index))

    return adata_cell_pop


def aggregate_with_replace(
    adata,
    cell_identity,
    donor_key="sample",
    cell_identity_key="cluster",
    obs_to_keep=[],# which additional metadata to keep, e.g. gender, age, etc.
):

    # subset adata to the given cell identity
    adata_cell_pop = adata[adata.obs[cell_identity_key] == cell_identity].copy()
    # check which donors to keep according to the number of cells specified with NUM_OF_CELL_PER_DONOR
    size_by_donor = adata_cell_pop.obs.groupby([donor_key]).size()
    donors_to_drop = [
        donor
        for donor in size_by_donor.index
        if size_by_donor[donor] <= NUM_OF_CELL_PER_DONOR
    ]
    if len(donors_to_drop) > 0:
        print("Dropping the following samples:")
        print(donors_to_drop)
    df = list()

    rs = np.random.RandomState(config.PROTO.SEED)
    adata_cell_pop.obs[donor_key] = adata_cell_pop.obs[donor_key].astype("category")
    for i, donor in enumerate(donors := adata_cell_pop.obs[donor_key].cat.categories):
        print(f"\tProcessing donor {i+1} out of {len(donors)}...", end="\r")
        if donor not in donors_to_drop:
            adata_donor = adata_cell_pop[adata_cell_pop.obs[donor_key] == donor]
            # create replicates for each donor
            indices = list(adata_donor.obs_names)            

            for i in range(config.PROTO.PSEUDO_BULK.REPS):
                rep_idx = rs.choice(np.array(indices), config.PROTO.PSEUDO_BULK.SIZE)
                adata_replicate = adata_donor[rep_idx]
                # specify how to aggregate: sum gene expression for each gene for each donor and also keep the condition information
                agg_dict = {gene: "sum" for gene in adata_replicate.var_names}
                for obs in obs_to_keep:
                    agg_dict[obs] = "first"
                # create a df with all genes, donor and condition info
                df_donor = pd.DataFrame(adata_replicate.X.A)
                df_donor.index = adata_replicate.obs_names
                df_donor.columns = adata_replicate.var_names
                #df_donor = df_donor.join(adata_replicate.obs[obs_to_keep])
                # aggregate
                df_donor = df_donor.sum()#.agg(agg_dict)
                df_donor = pd.concat([df_donor, adata_replicate.obs[obs_to_keep].iloc[0]])
                df_donor[donor_key] = donor
                df_donor.name = f"donor_{donor}_{i}"
                df.append(df_donor)

    df = pd.DataFrame(df, columns=[*adata_cell_pop.var_names, *obs_to_keep])
    print("\n")
    # create AnnData object from the df
    adata_cell_pop = sc.AnnData(
        df[adata_cell_pop.var_names],
        obs=df.drop(columns=adata_cell_pop.var_names),
    )
    adata_cell_pop.var['holder'] = np.arange(
        len(adata_cell_pop.var.index))

    return adata_cell_pop


def aggregate_wo_replace(
    adata,
    cell_identity,
    donor_key="sample",
    cell_identity_key="clusters",
    obs_to_keep=[],# which additional metadata to keep, e.g. gender, age, etc.
):

    # subset adata to the given cell identity
    adata_cell_pop = adata[adata.obs[cell_identity_key] == cell_identity].copy()
    # check which donors to keep according to the number of cells specified with NUM_OF_CELL_PER_DONOR
    size_by_donor = adata_cell_pop.obs.groupby([donor_key]).size()
    donors_to_drop = [
        donor
        for donor in size_by_donor.index
        if size_by_donor[donor] <= NUM_OF_CELL_PER_DONOR
    ]
    if len(donors_to_drop) > 0:
        print("Dropping the following samples:")
        print(donors_to_drop)

    df = list()

    rs = np.random.RandomState(config.SEED)
    adata_cell_pop.obs[donor_key] = adata_cell_pop.obs[donor_key].astype("category")
    for i, donor in enumerate(donors := adata_cell_pop.obs[donor_key].cat.categories):
        print(f"\tProcessing donor {i+1} out of {len(donors)}...", end="\r")
        if donor in donors_to_drop: continue

        adata_donor = adata_cell_pop[adata_cell_pop.obs[donor_key] == donor]
        # create replicates for each donor
        indices = list(adata_donor.obs_names)            

        #rep_idx = rs.choice(np.array(indices), config.SAMPLING.SIZE, replace=False)
        adata_replicate = adata_donor#[rep_idx]
        # specify how to aggregate: sum gene expression for each gene for each donor and also keep the condition information
        agg_dict = {gene: "sum" for gene in adata_replicate.var_names}
        for obs in obs_to_keep:
            agg_dict[obs] = "first"
        # create a df with all genes, donor and condition info
        df_donor = pd.DataFrame(adata_replicate.X.A)
        df_donor.index = adata_replicate.obs_names
        df_donor.columns = adata_replicate.var_names
        #df_donor = df_donor.join(adata_replicate.obs[obs_to_keep])
        # aggregate
        df_donor = df_donor.sum()#.agg(agg_dict)
        df_donor = pd.concat([df_donor, adata_replicate.obs[obs_to_keep].iloc[0]])
        df_donor[donor_key] = donor
        df_donor.name = f"donor_{donor}_{i}"
        df.append(df_donor)

    df = pd.DataFrame(df, columns=[*adata_cell_pop.var_names, *obs_to_keep])
    print("\n")
    # create AnnData object from the df
    adata_cell_pop = sc.AnnData(
        df[adata_cell_pop.var_names],
        obs=df.drop(columns=adata_cell_pop.var_names),
    )
    adata_cell_pop.var['holder'] = np.arange(
        len(adata_cell_pop.var.index))

    return adata_cell_pop


def gmt_to_decoupler(pth: Path) -> pd.DataFrame:
    """
    Parse a gmt file to a decoupler pathway dataframe.
    """
    from itertools import chain, repeat

    pathways = {}

    with Path(pth).open("r") as f:
        for line in f:
            name, _, *genes = line.strip().split("\t")
            pathways[name] = genes

    return pd.DataFrame.from_records(
        chain.from_iterable(zip(repeat(k), v) for k, v in pathways.items()),
        columns=["geneset", "genesymbol"],
    )


def prep_anndata(adata_):
    def fix_dtypes(adata_):
        df = pd.DataFrame(adata_.X.A, index=adata_.obs_names, columns=adata_.var_names)
        df = df.join(adata_.obs)
        return sc.AnnData(df[adata_.var_names], obs=df.drop(columns=adata_.var_names))

    adata_ = fix_dtypes(adata_)
    sc.pp.filter_genes(adata_, min_cells=3)
    return adata_


def merge_datasets(datasets: list):
    var = datasets[0].var.index
    symdiff = set()
    for data in datasets[1:]:
        symdiff = symdiff.union(var.symmetric_difference(data.var.index))
        var = var.intersection(data.var.index)

    for i, data in enumerate(datasets):
        data = data[:, data.var_names.isin(var)]
        datasets[i] = data

    return sc.concat(datasets)


def erase_log():
    with open(config.PATHS.LOGS/"logs.md", 'w') as file:
        file.write("")


def log(*message):
    print(*message)
    with open(config.PATHS.LOGS/"logs.md", 'a') as file:
        if image is not None:
            file.write(
                "![](./appendix/{}.png)\n".format(image))
        file.write("{}  \n".format(message))

def log_image(*message):
    print(*message)
    with open(config.PATHS.LOGS/"logs.md", 'a') as file:
        if image is not None:
            file.write(
                "![](./appendix/{}.png)\n".format(image))
        file.write("{}  \n".format(message))


def log(*message):
    print(*message)
    with open(config.PATHS.LOGS/"logs.md", 'a') as file:
        file.write("{}  \n".format(' '.join(list(map(str,message)))))
