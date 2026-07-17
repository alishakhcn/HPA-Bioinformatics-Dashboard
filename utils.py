import pandas as pd


# ==========================================================
# Load Data
# ==========================================================

def load_data():

    tissue = pd.read_csv(
        "data/rna_tissue_small.tsv",
        sep="\t",
        usecols=[
            "Gene",
            "Gene name",
            "Tissue",
            "nTPM"
        ]
    )

    cancer = pd.read_csv(
        "data/kidney_cellline.tsv",
        sep="\t",
        usecols=[
            "Gene",
            "Gene name",
            "Cell line",
            "nTPM"
        ]
    )

    return tissue, cancer

# ==========================================================
# Tissue List
# ==========================================================
def get_tissues(tissue):

    return sorted(
        tissue["Tissue"].dropna().unique()
    )

# ==========================================================
# Filter Tissue
# ==========================================================
def filter_tissue(tissue, tissue_name):

    return tissue[
        tissue["Tissue"] == tissue_name
    ]


# ==========================================================
# Biomarker Ranking
# ==========================================================
def biomarker_table(normal, cancer):

    cancer_avg = (

        cancer.groupby(
            ["Gene", "Gene name"]
        )["nTPM"]

        .mean()

        .reset_index()

    )

    cancer_avg.rename(

        columns={
            "nTPM": "Cancer nTPM"
        },

        inplace=True

    )

    normal_df = normal.rename(

        columns={
            "nTPM": "Normal nTPM"
        }

    )

    merged = pd.merge(

        normal_df,

        cancer_avg,

        on=[
            "Gene",
            "Gene name"
        ]

    )

    merged["Difference"] = (

        merged["Cancer nTPM"]

        -

        merged["Normal nTPM"]

    )

    merged["Fold Change"] = (

        (merged["Cancer nTPM"] + 1)

        /

        (merged["Normal nTPM"] + 1)

    )

    merged = merged.sort_values(

        by="Difference",

        ascending=False

    )

    return merged


# ==========================================================
# Top Biomarkers
# ==========================================================
def top_biomarkers(
        normal,
        cancer,
        top_n=10
):

    bio = biomarker_table(
        normal,
        cancer
    )

    return bio.head(top_n)


# ==========================================================
# Gene Statistics
# ==========================================================
def gene_statistics(
        normal,
        cancer,
        gene
):

    normal_gene = normal[
        normal["Gene name"] == gene
    ]

    cancer_gene = cancer[
        cancer["Gene name"] == gene
    ]

    if normal_gene.empty or cancer_gene.empty:

        return None

    normal_value = float(
        normal_gene["nTPM"].iloc[0]
    )

    cancer_value = float(
        cancer_gene["nTPM"].mean()
    )

    difference = (
        cancer_value -
        normal_value
    )

    fold_change = (

        (cancer_value + 1)

        /

        (normal_value + 1)

    )

    if fold_change >= 2:

        status = "Highly Overexpressed"

    elif fold_change >= 1:

        status = "Moderately Overexpressed"

    else:

        status = "Downregulated"

    return {

        "Normal":
        round(normal_value, 2),

        "Cancer":
        round(cancer_value, 2),

        "Difference":
        round(difference, 2),

        "Fold Change":
        round(fold_change, 2),

        "Status":
        status

    }


# ==========================================================
# Heatmap Data
# ==========================================================
def heatmap_data(cancer, normal, top_n=20):

    bio = biomarker_table(normal, cancer)

    top_genes = bio.head(top_n)["Gene name"]

    heat = cancer[
        cancer["Gene name"].isin(top_genes)
    ]

    heat = heat.pivot_table(
        index="Gene name",
        columns="Cell line",
        values="nTPM",
        aggfunc="mean"
    )

    return heat


# ==========================================================
# Dataset Summary
# ==========================================================
def dataset_summary(
        normal,
        cancer
):

    summary = {

        "Genes":
        normal["Gene name"].nunique(),

        "Cell Lines":
        cancer["Cell line"].nunique(),

        "Average nTPM":
        round(
            cancer["nTPM"].mean(),
            2
        ),

        "Maximum nTPM":
        round(
            cancer["nTPM"].max(),
            2
        ),

        "Minimum nTPM":
        round(
            cancer["nTPM"].min(),
            2
        )

    }

    return summary