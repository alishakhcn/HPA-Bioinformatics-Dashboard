import streamlit as st
import pandas as pd
import plotly.express as px

from utils import *

# ==========================================================
# PAGE CONFIG
# ==========================================================
st.set_page_config(
    page_title="Human Protein Atlas Dashboard",
    page_icon="🧬",
    layout="wide"
)

# ==========================================================
# LOAD DATA
# ==========================================================
tissue, cancer = load_data()

# ==========================================================
# HEADER
# ==========================================================
st.title("🧬 Human Protein Atlas Gene Expression Dashboard")

st.markdown("""
### Comparative Transcriptomic Analysis using Human Protein Atlas (HPA)

This dashboard enables interactive comparison of transcriptomic gene expression (nTPM)
between **normal human tissue** and **kidney cancer cell lines** using Human Protein Atlas data.
""")

st.divider()

# ==========================================================
# SIDEBAR
# ==========================================================


st.sidebar.header("🧬 Dashboard Filters")

# Keep only the tissues required for the project
allowed_tissues = [
    "breast",
    "colon",
    "kidney",
    "liver",
    "lung",
    "pancreas",
    "prostate",
    "skin"
]

# Filter the dataframe to these tissues only
tissue = tissue[tissue["Tissue"].isin(allowed_tissues)].copy()

selected_tissue = st.sidebar.selectbox(
    "Select Normal Tissue",
    sorted(tissue["Tissue"].unique())
)
normal = tissue[tissue["Tissue"] == selected_tissue]

genes = sorted(normal["Gene name"].dropna().unique())

selected_gene = st.sidebar.selectbox(
    "Select Gene",
    genes
)

top_n = st.sidebar.slider(
    "Number of Top Biomarkers",
    min_value=5,
    max_value=50,
    value=20,
    step=5
)

# ==========================================================
# BIOMARKERS
# ==========================================================
bio = biomarker_table(
    normal,
    cancer
)

# ==========================================================
# KPI CARDS
# ==========================================================
c1, c2, c3, c4 = st.columns(4)

c1.metric(
    "🧬 Total Genes",
    normal["Gene name"].nunique()
)

c2.metric(
    "🧫 Kidney Cell Lines",
    cancer["Cell line"].nunique()
)

c3.metric(
    "📍 Selected Tissue",
    selected_tissue
)

c4.metric(
    "🏆 Candidate Biomarkers",
    bio["Gene name"].nunique()
)

st.divider()

# ==========================================================
# GENE COMPARISON
# ==========================================================
st.subheader("🧬 Gene Expression Comparison")

ng = normal[
    normal["Gene name"] == selected_gene
]

cg = cancer[
    cancer["Gene name"] == selected_gene
]

if not ng.empty and not cg.empty:

    comparison = pd.DataFrame({

        "Group": [

            f"Normal ({selected_tissue})",

            "Kidney Cancer Cell Lines"

        ],

        "nTPM": [

            float(
                ng["nTPM"].iloc[0]
            ),

            float(
                cg["nTPM"].mean()
            )

        ]

    })

    fig = px.bar(

        comparison,

        x="Group",

        y="nTPM",

        text="nTPM",

        color="Group",

        title=f"{selected_gene}: Normal {selected_tissue.title()} vs Kidney Cancer Cell Lines"

    )

    fig.update_traces(
        texttemplate="%{text:.2f}",
        textposition="outside"
    )

    st.plotly_chart(
        fig,
        use_container_width=True
    )

else:

    st.warning(
        "No expression data available for this gene."
    )

st.divider()

# ==========================================================
# TOP BIOMARKERS
# ==========================================================
st.subheader("🔥 Top Differentially Expressed Biomarkers")

top = top_biomarkers(
    normal,
    cancer,
    top_n
)

fig = px.bar(

    top,

    x="Difference",

    y="Gene name",

    orientation="h",

    color="Fold Change",

    title=f"Top {top_n} Biomarkers"

)

fig.update_layout(
    yaxis=dict(
        autorange="reversed"
    )
)

st.plotly_chart(
    fig,
    use_container_width=True
)

st.divider()

# ==========================================================
# HEATMAP
# ==========================================================
st.subheader("🌡 Heatmap of Top Biomarkers")

heat = heatmap_data(
    cancer,
    normal,
    top_n
)

fig2 = px.imshow(

    heat,

    aspect="auto",

    color_continuous_scale="Viridis",

    labels=dict(
        color="nTPM"
    )

)

st.plotly_chart(
    fig2,
    use_container_width=True
)

st.divider()


# ==========================================================
# HISTOGRAM - Gene Expression Distribution
# ==========================================================

st.subheader("📈 Distribution of Kidney Cancer Gene Expression")

# Remove zero expression values for log scale
cancer_positive = cancer[cancer["nTPM"] > 0].copy()

# Log transform expression values
cancer_positive["log_nTPM"] = np.log2(cancer_positive["nTPM"] + 1)

fig3 = px.histogram(
    cancer_positive,
    x="log_nTPM",
    nbins=50,
    title="Log2 Transformed Gene Expression Distribution",
    labels={
        "log_nTPM": "Log2(nTPM + 1)",
        "count": "Number of Genes"
    },
    color_discrete_sequence=["royalblue"]
)

fig3.update_layout(
    xaxis_title="Log2(nTPM + 1)",
    yaxis_title="Number of Genes",
    bargap=0.1
)

st.plotly_chart(
    fig3,
    use_container_width=True
)

st.divider()
# ==========================================================
# BIOMARKER TABLE
# ==========================================================
st.subheader("🧬 Biomarker Ranking")

st.dataframe(

    bio[
        [
            "Gene name",
            "Normal nTPM",
            "Cancer nTPM",
            "Difference",
            "Fold Change"
        ]
    ].style.background_gradient(

        subset=["Difference"],

        cmap="Blues"

    ),

    use_container_width=True,

    height=500

)

st.download_button(

    "📥 Download Biomarker CSV",

    bio.to_csv(index=False),

    "kidney_biomarkers.csv",

    mime="text/csv"

)

st.divider()

# ==========================================================
# GENE INTERPRETATION
# ==========================================================
stats = gene_statistics(
    normal,
    cancer,
    selected_gene
)

if stats:

    st.subheader("🧠 Biological Interpretation")

    st.info(f"""

**Selected Gene : {selected_gene}**

• Normal Expression : **{stats['Normal']} nTPM**

• Cancer Expression : **{stats['Cancer']} nTPM**

• Fold Change : **{stats['Fold Change']}**

• Status : **{stats['Status']}**

This comparison highlights potential transcriptomic differences between normal tissue and kidney cancer cell lines and may assist in identifying candidate biomarkers for further biological investigation.

""")

st.divider()

# ==========================================================
# DATASET PREVIEW
# ==========================================================
with st.expander("📄 Preview Normal Tissue Dataset"):

    st.dataframe(
        normal.head(10)
    )

with st.expander("📄 Preview Kidney Cancer Cell Line Dataset"):

    st.dataframe(
        cancer.head(10)
    )

st.divider()

# ==========================================================
# ABOUT
# ==========================================================
st.subheader("ℹ️ About This Dashboard")

st.write("""
This dashboard was developed using **Python**, **Streamlit**, **Pandas**, and **Plotly**
to visualize Human Protein Atlas (HPA) transcriptomic data.

It enables users to compare normal tissue expression with kidney cancer cell lines,
explore differentially expressed genes, identify candidate biomarkers,
and interpret transcriptomic expression patterns through interactive visualizations.
""")

st.divider()

st.caption(
"Developed by **Alisha Khan (2409036)** | MSc Big Data Analytics II | Bioinformatics Project | Human Protein Atlas Dashboard"
)