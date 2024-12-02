import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sanbomics.tools import id_map

def load_count_matrix():
    """
    Load, merge, and preprocess count files into a single count matrix.

    Returns:
    pd.DataFrame: A preprocessed count matrix with gene IDs as columns and sample names as rows,
                  ready for downstream analysis.
    """
    counts_files = {
        "0M": "htseqResults/0M_output.csv",
        "0P": "htseqResults/0P_output.csv",
        "1A": "htseqResults/1a_output.csv",
        "2P": "htseqResults/2P_output.csv",
        "3L": "htseqResults/3L_output.csv",
        "4C": "htseqResults/4C_output.csv",
    }

    # Initialize an empty DataFrame
    count_matrix = None

    # Process each file
    for sample, file_path in counts_files.items():
        # Read the file with tab-separated values
        df = pd.read_csv(file_path, sep="\t", header=None, names=["geneID", sample])

        # Merge all files into one count matrix
        if count_matrix is None:
            count_matrix = df
        else:
            count_matrix = pd.merge(count_matrix, df, on="geneID", how="outer")

    # Set the geneID column as the index
    count_matrix.set_index("geneID", inplace=True)

    # Replace NaN values with zeros (common for gene count data)
    count_matrix.fillna(0, inplace=True)

    # Filter out rows (genes) with zero counts across all samples
    count_matrix = count_matrix[count_matrix.sum(axis=1) > 0]

    count_matrix = count_matrix[count_matrix.index.str.startswith('E')]

    # Transpose the count matrix to match DESeq2 format (samples as rows, genes as columns)
    count_matrix = count_matrix.T

    return count_matrix

def run_combined_test(
    count_matrix, metadata, species='human'
):
    """
    Run DESeq2 differential expression test for the combined group (RS) against the control group.
    Includes gene mapping using Sanbomics.

    Parameters:
    - count_matrix (pd.DataFrame): Transposed count matrix with samples as rows and genes as columns.
    - control_label (str): The label for the control group in the metadata.
    - combined_label (str): The label for the combined treatment group in the metadata.
    - padj_thresh (float): Threshold for adjusted p-value to consider significance.
    - lfc_thresh (float): Threshold for log2 fold change to consider significance.
    - species (str): The species for the gene mapping (default is 'human').

    Returns:
    - pd.DataFrame: A DataFrame with significant results (with mapped symbols) for the combined test.
    """

    # Initialize DeseqDataSet
    dds = DeseqDataSet(counts=count_matrix, metadata=metadata, design_factors="Condition")
    dds.deseq2()

    # Normalize total counts in the AnnData object
    sc.pp.normalize_total(dds, target_sum=1e4)

    sc.tl.pca(dds)

    sc.pl.pca(dds, color = 'Condition', size = 200)

    # Run differential expression test for the combined group
    stat_res = DeseqStats(dds, contrast=('Condition', 'RS', 'C'))

    # Run Wald tests
    stat_res.summary()

    # Return results
    res = stat_res.results_df

    # Map gene IDs to symbols
    mapper = id_map(species=species)
    res['Symbol'] = res.index.map(mapper.mapper)

    return res

def run_individual_tests(
    count_matrix, treatments, metadata, species='human'
):
    """
    Run DESeq2 differential expression tests for individual treatments against the control group.
    Includes gene mapping using Sanbomics.

    Parameters:
    - count_matrix (pd.DataFrame): Transposed count matrix with samples as rows and genes as columns.
    - treatments (list): List of treatment labels to test against the control.
    - control_label (str): The label for the control group in the metadata.
    - padj_thresh (float): Threshold for adjusted p-value to consider significance.
    - lfc_thresh (float): Threshold for log2 fold change to consider significance.
    - species (str): The species for the gene mapping (default is 'human').

    Returns:
    - dict: A dictionary with treatment labels as keys and significant results (with mapped symbols) as DataFrames.
    """

    # Initialize DeseqDataSet
    dds = DeseqDataSet(counts=count_matrix, metadata=metadata, design_factors="Condition")
    dds.deseq2()

    # Normalize total counts in the AnnData object
    sc.pp.normalize_total(dds, target_sum=1e4)

    sc.tl.pca(dds)

    sc.pl.pca(dds, color = 'Condition', size = 200)

    results = {}  # Store results for each treatment
    mapper = id_map(species=species)  # Initialize gene ID mapper

    for treatment in treatments:
        # Run differential expression test for each treatment vs control
        stat_res = DeseqStats(dds, contrast=('Condition', treatment, 'C'))

        # Run Wald tests
        stat_res.summary()

        # Return results
        res = stat_res.results_df

        # Map gene IDs to symbols
        res['Symbol'] = res.index.map(mapper.mapper)

        # Save the results to the dictionary
        results[treatment] = res

    return results

# Load the count matrix
count_matrix = load_count_matrix()

# Display the first few rows of the processed count matrix
print(count_matrix)

