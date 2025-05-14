from os import walk, path
import pandas as pd
import time
import psutil as ps
from concurrent.futures import ThreadPoolExecutor


def load_data(base_path: str = "data") -> pd.DataFrame:
    """
    Load data from specified path into a single dataframe.

    Parameters
    ----------
    base_path : str (optional)
        The path where the downloaded folders from GDC are located. Can be relative or absolute.

    Returns
    -------
    A single dataframe containing all datapoints.
    """
    dataframes = []
    for root, _, files in walk(base_path):
        # Skip top-level directory as it does not contain any files.
        if not files:
            continue
        for file in files:
            if file == "annotations.txt":
                continue
            if file == "MANIFEST.txt":
                continue
            filepath = path.join(root, file)
            data = pd.read_csv(filepath, encoding="utf-8", sep="\t", )
            dataframes.append(data)
    if not dataframes:
        raise ValueError("Could not find any matching files, wrong directory")
    dataset = pd.concat(dataframes)
    return dataset


def calculate_average_mirna_expression(dataset: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate average expression of each miRNA over the provided dataset.

    Parameters
    ----------
    dataset : pd.DataFrame
        A dataframe that contains the columns 'miRNA_ID' and 'reads_per_million_miRNA_mapped'.

    Returns
    -------
    A dataframe containing one row for each unique miRNA ID in the dataset provided as a parameter.
    Please note that no checks are provided if every miRNA has the exact same amount of samples as
    the others.
    """
    mirna_expressions = []
    counter = 0
    for mirna_id in dataset["miRNA_ID"].unique():  # type: ignore
        mirna_subset = dataset[dataset["miRNA_ID"] == mirna_id]
        mirna_expressions.append(
            {
                "mirna_id": mirna_id,
                "average_expression": mirna_subset[
                    "reads_per_million_miRNA_mapped"
                ].mean()
            }
        )
        counter += 1
        print(counter)
    mirna_expressions = pd.DataFrame(mirna_expressions)
    return mirna_expressions


def calculate_average_and_variance_mirna_expression(dataset: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate average expression and variance of each miRNA over the provided dataset.

    Parameters
    ----------
    dataset : pd.DataFrame
        A dataframe that contains the columns 'miRNA_ID' and 'reads_per_million_miRNA_mapped'.

    Returns
    -------
    A dataframe containing one row for each unique miRNA ID in the dataset provided as a parameter,
    including columns 'miRNA_ID', 'average_expression', and 'variance_expression'.
    Please note that no checks are provided if every miRNA has the exact same amount of samples as
    the others.
    """
    cpu_before = ps.cpu_percent(interval=None)
    ram_before = ps.virtual_memory().percent
    disk_before = ps.disk_usage('/').percent
    mirna_statistics = []
    for mirna_id in dataset["miRNA_ID"].unique():
        mirna_subset = dataset[dataset["miRNA_ID"] == mirna_id]
        average_expression = mirna_subset["reads_per_million_miRNA_mapped"].mean()
        variance_expression = mirna_subset["reads_per_million_miRNA_mapped"].var()
        stddev_expression = mirna_subset["reads_per_million_miRNA_mapped"].std()
        mirna_statistics.append({
            "miRNA_ID": mirna_id,
            "average_expression": average_expression,
            "variance_expression": variance_expression,
            "standard_deviation": stddev_expression
        })

    cpu_after = ps.cpu_percent(interval=None)
    ram_after = ps.virtual_memory().percent
    disk_after = ps.disk_usage('/').percent

    cpu_usage = cpu_after - cpu_before
    ram_usage = ram_after - ram_before
    disk_usage = disk_after - disk_before

    print("CPU Usage:", cpu_usage, "%")
    print("RAM Usage:", ram_usage, "%")
    print("Disk Usage:", disk_usage, "%")

    return pd.DataFrame(mirna_statistics)


start = time.time()
miRNA_data = load_data("./miRNA Files")
print(miRNA_data)
analysis = calculate_average_and_variance_mirna_expression(miRNA_data)
analysis.to_csv("miRNA_output.txt", sep='\t', index=True)
print(analysis)
end = time.time()
print("Analysis Time: %ds" % (end - start))
