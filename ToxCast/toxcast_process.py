'''
Processes the ToxCast Summary files to prepare it for ingestion in the target safety object. 
'''
from datetime import datetime
import glob
from typing import Tuple

import pandas as pd
import numpy as np

# Datasets paths
assays_path = "INVITRODB_V3_3_SUMMARY/assay_methods_invitrodb_v3_3.xlsx"
targets_path = "INVITRODB_V3_3_SUMMARY/gene_target_information_invitrodb_v3_3.xlsx"
conc_curves_paths = glob.glob("INVITRODB_V3_3_SUMMARY/EXPORT_*.csv")

# Columns of interest
assays_cols = [
    "aeid", "assay_component_endpoint_name",
    "assay_function_type", "intended_target_family",
    "intended_target_family_sub", "assay_component_desc", "assay_design_type",
    "biological_process_target", "organism", "tissue",
    "cell_format", "cell_short_name", "assay_format_type"]
targets_cols = ["official_symbol", "aenm"]
conc_curves_cols = ["aenm", "hitc", "flags"]


def load_data(
    assays_path: str,
    assays_cols: list,
    targets_path: str,
    targets_cols: list,
    conc_curves_paths: str,
    conc_curves_cols: list
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    assays = pd.read_excel(assays_path,
                            sheet_name="assay_merge",
                            engine = 'openpyxl')[assays_cols]

    targets = pd.read_excel(targets_path,
                            sheet_name="Sheet 1",
                            engine = 'openpyxl')[targets_cols]

    conc_curves = [pd.read_csv(f, usecols=conc_curves_cols) for f in conc_curves_paths]
    conc_curves = pd.concat(conc_curves, ignore_index=True)
    return assays, targets, conc_curves

def enrich_assays(
    assays: pd.DataFrame,
    targets: pd.DataFrame,
    conc_curves: pd.DataFrame
) -> pd.DataFrame:
    assays_enriched = (assays
        # Merge assays with targets and concentration fitting curves tables
        .merge(
            targets,
            left_on="assay_component_endpoint_name",
            right_on="aenm",
            how="inner")
        .merge(
            conc_curves,
            left_on="assay_component_endpoint_name",
            right_on="aenm",
            how="inner"
        )
    )
    return assays_enriched

def filter_assays(
    assays: pd.DataFrame
) -> pd.DataFrame:
    assays = (assays
        # Filter human assays
        .query('organism == "human"')
        # Filter out endpoints describing features associated with assay cytotoxic effects
        .query('intended_target_family_sub != "cytotoxicity" or biological_process_target != ["cell death", "cell proliferation"]')
        # Filter out endpoints describing features associated with assay background effects
        .query('assay_function_type != "background control" or assay_design_type != ["background reporter", "viability reporter"] or intended_target_family != "background measurement"')
        # Keep only active assays
        .query('hitc == 1')
        # Filter out potential false positives assays (flagged)
        #.query('flags ')
    )
    assays = assays[pd.isna(assays["flags"]) == True]
    return assays

def adjust_tissue(
    assays: pd.DataFrame
):
    '''
    The tissue value is dropped whenever the assay is cell based.
    The rationale to test in a cell line is often due to the availability and stability
    of the cell line, so the inference between a cell line and the tissue is incorrect
    '''
    assays.loc[~assays["cell_short_name"].isna(), "tissue"] = np.nan


if __name__ == '__main__':
    assays, targets, conc_curves = load_data(
        assays_path, assays_cols,
        targets_path, targets_cols,
        conc_curves_paths, conc_curves_cols
    )
    assays_enriched = enrich_assays(assays, targets, conc_curves)
    adjust_tissue(assays_enriched)
    assays_filtered = filter_assays(assays_enriched)
    out = (assays_filtered
        # Drop duplicates
        .drop_duplicates(
            subset=["assay_component_endpoint_name", "biological_process_target", "official_symbol"]
        )
        # Select features of interest
        [[
            "aeid", "assay_component_endpoint_name", "assay_function_type",
            "intended_target_family", "assay_component_desc",
            "biological_process_target", "tissue", "cell_format",
            "cell_short_name", "assay_format_type", "official_symbol"
        ]]
    )

    print(out["assay_component_endpoint_name"].nunique())
    print(out["official_symbol"].nunique())
    out.to_csv(
        f"ToxCast_{datetime.today().strftime('%Y-%m-%d')}.tsv",
        sep="\t",
        header=True, index=False)    
