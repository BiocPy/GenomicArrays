from typing import Tuple

import numpy as np
import pandas as pd
import pyBigWig as bw

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def extract_bw_values(
    bw_path: str,
    chrom: str,
) -> Tuple[np.ndarray, int]:
    bwfile = bw.open(bw_path)
    if chrom not in bwfile.chroms():
        return None, None

    chrom_length = bwfile.chroms(chrom)
    data = bwfile.values(chrom, 0, chrom_length)
    return np.array(data), chrom_length


# def extract_aspd():
#     df_data = pd.DataFrame(data, columns=["start", "end", "value"])
#     df_data["chrom"] = row.chrom
#     data_nnz = df_data[df_data["value"] > 0]

#     data_nnz["tmp_group"] = (
#         data_nnz["value"] != data_nnz["value"].shift()
#     ).cumsum()
#     results.append(
#         data_nnz.groupby("tmp_group").agg(
#             {"start": "first", "end": "last", "value": "first"}
#         )
#     )


def extract_bw_intervals_as_vec(
    bw_path: str,
    intervals: pd.DataFrame,
    agg_func: callable = np.nanmean,
    val_dtype: np.dtype = np.float32,
) -> np.ndarray:
    """Extract data from BigWig for a given region and apply the aggregate function.

    Args:
        bw_path:
            Path to the BigWig file.

        intervals:
            List of intervals to extract.

        agg_func:
            Aggregate function to apply.
            Defaults to np.nanmean.

        val_dtype:
            Dtype of the resulting array.

    Returns:
        A vector with length as the number of intervals,
        a value if the file contains the data for the corresponding
        region or ``np.nan`` if the region is not measured.
    """
    bwfile = bw.open(bw_path)

    results = []
    for row in intervals.itertuples():
        if row.chrom in bwfile.chroms():
            try:
                data = bwfile.values(row.chrom, row.start, row.end, numpy=True)
                if data is not None and len(data) != 0:
                    results.append(agg_func(data))
                else:
                    results.append(np.nan)
            except Exception as _:
                results.append(np.nan)
        else:
            results.append(np.nan)

    return np.array(results, dtype=val_dtype)
