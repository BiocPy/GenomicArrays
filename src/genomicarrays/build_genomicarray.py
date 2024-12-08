"""Build the `GenomicArrayDatset`.

This modules provides tools for converting genomic data from BigWig format to TileDB.
It supports parallel processing for handling large collection of genomic datasets.

Example:

    .. code-block:: python

        import pyBigWig as bw
        import numpy as np
        import tempfile
        from genomicarrays import build_genomicarray, MatrixOptions

        # Create a temporary directory
        tempdir = tempfile.mkdtemp()

        # Read BigWig objects
        bw1 = bw.open("path/to/object1.bw", "r")
        # or just provide the path
        bw2 = "path/to/object2.bw"

        features = pd.DataFrame({
            "seqnames": ["chr1", "chr1"],
            "starts": [1000, 2000],
            "ends": [1500, 2500]
        })

        # Build GenomicArray
        dataset = build_genomicarray(
            features=features
            output_path=tempdir,
            files=[bw1, bw2],
            matrix_options=MatrixOptions(dtype=np.float32),
        )
"""

import os
import warnings
from multiprocessing import Pool
from typing import Union

import pandas as pd
from cellarr import buildutils_tiledb_frame as utf
from pyfaidx import Fasta

from . import build_options as bopt
from . import buildutils_tiledb_array as uta
from . import utils_bw as ubw
from .GenomicArrayDataset import GenomicArrayDataset

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


# TODO: Accept files as a dictionary with names to each dataset.
def build_genomicarray(
    files: list,
    output_path: str,
    features: Union[str, pd.DataFrame],
    genome_fasta: str,
    # genome: Union[str, pd.DataFrame] = "hg38",
    sample_metadata: Union[pd.DataFrame, str] = None,
    sample_metadata_options: bopt.SampleMetadataOptions = bopt.SampleMetadataOptions(),
    matrix_options: bopt.MatrixOptions = bopt.MatrixOptions(),
    feature_annotation_options: bopt.FeatureAnnotationOptions = bopt.FeatureAnnotationOptions(),
    optimize_tiledb: bool = True,
    num_threads: int = 1,
):
    """Generate the `GenomicArrayDatset`.

    All files are expected to be consistent and any modifications
    to make them consistent is outside the scope of this function
    and package.

    Args:
        files:
            List of file paths to `BigWig` files.

        output_path:
            Path to where the output TileDB files should be stored.

        features:
            A :py:class:`~pandas.DataFrame` containing the input genomic
            intervals..

            Alternatively, may provide path to the file containing a
            list of intervals. In this case,
            the first row is expected to contain the column names,
            "seqnames", "starts" and "ends".

            Additionally, the file may contain a column `sequences`,
            to specify the sequence string for each region. Otherwise,
            provide a link to the fasta file using the ``genome_fasta``
            parameter.

        genome_fasta:
            Path to a fasta file containing the sequence information.

            Sequence information will be updated for each region
            in ``features`` if the DataFrame/path does not contain
            a column `sequences`.

        sample_metadata:
            A :py:class:`~pandas.DataFrame` containing the sample
            metadata for each file in ``files``. Hences the number of rows
            in the dataframe must match the number of ``files``.

            Alternatively, may provide path to the file containing a
            concatenated sample metadata across all BigWig files. In this case,
            the first row is expected to contain the column names.

            Additionally, the order of rows is expected to be in the same
            order as the input list of ``files``.

            Defaults to `None`, in which case, we create a simple sample
            metadata dataframe containing the list of datasets, aka
            each BigWig files. Each dataset is named as ``sample_{i}``
            where `i` refers to the index position of the object in ``files``.

        sample_metadata_options:
            Optional parameters when generating ``sample_metadata`` store.

        matrix_options:
            Optional parameters when generating ``matrix`` store.

        feature_annotation_options:
            Optional parameters when generating ``feature_annotation`` store.

        optimize_tiledb:
            Whether to run TileDB's vaccum and consolidation (may take long).

        num_threads:
            Number of threads.
            Defaults to 1.
    """
    if not os.path.isdir(output_path):
        raise ValueError("'output_path' must be a directory.")

    ####
    ## Process genome information
    ####
    # if isinstance(genome, str):
    #     chrom_sizes = pd.read_csv(
    #         "https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes",
    #         sep="\t",
    #         header=None,
    #         names=["chrom", "length"],
    #     )
    # elif isinstance(genome, pd.DataFrame):
    #     chrom_sizes = genome
    #     if "chrom" not in chrom_sizes:
    #         raise ValueError("genome does not contain column: 'chrom'.")

    #     if "length" not in chrom_sizes:
    #         raise ValueError("genome does not contain column: 'length'.")

    # else:
    #     raise TypeError(
    #         "'genome' is not an expected type (either 'str' or 'Dataframe')."
    #     )

    ####
    ## Writing the features aka interval regions
    ####
    if isinstance(features, str):
        input_intervals = pd.read_csv(features, header=0)
    elif isinstance(features, pd.DataFrame):
        input_intervals = features.copy()

        required_cols = {"seqnames", "starts", "ends"}
        if not required_cols.issubset(input_intervals.columns):
            missing = required_cols - set(input_intervals.columns)
            raise ValueError(f"Missing required columns: {missing}")
    else:
        raise TypeError("'input_intervals' is not an expected type (either 'str' or 'Dataframe').")

    # append start index for each interval
    input_intervals["widths"] = input_intervals["ends"] - input_intervals["starts"]
    total_length = None
    if feature_annotation_options.aggregate_function is None:
        total_length = int(input_intervals["widths"].sum())
        counter = input_intervals["widths"].shift(1)
        counter[0] = 0
        input_intervals["genarr_feature_start_index"] = counter.cumsum().astype(int)
    else:
        counter = [1] * len(input_intervals)
        total_length = len(input_intervals)
        counter[0] = 0
        input_intervals["genarr_feature_start_index"] = pd.Series(counter).cumsum().astype(int)

    ends = input_intervals["genarr_feature_start_index"].shift(-1)
    ends.iloc[-1] = total_length
    input_intervals["genarr_feature_end_index"] = ends.astype(int)

    if not feature_annotation_options.skip:
        if "sequence" not in input_intervals.columns:
            gen_fa = Fasta(genome_fasta, as_raw=True)
            sequences = []
            for row in input_intervals.itertuples():
                seq = gen_fa.get_seq(row.seqnames, row.starts, row.ends)
                sequences.append(seq)

            input_intervals["sequences"] = sequences

        _col_types = utf.infer_column_types(input_intervals, feature_annotation_options.column_types)

        if "genarr_feature_index" not in input_intervals.columns:
            input_intervals["genarr_feature_index"] = range(0, len(input_intervals))

        _feature_output_uri = f"{output_path}/{feature_annotation_options.tiledb_store_name}"
        utf.create_tiledb_frame_from_dataframe(_feature_output_uri, input_intervals, column_types=_col_types)

        if optimize_tiledb:
            uta.optimize_tiledb_array(_feature_output_uri)

    ####
    ## Writing the sample metadata file
    ####
    _samples = []
    for idx, _ in enumerate(files):
        _samples.append(f"sample_{idx + 1}")

    if sample_metadata is None:
        warnings.warn(
            "Sample metadata is not provided, each dataset in 'files' is considered a sample",
            UserWarning,
        )

        sample_metadata = pd.DataFrame({"genarr_sample": _samples})
    elif isinstance(sample_metadata, str):
        sample_metadata = pd.read_csv(sample_metadata, header=0)
        if "genarr_sample" not in sample_metadata.columns:
            sample_metadata["genarr_sample"] = _samples
    elif isinstance(sample_metadata, pd.DataFrame):
        if "genarr_sample" not in sample_metadata.columns:
            sample_metadata["genarr_sample"] = _samples
    else:
        raise TypeError("'sample_metadata' is not an expected type.")

    if not sample_metadata_options.skip:
        _col_types = utf.infer_column_types(sample_metadata, sample_metadata_options.column_types)

        _sample_output_uri = f"{output_path}/{sample_metadata_options.tiledb_store_name}"
        utf.create_tiledb_frame_from_dataframe(_sample_output_uri, sample_metadata, column_types=_col_types)

        if optimize_tiledb:
            uta.optimize_tiledb_array(_sample_output_uri)

    ####
    ## Writing the genomic ranges file
    ####
    if not matrix_options.skip:
        _cov_uri = f"{output_path}/{matrix_options.tiledb_store_name}"
        uta.create_tiledb_array(
            _cov_uri,
            matrix_attr_name=matrix_options.matrix_attr_name,
            x_dim_dtype=feature_annotation_options.dtype,
            y_dim_dtype=sample_metadata_options.dtype,
            matrix_dim_dtype=matrix_options.dtype,
            x_dim_length=total_length,
            y_dim_length=len(files),
            is_sparse=False,
        )

        all_bws_options = [
            (
                _cov_uri,
                input_intervals,
                bwpath,
                idx,
                feature_annotation_options.aggregate_function,
                total_length,
            )
            for idx, bwpath in enumerate(files)
        ]

        if num_threads > 1:
            with Pool(num_threads) as p:
                p.map(_wrapper_extract_bwinfo, all_bws_options)
        else:
            for opt in all_bws_options:
                _wrapper_extract_bwinfo(opt)

        if optimize_tiledb:
            uta.optimize_tiledb_array(_cov_uri)

    return GenomicArrayDataset(
        dataset_path=output_path,
        sample_metadata_uri=sample_metadata_options.tiledb_store_name,
        feature_annotation_uri=feature_annotation_options.tiledb_store_name,
        matrix_tdb_uri=matrix_options.tiledb_store_name,
    )


def _write_intervals_to_tiledb(outpath, intervals, bwpath, bwidx, agg_func, total_length):
    """Wrapper to extract the data for the given intervals from the bigwig file and write the output to the tiledb
    file."""
    data = ubw.wrapper_extract_bw_values(
        bw_path=bwpath,
        intervals=intervals,
        agg_func=agg_func,
        total_length=total_length,
    )

    if data is not None and len(data) > 0:
        uta.write_frame_intervals_to_tiledb(outpath, data=data, y_idx=bwidx)


def _wrapper_extract_bwinfo(args):
    """Wrapper for multiprocessing multiple files and intervals."""
    counts_uri, input_intervals, bwpath, idx, agg_func, total_length = args
    return _write_intervals_to_tiledb(counts_uri, input_intervals, bwpath, idx, agg_func, total_length)
