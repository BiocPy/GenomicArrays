from dataclasses import dataclass
from typing import Dict, Literal

import numpy as np

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


@dataclass
class MatrixOptions:
    """Optional arguments for the ``matrix`` store for :py:func:`~genomicarrays.build_genomicarray.build_genomicarray`.

    Attributes:
        matrix_attr_name:
            Name of the matrix to be stored in the TileDB file.
            Defaults to "data".

        skip:
            Whether to skip generating matrix TileDB.
            Defaults to False.

        dtype:
            NumPy dtype for the values in the matrix.
            Defaults to np.uint16.

            Note: make sure the matrix values fit
            within the range limits of chosen-dtype.

        tiledb_store_name:
            Name of the TileDB file.
            Defaults to `coverage`.

        chunk_size:
            Size of chunks for parallel processing.

        compression:
            TileDB compression filter (None, 'gzip', 'zstd', 'lz4').

        compression_level:
            Compression level (1-9).
    """

    skip: bool = False
    matrix_attr_name: str = "data"
    dtype: np.dtype = np.float32
    tiledb_store_name: str = "coverage"
    chunk_size: int = 1000
    compression: Literal["zstd", "gzip", "lz4"] = "zstd"
    compression_level: int = 4

    def __post_init__(self):
        """Validate configuration."""
        if self.compression not in {"zstd", "gzip", "lz4", None}:
            raise ValueError(f"Unsupported compression: {self.compression}")

        if not 1 <= self.compression_level <= 9:
            raise ValueError(f"Invalid compression level: {self.compression_level}")


@dataclass
class SampleMetadataOptions:
    """Optional arguments for the ``sample`` store for :py:func:`~genomicarrays.build_genomicarray.build_genomicarray`.

    Attributes:
        skip:
            Whether to skip generating sample TileDB.
            Defaults to False.

        dtype:
            NumPy dtype for the sample dimension.
            Defaults to np.uint32.

            Note: make sure the number of samples fit
            within the integer limits of chosen dtype.

        tiledb_store_name:
            Name of the TileDB file.
            Defaults to "sample_metadata".

        column_types:
            A dictionary containing column names as keys
            and the value representing the type to in
            the TileDB.

            If `None`, all columns are cast as 'ascii'.
    """

    skip: bool = False
    dtype: np.dtype = np.uint32
    tiledb_store_name: str = "sample_metadata"
    column_types: Dict[str, np.dtype] = None


@dataclass
class FeatureAnnotationOptions:
    """Optional arguments for the ``feature`` store for :py:func:`~genomicarrays.build_genomicarray.build_genomicarray`.

    Attributes:
        skip:
            Whether to skip generating sample TileDB.
            Defaults to False.

        dtype:
            NumPy dtype for the sample dimension.
            Defaults to np.uint32.

            Note: make sure the number of features fit
            within the integer limits of chosen dtype.

        tiledb_store_name:
            Name of the TileDB file.
            Defaults to "feature_annotations".

        column_types:
            A dictionary containing column names as keys
            and the value representing the type to in
            the TileDB.

            If `None`, all columns are cast as 'ascii'.

        aggregate_function:
            A callable to summarize the values in a given
            interval. Defaults to NumPy's `nanmean`.
    """

    skip: bool = False
    dtype: np.dtype = np.uint32
    tiledb_store_name: str = "feature_annotations"
    column_types: Dict[str, np.dtype] = None
    aggregate_function: callable = np.nanmean
