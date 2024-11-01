import tempfile

import numpy as np
import pandas as pd
import tiledb

from genomicarrays import build_genomicarray

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def test_ingest_bigwigs():
    tempdir = tempfile.mkdtemp()

    strts = np.arange(300, 600, 20)
    features = pd.DataFrame({"chrom": ["chr1"] * 15, "start": strts, "end": strts + 15})

    build_genomicarray(
        output_path=tempdir,
        files=["tests/data/test1.bw", "tests/data/test2.bw"],
        features=features,
    )

    cfp = tiledb.open(f"{tempdir}/coverage", "r")
    ffp = tiledb.open(f"{tempdir}/feature_annotations", "r")

    features = ffp.df[:]

    assert len(features) == 15

    assert np.allclose(
        cfp.multi_index[0:5, 0]["data"],
        np.repeat(1, 5),
    )
    assert np.allclose(
        cfp.multi_index[0:5, 1]["data"],
        np.repeat(0.5, 5),
    )

    sfp = tiledb.open(f"{tempdir}/sample_metadata", "r")
    samples = sfp.df[:]
    assert len(samples) == 2
