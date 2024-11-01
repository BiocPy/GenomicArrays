import numpy as np
import pyBigWig


def mock_bigwig(filename, max_range, value):
    bw = pyBigWig.open(filename, "w")
    bw.addHeader([("chr1", max_range)], maxZooms=0)

    starts = np.arange(0, max_range, 10, dtype=np.int64)
    ends = starts + 5
    values = np.array([value] * len(starts), dtype=np.float64)
    # np.array(np.random.random_sample(len(starts)), dtype=np.float64)
    chroms = np.array(["chr1"] * len(starts))

    bw.addEntries(chroms, starts, ends=ends, values=values)

    bw.close()


mock_bigwig("test1.bw", 1000, 1)
mock_bigwig("test2.bw", 500, 0.5)
