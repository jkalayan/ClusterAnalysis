#!/usr/bin/env python

import numpy as np
import tempfile
from pathlib import Path
from ClusterAnalysis.write.Writers import Writers


coords = np.array([[1,0,0],[0,0,0],[-1,0,0]])
atoms = [8,6,8]

def test_write_xyz():
    # write = Writers()
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir)
        output = output_dir / "test.xyz"
        # print(output)
        Writers.write_xyz([coords], atoms, output, "w")
        assert output.exists()
        with open(output, "r") as f:
            result = f.read()
        # print(repr(result))
        expected_result = ("3\n"
                            "1\n"
                            "   8    1.000000    0.000000    0.000000\n"
                            "   6    0.000000    0.000000    0.000000\n"
                            "   8   -1.000000    0.000000    0.000000\n"
                            )
        assert result == expected_result