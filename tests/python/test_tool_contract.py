
"""
Test for tool contract interface support.
"""

import unittest
import os.path as op

import pbcommand.testkit
import pbcore.data
from pbcore.io import ConsensusReadSet

CCS_DIR = op.dirname(op.dirname(op.dirname(__file__)))


class TestCCSApp(pbcommand.testkit.PbTestApp):
    # FIXME eventually the 'ccs' binary should handle TCI directly
    DRIVER_BASE = op.join(CCS_DIR, "bin", "task_pbccs_ccs")
    REQUIRES_PBCORE = True
    INPUT_FILES = [ pbcore.data.getUnalignedBam() ]
    TASK_OPTIONS = {
        "pbccs.task_options.min_snr": 3,
        "pbccs.task_options.min_read_score": 0.75,
        "pbccs.task_options.min_length": 10,
        "pbccs.task_options.min_passes": 3,
        "pbccs.task_options.min_zscore": -5,
        "pbccs.task_options.max_drop_frac": 0.33,
        "pbccs.task_options.no_polish": True,
    }

    def run_after(self, rtc, output_dir):
        with ConsensusReadSet(rtc.task.output_files[0]) as ds_out:
            self.assertEqual(len(ds_out), 3)


if __name__ == "__main__":
    unittest.main()
