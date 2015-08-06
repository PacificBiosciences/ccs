
"""
Test for tool contract interface support.
"""

import unittest
import os.path as op

import pbcommand.testkit

CCS_DIR = op.dirname(op.dirname(op.dirname(__file__)))

#XXX totally arbitrary input file - we just need subreads w/P6-C4 chemistry,
# and this one happens to be small and usable enough to test TCI support
DATA = "/mnt/secondary-siv/testdata/kineticsTools/Hpyl_1_2500.bam"


@unittest.skipUnless(op.isfile(DATA), "Missing %s" % DATA)
class CCSTestApp(pbcommand.testkit.PbTestApp):
    # FIXME eventually the 'ccs' binary should handle TCI directly
    DRIVER_BASE = op.join(CCS_DIR, "bin", "task_pbccs_ccs")
    REQUIRES_PBCORE = True
    INPUT_FILES = [ DATA ]
    TASK_OPTIONS = {
        "pbccs.task_options.min_snr": 4,
        "pbccs.task_options.min_read_score": 0.75,
        "pbccs.task_options.min_length": 10,
        "pbccs.task_options.min_passes": 4,
    }


if __name__ == "__main__":
    unittest.main()
