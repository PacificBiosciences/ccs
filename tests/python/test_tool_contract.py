#!/usr/bin/env python

"""
Test for tool contract interface support, using a SubreadSet with filters
applied as input.
"""

import logging
import tempfile
import unittest
import os.path as op
import os

import pbcommand.testkit
from pbcore.io import ConsensusReadSet, SubreadSet

import pbtestdata

logging.basicConfig(level=logging.INFO)
CCS_DIR = op.dirname(op.dirname(op.dirname(__file__)))
CHUNK_INDEX = 1

class TestCCSApp(pbcommand.testkit.PbTestApp):
    # FIXME eventually the 'ccs' binary should handle TCI directly
    DRIVER_BASE = os.environ['__PBTEST_CCS_EXE']
    REQUIRES_PBCORE = True
    INPUT_FILES = [tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name]
    TASK_OPTIONS = {
        "pbccs.task_options.min_snr": 3.0,
        "pbccs.task_options.min_read_score": 0.75,
        "pbccs.task_options.min_length": 10,
        "pbccs.task_options.min_passes": 3,
        "pbccs.task_options.min_zscore": -5.0,
        "pbccs.task_options.max_drop_frac": 0.33,
        "pbccs.task_options.polish": False
    }

    def setUp(self):
        BAM_IN = pbtestdata.get_file("subreads-bam")
        ds = SubreadSet(BAM_IN, strict=True)
        chunks = ds.split(zmws=True, chunks=2, targetSize=2)
        assert len(chunks) == 2
        self.zmw_range = chunks[CHUNK_INDEX].zmwRanges[0][1:3]
        logging.info("zmwRanges[CHUNK_INDEX] = {r}".format(
            r=str(chunks[CHUNK_INDEX].zmwRanges)))
        logging.info("SubreadSet = {f}".format(f=self.INPUT_FILES[0]))
        chunks[CHUNK_INDEX].write(self.INPUT_FILES[0])
        #tf = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name
        #print tf
        #chunks[int(not CHUNK_INDEX)].write(tf)

    def run_after(self, rtc, output_dir):
        with ConsensusReadSet(rtc.task.output_files[0]) as ds_out:
            zmws = set(ds_out.resourceReaders()[0].holeNumber)
            logging.info("ZMWs = {z}".format(z=zmws))
            for z in zmws:
                self.assertTrue(self.zmw_range[0] < z < self.zmw_range[1])


if __name__ == "__main__":
    unittest.main()
