import os
# import sys
# import tempfile
# import json

from test_util import WDLTestCase


class TestPhylotree(WDLTestCase):
    wdl = os.path.join(os.path.dirname(__file__), "..", "run.wdl")

    def test_phylotree(self):
        res = self.run_miniwdl(args=[])
        self.assertEqual(res, {})
