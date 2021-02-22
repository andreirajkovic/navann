"""
vcf annotator tests: version 0.0.1
Author: Andrei Rajkovic
Contact: rajkovic.1@osu.edu
"""

from annotator.VariantFileAnnotator import parseDepthData
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.pandas import column, data_frames, range_indexes



@given(
    data_frames([column('A', elements=st.text('01235', min_size=1, max_size=10)), column('B', elements=st.text('01235', min_size=1, max_size=10))], index=range_indexes(min_size=1, max_size=100))
    )
def test_parseDepthData(df):
    """
    test the generation of sample_pd
    """
    df['DPR'] = df['A'] + ',' + df['B']  # emulate the DPR data
    parseDepthData(df)