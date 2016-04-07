#!/usr/bin/env python
from __future__ import print_function

import pandas as pd

class BestHits(object):

    def __init__(self, comparison_cols=['E'], query_name_col='q_name', 
                 subject_name_col='s_name', query_length_col='q_len',
                 subject_length_col='q_len'):
        self.comparison_cols = comparison_cols
        self.query_name_col = query_name_col
        self.subject_name_col = subject_name_col
        self.query_length_col = query_length_col
        self.subject_length_col = subject_length_col

    def best_hits(self, aln_df, inplace=True):
        '''Get the best hit for each query in the alignment DataFrame.

        Operates in-place. Sorts the hits by query name and then comparison_cols,
        then uses the drop_duplicates() function to remove all but the
        first hit for each query.

        Args:
            aln_df (DataFrame): The MAF alignment DataFrame.
            inplace (bool): If True, perform the operation in-place and
                return the same DataFrame. If False, return a copy.
        Returns:
            DataFrame with the best hits.
        '''

        if inplace:
            aln_df.sort_values(by=[self.query_name_col] + self.comparison_cols,
                               inplace=True)
            aln_df.drop_dupliates(subset=self.query_name_col, inplace=True)
            return aln_df
        else:
            return aln_df.sort_values(
                       by=[self.query_name_col] + self.comparison_cols
                   ).drop_duplicates(subset=self.query_name_col)

    def reciprocal_best_hits(self, aln_df_A, aln_df_B, inplace=False):
        '''Given to DataFrames with reciprocal MAF alignments, get the
        reciprocal best hits.

        Uses the best_hits function to get best hits for each DataFrame,
        and then does an inner join to find the reciprocals.

        Args:
            aln_df_A (DataFrame): The query hits.
            aln_df_B (DataFrame): The subject hits.
            inplace (bool): Passed to the best_hits calls.
        Returns:
            DataFrame with the reciprocal best hits.
        '''

        aln_df_A = self.best_hits(aln_df_A, inplace=inplace)
        aln_df_B = self.best_hits(aln_df_B, inplace=inplace)

        # Join between subject A and query B
        rbh_df = pd.merge(aln_df_A, aln_df_B, how='inner', 
                          left_on=self.subject_name_col, 
                          right_on=self.query_name_col,
                          suffixes=('_A', '_B'))

        # Renamed columns after join
        left_q_col = self.query_name_col + '_A'
        left_s_col = self.subject_name_col + '_A'
        left_len_col = self.query_length_col + '_A'
        #left_comp_col = self.comparison_col + '_A'
        right_q_col = self.query_name_col + '_B'
        right_s_col = self.subject_name_col + '_B'
        right_len_col = self.query_length_col + '_B'

        # Select those where query A is the same as subject B
        rbh_df = rbh_df[(rbh_df[left_q_col] == rbh_df[right_s_col])]

        # Drop extra columns and rename for clarity
        del rbh_df[left_s_col]
        del rbh_df[right_s_col]
        rbh_df.rename(columns={left_q_col: self.query_name_col, 
                               right_q_col: self.subject_name_col,
                               #left_comp_col: self.comparison_col,
                               left_len_col: self.query_length_col,
                               right_len_col: self.subject_length_col},
                      inplace=True)
        
        return rbh_df


