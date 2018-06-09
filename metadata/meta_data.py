import os,sys
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np

# Gets values that arenot null
def get_notnull(df, key):
    return pd.DataFrame(df[df[key].notnull()][key]).rename_axis(None)

# Takes in key and see if it is in the attribute list 
def get_attributes(keys, file='attributes.txt'):
    attributes = set()
    with open('attributes.txt','r') as f:
        for line in f.readlines():
            if all(x in line for x in keys):
                attributes.add(line.strip())
                return list(attributes)

class MetaData(object):
    
    def __init__(self, input_file_name='chip_atlas_experiment_list.csv'):
        
        raw_df =  pd.read_csv(input_file_name)
        
        self.meta_df = self.extract_meta_data(raw_df)
        self.exp_df = self.meta_df.join(raw_df)
        
        
        
    def extract_meta_data(self, raw_df):
        new_df = raw_df.set_index('Experimental ID')

        exp_list = []
        for idx,row in new_df.iterrows():

            meta_dict = []
            try: 
                meta_list = row['Meta data'].strip().split('||')
                for att in meta_list: 
                    tup = att.strip().split('=')
                    if len(tup) == 1:
                        tup = ['misc', str(tup[0])]
                    elif len(tup) > 2:
                        tup[1:] = [' '.join(tup[1:])]
                    meta_dict.append(tuple(tup))

            # No Metadata Available
            except:
                pass

            meta_dict.append(('exp_ID', idx))        
            exp_list.append(dict(meta_dict))

        return pd.DataFrame(exp_list).set_index('exp_ID')
        
        # Takes in a list of attributes and outputs 
        def filter_attributes(self, attributes):
            if attributes == []:
                raise ValueError("NO ATTRIBUTES PROVIDED")

            for i,att in enumerate(attributes):
                if i == 0:
                    left_df = self.exp_df[self.exp_df[att].notnull()][att]
                else:
                    right_df = self.exp_df[self.exp_df[att].notnull()][att]
                    left_df = pd.concat([left_df, right_df], axis=1)
            return left_df
        
        # Take in key value and see if it exists in attribute list

        

        def combine_attributes(df, combined_col_name, verbose=True):
            ids, attributes = [],[]
            multi_sets = set()

            for index, row in df.iterrows():
                flag=True
                for col in list(df):
                    if not pd.isnull(row[col]) and not flag:
                        multi_sets.add(index)
                    if not pd.isnull(row[col]) and flag:
                        ids.append(index)
                        attributes.append(row[col])
                        flag=False

            value_dict = {}
            value_dict['exp_ID'], value_dict[combined_col_name] = ids, attributes

            if verbose:
                print("Num multi-rows:",len(multi_sets))
                print("Total experiments: ",len(ids))

            return pd.DataFrame.from_dict(value_dict).set_index('exp_ID').rename_axis(None)

        def extract_meta_columns(keys, df, col_name = 'cell_line_X', relevant=None):
            attributes = get_attributes(keys)
            filtered_df = filter_attributes(df, attributes)

            # Filtering step if some attributes are not relevant
            if relevant != None:
                filtered_df = filtered_df[relevant]

            return combine_attributes(filtered_df, combined_col_name = col_name)

        def print_nans(keys, df):
            attributes = get_attributes(keys)
            filtered_df = filter_attributes(df, attributes)   

            for att in attributes:
                val_df = get_notnull(filtered_df, att)
                print(att,': ',val_df.shape[0])

        def get_attribute_frequency(df, num_attributes=25):       
            attribute_dict = {}

            for att in list(df):
                val_df = get_notnull(df, att)
                attribute_dict[att] = val_df.shape[0]

            sorted_d = sorted(attribute_dict.items(), key=lambda x: x[1])
            return sorted_d[-25:]
       