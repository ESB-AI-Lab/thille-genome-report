#!/usr/bin/env python
# coding: utf-8

# In[3]:


import pandas as pd


# In[14]:


# PAF file format
# Col	Type	Description
# 1	string	Query sequence name
# 2	int	Query sequence length
# 3	int	Query start (0-based; BED-like; closed)
# 4	int	Query end (0-based; BED-like; open)
# 5	char	Relative strand: "+" or "-"
# 6	string	Target sequence name
# 7	int	Target sequence length
# 8	int	Target start on original strand (0-based)
# 9	int	Target end on original strand (0-based)
# 10	int	Number of residue matches
# 11	int	Alignment block length
# 12	int	Mapping quality (0-255; 255 for missing)


# In[4]:


def read_paf(paf_file):
    """Read a PAF file with a more flexible approach to handle variable number of fields."""
    data = []
    with open(paf_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            # Truncate or extend the parts to match the expected number of columns
            if len(parts) > 11:
                parts = parts[:11]  # Keeping only the first 11 columns as standard PAF
            elif len(parts) < 11:
                # If there are fewer columns, extend with None (or appropriate default)
                parts.extend([None] * (11 - len(parts)))
                
            data.append(parts)
    
    # Convert to DataFrame
    columns = ['query_name', 'query_length', 'query_start', 'query_end',
               'strand', 'ref_name', 'ref_length', 'ref_start', 'ref_end',
               'residue_matches', 'alignment_block_length']
    df = pd.DataFrame(data, columns=columns[:len(data[0])])
    
    # Convert numeric columns
    numeric_cols = ['query_length', 'query_start', 'query_end',
                    'ref_length', 'ref_start', 'ref_end',
                    'residue_matches', 'alignment_block_length']
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    return df


# In[11]:


def map_positions(intersect_file, paf_file):
    # Read intersected genes
    intersect_df = pd.read_csv(intersect_file, sep='\t', header=None,
                               names=['gene', 'gene_start', 'gene_end',
                                      'match_name', 'match_start', 'match_end', 
                                      'strand', 'overlap'],
                               comment='#')  # Skips lines starting with '#'
    
    # Convert to numeric types where necessary
    intersect_df[['gene_start', 'gene_end', 'match_start', 'match_end']] =         intersect_df[['gene_start', 'gene_end', 'match_start', 'match_end']].apply(pd.to_numeric, errors='coerce')
    
    # Read PAF alignments
    paf_df = read_paf(paf_file)
    
    # print("Intersect DataFrame:\n", intersect_df.head())  # Debug output
    # print("\n")
    # print("PAF DataFrame:\n", paf_df.head())              # Debug output
    
    results = []

    # Process each intersected gene
    for idx, row in intersect_df.iterrows():
        # Ensure valid numeric data to avoid comparison failures
        if pd.isna(row['gene_start']) or pd.isna(row['gene_end']):
            #print(f"Skipping row due to NaN in start/end: {row}")
            continue
        
        # Find matching PAF alignments
        paf_row = paf_df[(paf_df['ref_name'] == row['match_name']) &
                         (paf_df['ref_start'] <= row['gene_start']) & 
                         (paf_df['ref_end'] >= row['gene_end'])]
        
        # print(paf_row)

        if not paf_row.empty:
            paf_row = paf_row.iloc[0]
            if row['strand'] == '+':
                query_start = paf_row['query_start'] + (row['gene_start'] - paf_row['ref_start'])
                query_end = paf_row['query_start'] + (row['gene_end'] - paf_row['ref_start'])
            elif row['strand'] == '-':
                query_start = paf_row['query_end'] - (row['gene_end'] - paf_row['ref_start'])
                query_end = paf_row['query_end'] - (row['gene_start'] - paf_row['ref_start'])
            else:
                print(f"Skipping row due to unknown strand: {row}")
                continue
            
            results.append({
                'ref': row['gene'],
                'ref_start': row['gene_start'],
                'ref_end': row['gene_end'],
                'qry': paf_row['query_name'],
                'query_start': int(query_start),
                'query_end': int(query_end),
                'strand': row['strand'],
                'overlap': row['overlap']
            })
        else:
            print("No matching alignment found for:", row['gene'])  # Debugging output

    return pd.DataFrame(results)


# In[13]:


intersect_file = 'gwen_to_UoQ_intersect.txt'  
paf_file = 'gwen_to_UoQ.paf' 
output = 'gwen_to_UoQ_positions.txt'

mapped_positions = map_positions(intersect_file, paf_file)

with open(output, 'w') as file:
    file.write("Mapped Positions (Gwen to UoQ):\n")
    file.write(str(mapped_positions))


# In[ ]:




