import pyreadr

# Load the rds file
result = pyreadr.read_r('06_Cb_BTM_table.rds')

# pyreadr returns a dictionary; the data is usually under the 'None' key 
# because rds files contain a single object without a specific name.
df = result[None]

# Now you can use it like a regular pandas DataFrame

for col in df.columns:
    if(len(df[col].unique())<100 and col not in ["gene", "n_mutations",""]):
        print(col + ": " + str(df[col].unique()))
        

print(df.columns)

print(df[["segment_id","sample_id","gene","clock_rank", "best_K"]][0:4])

test = df[df["sample_id"]=="ffb4f42b-58e9-40c3-8963-11804f041375"]

print(test["clock_mean"].unique())