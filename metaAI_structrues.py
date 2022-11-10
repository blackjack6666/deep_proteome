import pandas as pd
import pyarrow.parquet as pq

metadata = 'F:/Meta_predict/stats.parquet'

parquet_f = pq.ParquetFile(metadata)
count = 1
for batch in parquet_f.iter_batches():
    batch_df = batch.to_pandas()
    batch_df.to_csv('F:/Meta_predict/'+str(count)+'.csv', sep=',')
    print (count)
    count+=1


# sample_df = df.iloc[:1000,:]
# sample_df.to_csv('F:/Meta_predict/sample_stats.csv',sep=',')