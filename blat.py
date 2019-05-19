from scoring import score_nt
from index import database_index,query_seed_preparing
from findHomoReg import scaning_and_extending

one_query = 'ACAGACTGACCGAGCGAACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGACGGTTCTCACACCATCCAG'
query_file = 'query.txt'
database_file = 'A_nuc.fasta'
db,genes = database_index(database_file)
seed = query_seed_preparing(query_file)
one_seed = seed[one_query]
segment_hits = scaning_and_extending(one_seed,db)
print(segment_hits)
