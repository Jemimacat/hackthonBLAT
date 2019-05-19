from scoring import score_nt
from index import database_index,query_seed_preparing
from findHomoReg import scaning_and_extending

one_query = 'ACAGACTGACCGAGCGAACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGACGGTTCTCACACCATCCAG'
query_file = 'query.txt'
database_file = 'A_nuc.fasta'
db,genes = database_index(database_file)
seeds = query_seed_preparing(query_file)

segment_hits = {}

for query in seeds.keys():
    one_seed = seeds[query]
    one_segment_hits = scaning_and_extending(one_seed, db, word_size=11, threshold=11, max_gap=5, max_dist=300)
    segment_hits.update({query:one_segment_hits})

result = open('result.txt','w+')

for query in segment_hits.keys():
    one_segment_hits = segment_hits[query]
    result.writelines('query: '+query+'\n')
    count = 1
    for score in sorted(one_segment_hits.keys(),reverse=True):
        for consensus in sorted(one_segment_hits[score].keys(),reverse=True):
            for q_length in sorted(one_segment_hits[score][consensus].keys()):
                if q_length < 0:
                    continue
                for gene in sorted(one_segment_hits[score][consensus][q_length].keys()):
                    for q_segment in sorted(one_segment_hits[score][consensus][q_length][gene].keys()):
                        for db_segment in sorted(one_segment_hits[score][consensus][q_length][gene][q_segment].keys()):
                            for items in sorted(one_segment_hits[score][consensus][q_length][gene][q_segment][db_segment].keys()):
                                (layer,q_start,q_end,db_start,db_end) = items
                                result.writelines(str(count) + '\tGene:' + gene + '\tLength:' + str(q_length) + '\tCon:' + str(consensus) + '\tScore:' + str(score) + '\tQuery:' + str(q_start) + ',' + str(q_end) + '\tRef:' + str(db_start) + ',' + str(db_end) +'\n')
                                result.writelines('Ref: '+db_segment +'\n')
                                result.writelines('     '+layer +'\n')
                                result.writelines('Que: '+q_segment +'\n')
                                result.writelines('\n')
                                count += 1
result.close()