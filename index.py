from Bio import SeqIO
import fileinput

# Database indexing (non-overlap)
def database_index(fasta,w=11):
    database = {}
    genes = {}
    for db_record in SeqIO.parse(fasta, "fasta"):
        db_seq = str(db_record.seq)
        db_id = str(db_record.id).split(':')[1]
        genes[db_id] = db_seq
        database[db_id] = {}
        for i in range(0,len(db_seq),w):
            if i+w <= len(db_seq):
                db_word = db_seq[i:i+w]
            else:
                db_word = db_seq[i:len(db_seq)] + 'N'*(w+i-len(db_seq))
            if not db_word in database[db_id].keys():
                database[db_id][db_word] = []
            database[db_id][db_word].append(i)
      
    return database,genes

# Generating seeds from query

def seed_list_of_query_generating(query_seq,w=11):
    one_seed = {}
    for i in range(len(query_seq)):
        if i+w <= len(query_seq):
            word = query_seq[i:i+w]
        else:
            word = query_seq[i:len(query_seq)] + 'N'*(w+i-len(query_seq))
        if not word in one_seed:
            one_seed[word] = []
        one_seed[word].append(i)
    return one_seed

def query_seed_preparing(query_file,w=11):
    seeds = {}
    for line in fileinput.input(query_file):
        one_query = line.rstrip('\n')
        one_seed = seed_list_of_query_generating(one_query,w)
        seeds.update({one_query:one_seed})
    return seeds

#one_query = 'GACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATATGTTGGCCCACTCACAGACTGACCGAGCGAACCTGGGGACCCTGCGCTACTACTACAACCAGAGCGAGGACGGTTCTCACACCATCCAGATAATGTATGGCTGCGACGTGGGGCCGGACGGGCGCTTCCTCCGCGTACCGGCAGG'
#one_seed = seed_list_of_query_generating(one_query)
#db,gene = database_index('A_nuc.fasta')
#print(db)
