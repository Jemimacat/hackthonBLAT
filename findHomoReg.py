import numpy as np
from scoring import score_nt_seq
from swalgorithm import smith_waterman

def mode(list):
    counts = np.bincount(list)
    mode = np.argmax(counts)
    return mode

def scaning_and_extending(one_seed,database,word_size=11,threshold=11,max_gap=5,max_dist=300):
    segment_hits = {}
    ## search by genes
    for gene in database.keys():
        scores = {}
        pos_diff = {}
        db_words = {}
        q_words = {}
        for db_word in database[gene].keys():
            db_pos = sorted(database[gene][db_word])
            for q_word in one_seed.keys():
                this_score = score_nt_seq(db_word,q_word)
                if this_score >= threshold:
                    q_pos = sorted(one_seed[q_word])
                    for i in db_pos:
                        if not i in scores.keys():
                            scores[i] = {}
                        for j in q_pos:
                            scores[i][j] = i-j ## pos diff
                            db_words[i] = db_word
                            q_words[j] = q_word
                            if not i-j in pos_diff.keys():
                                pos_diff[i-j] = 0
                            pos_diff[i-j] += 1
        (diff,count) = sorted(pos_diff.items(),key=lambda d:d[1], reverse=True)[0]
        ## find homologous regions
        db_pos = []
        q_pos = []
        for i in sorted(scores.keys()):
            (j,this_diff) = sorted(scores[i].items(),key=lambda d:d[1],reverse=False)[0]
            if this_diff <= max_gap:
                if db_pos[-1]:
                    if i-db_pos[-1]<=word_size+max_dist:
                        db_pos.append(i)
                        q_pos.append(j)
                    else:
                        db_pos = [i]
                        q_pos = [j]
                else:
                    db_pos = [i]
                    q_pos = [j]

        ## extending
        db_seq = ''
        q_seq = ''
        for k in range(len(db_pos)):
            if k == 0:
                db_seq += db_words[db_pos[k]]
                q_seq += q_words[q_pos[k]]
            else:
                if q_pos[k] - q_pos[k-1] == word_size:
                    ## continue sequence
                    if db_pos[k] - db_pos[k-1] == word_size:
                        db_seq += db_words[db_pos[k]]
                        q_seq += q_words[q_pos[k]]
                    ## deletion
                    else:
                        p = db_pos[k-1] + word_size
                        while p < db_pos[k]:
                            db_seq += db_words[p]
                            q_seq += '-'*word_size
                            p += word_size
                        db_seq += db_words[db_pos[k]]
                        q_seq += q_words[q_pos[k]]
                
                elif q_pos[k] - q_pos[k-1] > word_size:
                    ## insertion
                    if db_pos[k] - db_pos[k-1] == word_size:
                        for x in range(q_pos[k-1]+1,q_pos[k]):
                            q_seq += q_words[x][-1]
                            db_seq += '-'
                        q_seq += q_words[q_pos[k]]
                        db_seq += db_words[db_pos[k]]
                    else:
                        break
                else:
                    ## insertion
                    if db_pos[k] - db_pos[k-1] == word_size:
                        for x in range(q_pos[k] - q_pos[k-1],word_size):
                            db_seq += '-'
                        q_seq += q_words[q_pos[k]]
                        db_seq += db_words[db_pos[k]]
                    else:
                        break

        layer = ''
        score = score_nt_seq(q_seq,db_seq)
        consensus = 0
        q_length = 0
        for c in range(0,len(q_seq)):
            if q_seq[c] != '-':
                q_length += 1
            if q_seq[c] == db_seq[c]:
                consensus += 1
                layer += ' '
            else:
                layer += '+'
        
        if not score in segment_hits.keys():
            segment_hits[score] = {}
        if not consensus in segment_hits[score].keys():
            segment_hits[score][consensus] = {}
        if not q_length in segment_hits[score][consensus].keys():
            segment_hits[score][consensus][q_length] = {}

        if not gene in segment_hits[score][consensus][q_length].keys():
            segment_hits[score][consensus][q_length][gene] = {}
        if not q_seq in segment_hits[score][consensus][q_length][gene].keys():
            segment_hits[score][consensus][q_length][gene][q_seq] = {}
        if not db_seq in segment_hits[score][consensus][q_length][gene][q_seq].keys():
            segment_hits[score][consensus][q_length][gene][q_seq][db_seq] = {}
    
        segment_hits[score][consensus][q_length][gene][q_seq][db_seq][(layer,q_pos[0],q_pos[-1]+word_size-1,db_pos[0],db_pos[-1]+word_size)] = 1

    return segment_hits



        


