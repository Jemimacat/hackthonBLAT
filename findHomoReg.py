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
                q_pos = sorted(one_seed[q_word])
                for i in db_pos:
                    db_words[i] = db_word
                    if not i in scores.keys():
                        scores[i] = {}
                    for j in q_pos:
                        q_words[j] = q_word
                        if this_score >= threshold:
                            scores[i][j] = i-j ## pos diff                       
                            if not i-j in pos_diff.keys():
                                pos_diff[i-j] = 0
                            pos_diff[i-j] += 1
        (diff,count) = sorted(pos_diff.items(),key=lambda d:d[1], reverse=True)[0]
        ## find homologous regions
        scores2 = {}
        db_pos = []
        q_pos = []
        ## first filter
        for i in sorted(scores.keys()):
            for (j,this_diff) in sorted(scores[i].items(),key=lambda d:d[1],reverse=False):
                if abs(this_diff-diff) <= max_gap:
                    if not i in scores2.keys():
                        scores2[i] = {}
                    scores2[i][j] = abs(this_diff - diff)

        for i in sorted(scores2.keys()):
            min_abs_dist = sorted(scores2[i].values())[0]
            for (j,abs_dist) in sorted(scores2[i].items(),key=lambda d:d[1],reverse=False):
                if abs_dist == min_abs_dist:
                    if len(db_pos) > 0:
                        if i-db_pos[-1]<=word_size+max_dist:
                            db_pos.append(i)
                            q_pos.append(j)
                        else:
                            db_pos = [i]
                            q_pos = [j]
                    else:
                        db_pos = [i]
                        q_pos = [j]
                    break ## only fetch the first nearest hit

        ## extending
        db_seq = ''
        q_seq = ''
        for k in range(len(db_pos)):
            if k == 0:
                db_seq += db_words[db_pos[k]]
                q_seq += q_words[q_pos[k]]
            else:
                if db_pos[k] - db_pos[k-1] == word_size:
                    ## continue sequence
                    if q_pos[k] - q_pos[k-1] == word_size:
                        db_seq += db_words[db_pos[k]]
                        q_seq += q_words[q_pos[k]]
                    ## insertion
                    elif q_pos[k] - q_pos[k-1] > word_size:
                        for x in range(q_pos[k-1]+1,q_pos[k]-word_size+1):
                            q_seq += q_words[x][-1]
                            db_seq += '-'
                        q_seq += q_words[q_pos[k]]
                        db_seq += db_words[db_pos[k]]
                    else:
                        break
                else:
                    ## deletion
                    if q_pos[k] - q_pos[k-1] == word_size:
                        for x in range(db_pos[k-1]+1,db_pos[k]-word_size+1):
                            q_seq += '-'
                            if x in db_words.keys():
                                db_seq += db_words[x]
                        q_seq += q_words[q_pos[k]]
                        db_seq += db_words[db_pos[k]]
                    else:
                        break
        if not q_seq:
            continue

        ## check forward
        if db_pos[-1]+word_size in db_words.keys():
            next_db_word = db_words[db_pos[-1]+word_size]
            for i in range(1,word_size+1):
                if q_pos[-1]+i in q_words.keys() and q_words[q_pos[-1]+i][-1] == next_db_word[i-1]:
                    db_seq += next_db_word[i-1]
                    q_seq += q_words[q_pos[-1]+i][-1]
                else:
                    break
        ## check backward
        if db_pos[0]-word_size in db_words.keys():
            prior_db_word = db_words[db_pos[0]-word_size]
            for i in range(1,word_size+1):
                if q_pos[0]-i in q_words.keys() and q_words[q_pos[0]-i][0] == prior_db_word[-i]:
                    db_seq = prior_db_word[-i] + db_seq
                    q_seq = q_words[q_pos[0]-i][0] + q_seq
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



        


