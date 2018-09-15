import sys
import string
import csv
from itertools import imap, ifilter
from scipy.stats import norm

def fasta_read(f):
    name, seq = '', ''
    state = 0
    for line in ifilter(lambda x: len(x)>0, imap(string.strip, f)):
        if state == 0:
            if line.startswith('>'):
                name = line[1:].split()[0]
                state = 1
            else:
                raise Exception('invalid data')
        elif state == 1:
            if line.startswith('>'):
                yield name, seq
                name = line[1:].split()[0]
                seq = ''
            else:
                seq += line
    if len(seq) > 0:
        yield name, seq

def seq_warping(seq, alignform):
    count = 0
    pos = 0
    wseq = []
    for i in range(len(alignform)):
        if alignform[i] >= '0' and alignform[i] <= '9':
            count *= 10
            count += int(alignform[i])
        else:
            if alignform[i] == 'M':
                #print "count:%d"%count
                wseq.extend(seq[pos : pos + count])
                #print("".join(wseq))
                pos += count
            elif alignform[i] == 'D':
                #print "count:%d"%count
                wseq.extend(['D'] * count)
                #print("".join(wseq))
            elif alignform[i] == 'S' or alignform[i] == 'I':
                pos += count
            else:
                raise Exception('invalid data')
            count = 0

    return ''.join(wseq)

def get_SAM(f):
    for line in ifilter(lambda x: len(x) > 0 and x[0] != '@', imap(string.strip, f)):
        data = line.split()
        moleId = data[0]
        label = data[1]
        refName = data[2]
        pos = int(data[3])
        alignform = data[5]
        seq = data[9]
        if alignform == '*':
            continue
        yield label, refName, pos, alignform, seq

def seq_condicate(alignform):
    count = 0
    pos = 0
    condicates = []
    for i in range(len(alignform)):
        if alignform[i] >= '0' and alignform[i] <= '9':
            count *= 10
            count += int(alignform[i])
        else:
            if alignform[i] == 'M':
                if count > 7:
                    condicates.append((pos, pos + count))
                pos += count
            elif alignform[i] == 'D':
                pos += count
            elif alignform[i] == 'S' or alignform[i] == 'I':
                pass
            else:
                raise Exception('invalid data')
            count = 0
    return condicates

def score(data, mean, stdv):
    x = (data-mean)/stdv
    delta = norm.ppf(0.9)
    return norm.cdf(x+delta) - norm.cdf(x-delta)

def integrate_eventalign(path):
    reader = csv.reader(open(path), delimiter = '\t')
    position = -1
    count = 0
    summary = []
    product = 1
    for row in reader:
        if not row[3].isdigit():
            continue
        if row[3] == '0':
            if int(row[1]) == position:
                count = count+1
                event_level_mean = float(row[6])
                model_mean = float(row[10])
                model_stdv = float(row[11])
                if row[9] == 'NNNNNN':
                    product = 0
                else:
                    product = product * score(event_level_mean, model_mean, model_stdv)
            else:
                if not position == -1:
                    summary.append([position, count, product])
                if int(row[1]) > position + 1:
                    position = position + 1
                    while position < int(row[1]):
                        summary.append([position, 0, 0])
                        position = position + 1
                position = int(row[1])
                count = 1
                if row[9] == 'NNNNNN':
                    product = 0
                else:
                    event_level_mean = float(row[6])
                    model_mean = float(row[10])
                    model_stdv = float(row[11])
                    product = score(event_level_mean, model_mean, model_stdv)
        else:
            summary.append([position, count, product])
            break
    return summary

def get_events(f):
    pass

def get_eventAlign(f):
    pass

def condicate_score(condicate_range, eventalign):
    count = 0
    score = 1
    i = condicate_range[0]
    while i < condicate_range[1]:
        score = score * eventalign[i][2]
        count = count + eventalign[i][1]
        i = i+1
    return score**(1.0/count)

def signal_split(eventAlign):
    pass

def base_call(refSeq, events):
    pass

def snp_call(refSeq, events):
    pass

if __name__ == '__main__':
    contigs ={}
    with file(sys.argv[1]) as f:
        for name, seq in fasta_read(f):
            contigs[name] = seq
    samData = []
    with file(sys.argv[2]) as f:
        for label, refName, pos, alignform, seq in get_SAM(f):
            samData.append((label, refName, pos, alignform, seq))
    integrated_eventalign = integrate_eventalign(sys.argv[3])
    for sam in samData:
        label, refName, pos, alignform, seq = sam
        if refName == 'KF278742.1':
            #wseq = seq_warping(pos, seq, alignform)
            #condicates = seq_condicate(alignform)
            wseq, condicates = seq_warp_condicates(pos, seq, contigs[refName], alignform)
            print(contigs[refName][pos - 1: pos + len(wseq)])
            print(wseq)
            print(condicates)
            scores = []
            for pair in condicates:
                temp = condicate_score(pair, integrated_eventalign)
                scores.append(temp)
            print(scores)
            print '\n\n\n'
