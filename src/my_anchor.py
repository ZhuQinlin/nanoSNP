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
        yield moleId, label, refName, pos, alignform, seq

#pos seq and alignform both from output.sam
def seq_warp_condicates(pos_ref, read, ref, alignform):
    count = 0
    pos_ref = pos_ref - 1
    pos_read = 0
    wseq = []
    condicates = []
    for i in range(len(alignform)):
        if alignform[i] >= '0' and alignform[i] <= '9':
            count *= 10
            count += int(alignform[i])
        else:
            if alignform[i] == 'M':
                #print "count:%d"%count
                flag = 0
                c = 0
                for i in range(count):
                    if flag == 0:
                        if read[pos_read + i] == ref[pos_ref + i]:
                            pos = pos_ref + i
                            flag = 1
                            c = 1
                        else:
                            flag = 0
                    else:   #flag == 1
                        if read[pos_read + i] == ref[pos_ref + i]:
                            c = c+1
                            flag = 1
                        else:
                            flag = 0
                            if c > 7:
                                condicates.append((pos, pos+c))
                if flag == 1:
                    if c > 7:
                        condicates.append((pos, pos+c))
                #condicates.append((pos_ref, pos_ref + count))
                wseq.extend(seq[pos_read : pos_read + count])
                #print("".join(wseq))
                pos_ref += count
                pos_read += count
            elif alignform[i] == 'D':
                #print "count:%d"%count
                wseq.extend(['D'] * count)
                #print("".join(wseq))
                pos_ref += count
            elif alignform[i] == 'S' or alignform[i] == 'I':
                pos_read += count
            else:
                raise Exception('invalid data')
            count = 0

    return ''.join(wseq), condicates

def score(data, mean, stdv):
    x = (data-mean)/stdv
    delta = norm.ppf(0.9)
    return norm.cdf(x+delta) - norm.cdf(x-delta)

def integrate_eventalign(path):
    reader = csv.reader(open(path), delimiter = '\t')
    position = -1
    count = 0
    dict = {}
    integrated = []
    product = 1
    readID = ''
    for row in reader:
        if not row[1].isdigit():
            continue
        if not row[3] == readID:
            if not readID == '':
                integrated.append([position, count, product])
                dict[readID] = integrated
                integrated = []
                position = -1
        readID = row[3]
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
                integrated.append([position, count, product])
            if int(row[1]) > position + 1:
                position = position + 1
                while position < int(row[1]):
                    integrated.append([position, 0, 0])
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
    return dict

def simplify_eventalign(f):
    cover = set()
    summary = []
    for line in f:
        line = line.split()
        summary.append(line)
        cover.add(line[0])
    print summary
    return cover, summary

def find_in_eventsummary(pair, anchor, eventsummary):
    result = []
    for read in eventsummary:
        if read[1][pair[0]:pair[1]] == anchor:
            result.append(read[0])
    return result

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
    integrated_eventalign = integrate_eventalign(sys.argv[3])
    print 'eventalign read'
    with file(sys.argv[1]) as f:
        for name, seq in fasta_read(f):
            contigs[name] = seq
    print 'fasta read'
    samData = []
    with file(sys.argv[2]) as f:
        for moleId, label, refName, pos, alignform, seq in get_SAM(f):
            samData.append((moleId, label, refName, pos, alignform, seq))
    print 'sam read'
    #integrated_eventalign = integrate_eventalign(sys.argv[3])
    #print 'eventalign read'
    for sam in samData:
        moleId, label, refName, pos, alignform, seq = sam
        print moleId
        if refName == 'KF278742.1':
            wseq, condicates = seq_warp_condicates(pos, seq, contigs[refName], alignform)
            print moleId
            #print(contigs[refName][pos - 1: pos + len(wseq)])
            #print(wseq)
            #print alignform
            for pair in condicates:
                if pair[1] > 2274:
                    continue
                print pair
                print contigs[refName][pair[0]: pair[1]]
                print wseq[pair[0]-pos+1: pair[1]-pos+1]
                temp = condicate_score(pair, integrated_eventalign[moleId])
                print temp
                #scores.append(temp)
            print '\n\n\n'
