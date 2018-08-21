import sys
import string
from itertools import imap, ifilter

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

def get_events(f):
    pass

def get_eventAlign(f):
    pass

def candicate_socre(candicate_range):
    pass

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
    for sam in samData:
        label, refName, pos, alignform, seq = sam
        if refName == 'KF278742.1':
            wseq = seq_warping(seq, alignform)
            condicates = seq_condicate(alignform)
            print(contigs[refName][pos - 1: pos + len(wseq)])
            print(wseq)
            print(condicates)
            break
