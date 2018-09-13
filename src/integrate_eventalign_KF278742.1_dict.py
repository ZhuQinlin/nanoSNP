import sys
import string
import csv
from itertools import imap, ifilter
from scipy.stats import norm

def score(data, mean, stdv):
    x = (data-mean)/stdv
    delta = norm.ppf(0.9)
    return norm.cdf(x+delta) - norm.cdf(x-delta)

reader = csv.reader(open(sys.argv[1]), delimiter = '\t')
position = -1
count = 0
dict = {}
integrated = []
product = 1
readID = ''
for row in reader:
    if not row[1].isdigit():
        continue
    if not row[0] == 'KF278742.1':
        break
    if not cmp(row[3], readID) == 0:
        if not readID == '':
            integrated.append([position, count, product])
            dict[readID] = integrated
            print readID
            #for item in integrated:
            #    print item
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

f = open('integrated_eventalign_KF_dict.txt', 'w')
f.write(str(dict))
f.close()
