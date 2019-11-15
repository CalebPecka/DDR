import codecs
import csv
import pandas
import sys, operator

def getRefExpLevel(lowerTh, upperTh, outputfile):
    with codecs.open("final_out.csv", "r", encoding='utf-8', errors='ignore') as f:
        with open(outputfile, 'w') as csvfile:
          reader = csv.DictReader(f)
          writer = csv.writer(csvfile)
          writer.writerow(["genes", "cov", "mean", "std", "MFC", "covMFC"])
          for row in reader:
            #try:

                ens = []
                if lowerTh <= pandas.to_numeric(row["mean"]) <= upperTh:
                  covMFC = float(row["MFC"])*float(row["cov"])
                  Account = [row["genes"], row["cov"], row["mean"], row["std"], row["MFC"], covMFC]
                #print(row["Gene_ID"] + " " + mg.getgene(row["Gene_ID"])['symbol'])
                  writer.writerow(Account)

def findTop(resultFile):
    df  = pandas.read_csv(resultFile)
    sortedlist = df.sort_values(by='covMFC')
    for index, row in sortedlist.iterrows():
        if row[0] in lines:
            print(row[0])
            break

text_file = open("data/housekeeping_symbol.txt", "r")
lines = text_file.read().splitlines()

refRange = open("data/referenceRanges.csv", "r")
ULrange = refRange.read().split()
del ULrange[:4]
ULrange = [float(i) for i in ULrange]

lows = []
ups = []

for i in range(4):
    lows.append(int(ULrange[i] - ULrange[i]*0.2))
    ups.append(int(ULrange[i] + ULrange[i]*0.2) + 1)

outfiles = ['ref1.csv','ref2.csv','ref3.csv','ref4.csv']

for i in range(len(outfiles)):
    getRefExpLevel(lows[i],ups[i],outfiles[i])

for outfile in outfiles:
    findTop(outfile)

sys.exit(0)
sys.exit(1)
