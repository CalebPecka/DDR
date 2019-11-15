import pandas as pd
from sklearn.model_selection import train_test_split
import numpy as np
from sklearn import datasets
from sklearn import svm
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, GradientBoostingClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import LinearSVC, NuSVC
from sklearn.naive_bayes import GaussianNB
from imblearn.over_sampling import SMOTE

from sklearn.model_selection import GridSearchCV

import random
import csv
import sys

###########MEAN###############
def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

########The SELECTION of BIOMARKERS#############UNUSED
#def biomarker_selection(data, size):
#    top_biomarker = list(data.iloc[:,0:size]) + ['label']
#    data_marker = data[top_biomarker]
#    return data_marker

###classifer Methods########UNUSED
#def classifer (train_X, train_Y, test_X,test_Y, method):
#   cl = method.fit(train_X, train_Y)
#   score = cl.score(test_X,test_Y)
#   return score

#####cross validation #####################UNUSED
#def cross_validation(feature, annotation, method,test_size, cv):
#   accuracy = []
#   for x in range(cv):
#      train_X, test_X, train_Y, test_Y = train_test_split(feature, annotation, test_size = test_size)
#      score = classifer (train_X,train_Y, test_X,test_Y, method = method)
#      accuracy.append(score)
#   return accuracy
##############Sample Size Selection ################UNUSED
#def sample_selection(data,len1,size):
#    idx = np.random.randint(len(data)-len1, size = size)
#    idx1 = np.random.randint(len1, size = size)
#    s1 = data.ix[idx]
#    s2 = data.ix[idx1+len(data)-len1]
#    s = pd.concat([s1, s2])
#    return s

############Evaluation Metrics##############
###classifer Methods########
def classify (train_X, train_Y, test_X,test_Y, method):
   cl = method.fit(train_X, train_Y)
   return cl



#def accuracy(y_pred,y_true):
#    accuracy = accuracy_score(y_true, y_pred, average='macro')
#    return accuracy

#def precision(y_pred,y_true):
#    precision = precision_score(y_true, y_pred, average='macro')
#    return precision

#def cross_recall(feature, annotation, method,test_size, cv):#UNUSED
#   recall = []
#   for x in range(cv):
#      train_X, test_X, train_Y, test_Y = train_test_split(feature, annotation, test_size = test_size)
#      clf = classify (train_X,train_Y, test_X,test_Y, method = method)
#      y_pred = clf.predict(test_X)
#      rec = recall_score(test_Y, y_pred, average='macro')
#      recall.append(rec)
#   return recall
#
#def cross_accuracy(feature, annotation, method,test_size, cv):#UNUSED
#   accuracy = []
#   for x in range(cv):
#      train_X, test_X, train_Y, test_Y = train_test_split(feature, annotation, test_size = test_size)
#      clf = classify (train_X,train_Y, test_X,test_Y, method = method)
#      y_pred = clf.predict(test_X)
#      rec = accuracy_score(test_Y, y_pred)
#      accuracy.append(rec)
#   return accuracy

pattern = 'feature_and_labels.csv'
data_r = pd.read_csv(pattern)
columns = list(data_r)
titles = ["Gene1", "Gene2", "Accuracy", "Recall", "Precision", "f1 score", "MU acc", "WT acc"]
resultCollection = [titles]
group1SampleSize = int(sys.argv[1])
group2SampleSize = len(data_r) - group1SampleSize

for i in columns[:-2]:

   for j in columns[columns.index(i)+1: -1]:

      accuracy = []
      recall = []
      precision = []
      f1score = []
      MU_acc = []
      WT_acc = []

      #EACH y CREATES A RANDOM SAMPLE
      for y in range(1):
         #healthySample = []
         #data_r = biomarker_selection(data,18)
         #if group1SampleSize > group2SampleSize:
         #   healthySample = random.sample(range(0, group1SampleSize - 1), group1SampleSize - group2SampleSize)
         #else:
         #   healthySample = random.sample(range(group1SampleSize, len(data_r) - 1), group2SampleSize - group1SampleSize)

         #data_temp = data_r.drop(data_r.index[healthySample])

         #data = data_temp[[i, 'label']]
         data = data_r[[i, j, 'label']]

         MU = data['label'] == "group1"
         WT = data['label'] == "group2"

         y_MU = data[MU] .label
         X_MU = data[MU].drop('label', axis=1)

         y_WT = data[WT] .label
         X_WT = data[WT].drop('label', axis=1)

         for x in range(100):
            X_MU_train, X_MU_test, y_MU_train, y_MU_test = train_test_split(X_MU, y_MU,test_size=0.3, random_state=12)
            X_WT_train, X_WT_test, y_WT_train, y_WT_test = train_test_split(X_WT, y_WT,test_size=0.3, random_state=12)

            #TRAINING
            y_train = pd.concat([y_MU_train, y_WT_train])
            X_train = pd.concat([X_MU_train, X_WT_train])
            y_test = pd.concat([y_MU_test, y_WT_test])
            X_test = pd.concat([X_MU_test, X_WT_test])
            
            if group1SampleSize == group2SampleSize:
               X_train_res = X_train
               y_train_res = y_train
            else:
              sm = SMOTE(random_state=12, ratio = 1.0)
              X_train_res, y_train_res = sm.fit_sample(X_train, y_train)

            #DEFINING PARAMETER RANGE FOR CLF OPTIMIZATION
            param_grid = {'C': [0.1, 10, 100, 1000], 'gamma': [1, 0.1, 0.01, 0.001, 0.0001], 'kernel': ['rbf']}
            grid = GridSearchCV(svm.SVC(), param_grid, refit = True, verbose = False, cv = 5)
            grid.fit(X_train, y_train)

            #CLF RESULTS
            clf = svm.SVC(kernel='rbf', C=grid.best_params_['C'], gamma = grid.best_params_['gamma']).fit(X_train_res, y_train_res)
            y_pred = clf.predict(X_test)
            acc = accuracy_score(y_test, y_pred)
            rec = recall_score(y_test, y_pred, average='macro')
            prec = precision_score(y_test, y_pred, average='macro')
            f1 = f1_score(y_test, y_pred, average='macro')

            #APPEND THESE VALUES
            accuracy.append(acc)
            recall.append(rec)
            precision.append(prec)
            f1score.append(f1)

            MU_pred = y_pred[0:len(y_MU_test)]
            WT_pred = y_pred[len(y_MU_test):len(y_test)]
            #print(MU_pred)
            MU_corr = MU_pred == "group1"
            #print(MU_pred[MU_corr])
            WT_corr = WT_pred == "group2"
            MU_accuracy = len(MU_pred[MU_corr])/len(y_MU_test)
            WT_accuracy = len(WT_pred[WT_corr])/len(y_WT_test)
            MU_acc.append(MU_accuracy)
            WT_acc.append(WT_accuracy)
            #print(y_test)

      resultCollection.append([i, j, mean(accuracy), mean(recall), mean(precision), mean(f1score), mean(MU_acc), mean(WT_acc)])

with open('combinationGeneOutput.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(resultCollection)
csvFile.close()
