import pandas as pd
from sklearn.model_selection import train_test_split
import numpy as np
from sklearn import datasets
from sklearn import svm
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from imblearn.over_sampling import SMOTE

from sklearn.model_selection import GridSearchCV

import random
import csv
import sys

###########MEAN###############
def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

############Evaluation Metrics##############
###classifer Methods########
def classify (train_X, train_Y, test_X,test_Y, method):
	cl = method.fit(train_X, train_Y)
	return cl

pattern = sys.argv[1]
data_r = pd.read_csv(pattern)
columns = list(data_r)
titles = ["Gene", "Accuracy", "Recall", "Precision", "f1 score", "MU acc", "WT acc"]
resultCollection = [titles]
group1SampleSize = int(sys.argv[2])
group2SampleSize = len(data_r) - group1SampleSize

for i in columns[:-1]:

	accuracy = []
	recall = []
	precision = []
	f1score = []
	MU_acc = []
	WT_acc = []

	#healthySample = []
	#data_r = biomarker_selection(data,18)
	#if group1SampleSize > group2SampleSize:
	#	healthySample = random.sample(range(0, group1SampleSize - 1), group1SampleSize - group2SampleSize)
	#else:
	#	healthySample = random.sample(range(group1SampleSize, len(data_r) - 1), group2SampleSize - group1SampleSize)
	#data_temp = data_r.drop(data_r.index[healthySample])
	#data = data_temp[[i, 'label']]
	
	data = data_r[[i, 'label']]

	MU = data['label'] == "group1"
	WT = data['label'] == "group2"

	y_MU = data[MU] .label
	X_MU = data[MU].drop('label', axis=1)

	y_WT = data[WT] .label
	X_WT = data[WT].drop('label', axis=1)

	####################################
	#FAILED HYPERPARAMETER OPTIMIZATION
	####################################
	#N_ROWS = 10000
	#STATIC_PARAMS = {'boosting': 'gbdt', 'objective':'binary', 'metric': 'auc', 'num_threads': 12,}
	#N_CALLS = 10

	#space = [skopt.space.Real(0.01, 0.5, name='learning_rate', prior='log-uniform'), skopt.space.Integer(1, 30, name='max_depth'), skopt.space.Integer(1, 100, name='num_leaves'), skopt.space.Integer(10, 1000, name='min_data_in_leaf'), skopt.space.Real(0.1, 1.0, name='feature_fraction', prior='uniform'), skopt.space.Real(0.1, 1.0, name='subsample', prior='uniform'),]

	#def objective():
		#all_params = {**STATIC_PARAMS}
		#return -1.0 * train_evaluate(X_MU, y_MU, all_params)
		
	#results = skopt.forest_minimize(objective(), space, n_calls = N_CALLS)
	#print('Best VAlidation AUC: {}'.format(-1.0 * results.fun))
	#print('Best Params: {}',format(results.x))

	for x in range(100):
		X_MU_train, X_MU_test, y_MU_train, y_MU_test = train_test_split(X_MU, y_MU,test_size=0.3, random_state=12)
		X_WT_train, X_WT_test, y_WT_train, y_WT_test = train_test_split(X_WT, y_WT,test_size=0.3,random_state=12)
			
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
		param_grid = {'C': [0.1, 1, 10, 100, 1000], 'gamma': [1, 0.1, 0.01, 0.001, 0.0001], 'kernel': ['rbf']}
		grid = GridSearchCV(svm.SVC(), param_grid, refit = True, verbose = False, cv = 5)
		grid = grid.fit(X_train, y_train)	

		#CLF RESULTS
		clf = svm.SVC(kernel = 'rbf', C = grid.best_params_['C'], gamma = grid.best_params_['gamma']).fit(X_train_res, y_train_res)
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

	resultCollection.append([i, mean(accuracy), mean(recall), mean(precision), mean(f1score), mean(MU_acc), mean(WT_acc)])

with open('singleGeneOutput.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(resultCollection)
csvFile.close()
