import argparse
import itertools
import pickle
import os
import warnings

import shap
from sklearn.decomposition import PCA
import pandas as pd
import seaborn as sns
import numpy as np
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score, average_precision_score
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score
from sklearn.preprocessing import LabelEncoder, label_binarize, StandardScaler
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from xgboost import XGBClassifier
import matplotlib.pyplot as plt

def import_files(filename):  # import the csv file

	chunks = pd.read_csv(filename, sep='\t', low_memory=False, chunksize=500000)
	data = pd.concat(chunks, ignore_index=True)

	#data = change_sequence_id_column(data, 3)

	data['Sequence_id'] = data['Sequence_id'].str.replace('>', '', regex=False)

	if args.s != None:

		data['Sequence_id'] = data['Sequence_id'].str.split('|').str[args.s]

 	# Remove rows with empty values on the target feature
	data = data[data['Sequence_id'].notna() & (data['Sequence_id'].str.strip() != '')]

	return data


def drop_columns(data):
	le = LabelEncoder()
	X = data.drop(columns=['Sequence_id'])
	y = le.fit_transform(data['Sequence_id'])

	return X, y


def print_to_files(info):
	f_tex = open("performance_model.tex", "a")
	f_tsv = open("performance_model.tsv", "a")

	for elem in info:
		i = str(elem)
		f_tex.write(i + " & ")
		f_tsv.write(i + "\t")
	f_tex.write("\\\\\\hline\n")
	f_tsv.write("\n")

	f_tex.close()
	f_tsv.close()


def write_to_file(string_to_write):

	with open(file_tsv, 'a') as file:
		file.write(string_to_write)


def fit_and_predict(model, name, is_test):
	if is_test == True:

		scores, mean_score, std_dev = calculate_cv_scores(model, name)


		# Train the model
		model.fit(X_train, y_train)

		# Make predictions
		y_pred = model.predict(X_test)
		y_scores = model.predict_proba(X_test)

		# Binarize true labels using only the model's classes
		y_test_bin = label_binarize(y_test, classes=model.classes_)

		# Accuracy
		acc = accuracy_score(y_test, y_pred)

		# F1 Score (macro averages across classes)
		f1 = f1_score(y_test, y_pred, average='macro')

		# AUROC (macro-averaged across classes)
		auroc = roc_auc_score(y_test_bin, y_scores, multi_class='ovr', average='macro')

		# AUPRC (macro-averaged across classes)
		auprc = average_precision_score(y_test_bin, y_scores, average='macro')

		# Print results
		print(f"Accuracy: {acc:.4f}")
		print(f"F1 Score: {f1:.4f}")
		print(f"AUROC: {auroc:.4f}")
		print(f"AUPRC: {auprc:.4f}")

		string_to_write = (
				str(name) + "\t" +
				str(scores) + "\t" +
				str(mean_score) + "\t" +
				str(std_dev) + "\t" +
				str(acc) + "\t" +
				str(f1) + "\t" +
				str(auroc) + "\t" +
				str(auprc) + "\n"
		)

		write_to_file(string_to_write)

	else:
		print("Saving the " + name + "...")

		model.fit(X, Y)

		# save model
		filename = name + '.sav'
		pickle.dump(model, open(filename, 'wb'))


def balance_data(df):

	target = 'Sequence_id'


	for value in df[target].unique():

		if (df[target] == value).sum() <= 15:
			df = df[df[target] != value]

	# Count instances per label
	counts = df[target].value_counts()
	factor = 2 * counts.mean()

	for value in df[target].unique():

		if (df[target] == value).sum() > factor :

			instances_to_remove = (df[target] == value).sum() - factor
			indices_to_remove = df[df[target] == value].sample(n=int(instances_to_remove), random_state=42).index

			# Drop those rows
			df = df.drop(index=indices_to_remove)

	return df

def additional_removals (df):

	df = df[df['Sequence_id'] != "unknown"]
	df = df[df['Sequence_id'] != "RNA"]
	df = df[df['Sequence_id'] != "DNA"]
	df = df[df['Sequence_id'] != "dsDNA; ssDNA"]
	df = df[df['Sequence_id'] != "ssRNA"]


	return df


def PCA_feature_analysis():
	X_scaled = StandardScaler().fit_transform(X)

	# Run PCA
	pca = PCA(n_components=5)
	pca.fit(X_scaled)

	# Loadings: how much each original feature contributes to each PC
	loadings = pd.DataFrame(
		pca.components_.T,
		columns=[f'PC{i + 1}' for i in range(pca.n_components_)],
		index=X.columns
	)

	# Optional: plot loading magnitudes for PC1
	print(loadings['PC1'].abs().sort_values(ascending=False).head(10))


def calculate_cv_scores(model, name):
	# Print model name
	print("\n\nTesting the " + name + "...")

	scores = cross_val_score(model, X_train, y_train, cv=4)

	# Mean and standard deviation
	mean_score = np.mean(scores)
	std_dev = np.std(scores)


	print(f"Cross-validation scores: {scores}")
	print(f"Mean accuracy: {mean_score:.3f}")
	print(f"Standard deviation: {std_dev:.3f}")

	return scores, mean_score, std_dev


if __name__ == '__main__':

	warnings.filterwarnings("ignore")

	parser = argparse.ArgumentParser(description="Index", usage="Training and testing\n\npython3 genomeclass.py -f <input multi-FASTA file> -i <input (multi-)FASTA file>\n"
	                                                            "python3 genomeclass.py -t <input TSV file> -i <input (multi-)FASTA file>\n")

	parser.add_argument("-f", help="Input multi-FASTA file", type=str)
	parser.add_argument("-t", help="Input TSV file", type=str)
	parser.add_argument("-i", help="Input FASTA file containing the sequences to be classified", nargs="+", type=str, required=False)
	parser.add_argument("-s", help="Part of the Sequence_id that will become the target feature", type=int)
	parser.add_argument("-m", help="Machine learning model to be used. Default: RandomForestClassifier", type=str)
	parser.add_argument("-o", help="Options for the execution of the C file. Please surround the options with \"\"", type=str)
	parser.add_argument("-p", help="Add permutations of a certain number of characters", action='append', type=int, required=False)
	parser.add_argument("-a", help="Uses an auto balancer on the dataset", action='store_true', required=False)
	args = parser.parse_args()

	permutations_added = ""
	if args.p != None:
		chars = ['A', 'C', 'T', 'G']

		for i in args.p:

			permutations = [''.join(p) for p in itertools.product(chars, repeat=i)]
			for j in permutations:
				permutations_added += " -d "
				permutations_added += j

	if args.o != None:
		options = args.o + permutations_added
	else:
		options = "-s -g -c -e -m -t 7" + permutations_added

	if args.f is not None and os.path.exists(args.f) and args.i is not None and os.path.exists(args.i):

		print("Using " + args.f + " as the input file.\n")
		print("File to be classified: " + args.i[0] + "\n")
		os.system("make clean")
		os.system("make")
		print("\n./genomeclass -i " + args.f + " " + options + "\n\n")
		os.system("./genomeclass -i " + args.f + " " + options)
		dataset_name = "output.tsv"

		print("\n./genomeclass -i " + args.i + " " + options + " -o test.tsv\n\n")
		os.system("./genomeclass -i " + args.i + " " + options + " -o test.tsv")
		test_irl = "test.tsv"

	elif args.t is not None and os.path.exists(args.t):
		print("Using " + args.t + " as the input file.\n")
		dataset_name = args.t

	else:
		print("Invalid input files. Exiting.")
		exit(1)



	pd.set_option('display.max_columns', 30)

	file_tsv = "performance_genomeclass_testing.tsv"

	if os.path.exists(file_tsv):
		os.remove(file_tsv)


	write_to_file("Name_model\tScores_CV\tMean_score_CV\tStd_deviation\tAccuracy\tF1-score\tAUROC\tAUPRC\n")

	data = import_files(dataset_name)

	#print(data.dtypes)
	print(data.shape)
	#print(data.head(10))

	if args.a:
		data = balance_data(data)
	data = additional_removals(data)

	value_counts = data['Sequence_id'].value_counts()
	
	for seq_id, count in value_counts.items():
		print(f"{seq_id}: {count}")

	print(data['Sequence_id'].value_counts())

	print(f"Number of unique values: {data['Sequence_id'].nunique()}")
	data = data.loc[:, ~data.columns.str.startswith('Prob_sequence_')]


	# Start training process
	print(data.shape)

	# Separate features and target first
	X, Y = drop_columns(data)

	PCA_feature_analysis()

	# Train/test split (stratified)
	X_train, X_test, y_train, y_test = train_test_split(X, Y, stratify=Y, test_size=0.2, random_state=42)

	# Set models
	xgboost = XGBClassifier(random_state=42)
	random_forest = RandomForestClassifier(random_state=42)
	knn = KNeighborsClassifier()
	mlp = MLPClassifier(random_state=42)

	# Get performance in the cross validation and train the models on the entire training set

	fit_and_predict(xgboost, "XGBClassifier", True)
	fit_and_predict(random_forest, "RandomForestClassifier", True)
	fit_and_predict(knn, "KNeighborsClassifier", True)
	fit_and_predict(mlp, "MLPClassifier", True)


	# Save the models

	'''fit_and_predict(xgboost, "XGBClassifier", False)
	fit_and_predict(random_forest, "RandomForestClassifier", False)
	fit_and_predict(knn, "KNeighborsClassifier", False)
	fit_and_predict(mlp, "MLPClassifier", False)'''



	'''if args.i is not None:
		X_train = X
		y_train = Y

		test_data = import_files(args.i)

		X_test, y_test = drop_columns(test_data)

		#xgboost = XGBClassifier(random_state=42)
		#fit_and_predict(xgboost, "XGBClassifier", True)

		random_forest = RandomForestClassifier(random_state=42)
		fit_and_predict(random_forest, "RandomForestClassifier", True)

		knn = KNeighborsClassifier()
		fit_and_predict(knn, "KNeighborsClassifier", True)

		mlp = MLPClassifier(random_state=42)
		fit_and_predict(mlp, "MLPClassifier", True)
'''

