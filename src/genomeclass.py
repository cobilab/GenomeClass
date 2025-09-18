import argparse
import csv
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

def import_files(filename, sep):  # import the csv file

	#Read the tsv file
	chunks = pd.read_csv(filename, sep='\t', low_memory=False, chunksize=500000)
	data = pd.concat(chunks, ignore_index=True)

	#data = change_sequence_id_column(data, 3)

	data['Sequence_id'] = data['Sequence_id'].str.replace('>', '', regex=False)

	if sep != None:

		data['Sequence_id'] = data['Sequence_id'].str.split('|').str[sep]

 	# Remove rows with empty values on the target feature
	data = data[data['Sequence_id'].notna() & (data['Sequence_id'].str.strip() != '')]

	return data


def drop_columns(data):
	le = LabelEncoder()
	X = data.drop(columns=['Sequence_id'])

	try:
		y = le.fit_transform(data['Sequence_id'])
	except:
		y = None
		le = None

	return X, y, le


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


def fit_and_predict(model, name, is_test, X_train, y_train, X_test = None, y_test = None):

	if is_test == True:

		scores, mean_score, std_dev = calculate_cv_scores(model, name, X_train, y_train)


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

		ml_model = model.fit(X_train, y_train)

		# save model
		filename = name + '.sav'
		pickle.dump(model, open(filename, 'wb'))
		
		return ml_model


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

def min_balancing(df):

	target = 'Sequence_id'

	for value in df[target].unique():

		if (df[target] == value).sum() <= 2:
			df = df[df[target] != value]

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


def calculate_cv_scores(model, name, X_train, y_train):
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

def write_predictions(X_column, predictions, name_file):
	results_df = pd.DataFrame(X_column)
	results_df.insert(0, 'predicted_value', predictions)

	# --- Save to TSV file ---
	results_df.to_csv(name_file, sep='\t', index=False)

	print("Predictions saved to '" + name_file + "'")




if __name__ == '__main__':

	warnings.filterwarnings("ignore")

	parser = argparse.ArgumentParser(description="Index", usage="Training and testing\n\npython3 genomeclass.py -f <input multi-FASTA file> -i <input (multi-)FASTA file>\n"
																"python3 genomeclass.py -t <input TSV file> -i <input (multi-)FASTA file>\n")

	parser.add_argument("-tf", help="Input training multi-FASTA file", type=str)
	parser.add_argument("-tt", help="Input training TSV file", type=str)
	parser.add_argument("-cf", help="Input FASTA file containing the sequences to be classified", type=str, required=False)
	parser.add_argument("-ct", help="Input TSV file containing the sequences to be classified", type=str,
						required=False)
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

	if args.tf is not None and os.path.exists(args.tf) :

		print("Using " + args.tf + " as the training file and " + args.cf + " as the file to be classified.\n")
		#print("File to be classified: " + args.cf[0] + "\n")
		os.system("make clean")
		os.system("make")

		print("Analysing the training file.\n")
		dataset_name = "output_test.tsv"
		print("\n./genomeclass -i " + args.tf + " " + options + " -o " + dataset_name + "\n\n")
		os.system("./genomeclass -i " + args.tf + " " + options + " -o " + dataset_name)

	elif args.tt is not None and os.path.exists(args.tt):
		print("Using " + args.tt + " as the input file.\n")
		dataset_name = args.tt

	else:
		print("Invalid input files. Exiting.")
		exit(1)

	# Create the file that will hold the statistics for the training set
	file_tsv = "performance_genomeclass_testing.tsv"

	if os.path.exists(file_tsv):
		os.remove(file_tsv)

	write_to_file("S" + str(args.s) + "_Name_model\tScores_CV\tMean_score_CV\tStd_deviation\tAccuracy\tF1-score\tAUROC\tAUPRC\n")

	# Import training data
	data_tds = import_files(dataset_name, args.s)

	# Minimal balancing, removes classes with 2 or less instances - avoids execution errors
	data_tds = min_balancing(data_tds)
	
	# Balances the data if wanted
	if args.a:
		data_tds = balance_data(data_tds)
	data_tds = additional_removals(data_tds)
	
	# Removes the columns that begin with Prob_sequence as they hold information similar to those that begin with Avg_distance
	data_tds = data_tds.loc[:, ~data_tds.columns.str.startswith('Prob_sequence_')]
	
	# -> Train the models on the training dataset
	print(data_tds.shape)

	# Separate features and target first
	X_tds, Y_tds, le = drop_columns(data_tds)

	# Train/test split (stratified)
	X_train_tds, X_test_tds, y_train_tds, y_test_tds = train_test_split(X_tds, Y_tds, stratify=Y_tds, test_size=0.2, random_state=42)

	# Set models
	xgboost = XGBClassifier(random_state=42)
	random_forest = RandomForestClassifier(random_state=42)
	knn = KNeighborsClassifier()
	mlp = MLPClassifier(random_state=42)

	# Get performance in the cross validation and train the models on the entire training set
	fit_and_predict(xgboost, "XGBClassifier", True, X_train_tds, y_train_tds, X_test_tds, y_test_tds)
	fit_and_predict(random_forest, "RandomForestClassifier", True, X_train_tds, y_train_tds, X_test_tds, y_test_tds)
	fit_and_predict(knn, "KNeighborsClassifier", True, X_train_tds, y_train_tds, X_test_tds, y_test_tds)
	fit_and_predict(mlp, "MLPClassifier", True, X_train_tds, y_train_tds, X_test_tds, y_test_tds)


	# Save the models
	xgbmodel_trained = fit_and_predict(xgboost, "XGBClassifier", False, X_tds, Y_tds)
	rfmodel_trained = fit_and_predict(random_forest, "RandomForestClassifier", False, X_tds, Y_tds)
	knnmodel_trained = fit_and_predict(knn, "KNeighborsClassifier", False, X_tds, Y_tds)
	mlpmodel_trained = fit_and_predict(mlp, "MLPClassifier", False, X_tds, Y_tds)


	# Begining classification of outside sequences

	test_irl = "test.tsv"

	# If there is real data to be classified - FASTA format
	if args.cf is not None and os.path.exists(args.cf):

		print("Analysing the file to be classified.\n")
		print("\n./genomeclass -i " + args.cf + " " + options + " -o " + test_irl + "\n\n")
		os.system("./genomeclass -i " + args.cf + " " + options + " -o " + test_irl)

	# If there is real data to be classified - TSV format
	elif args.ct is not None and os.path.exists(args.ct):
		print("Using " + args.ct + " as the input file.\n")
		test_irl = args.ct

	else:
		print("File to be classified is not available. Exiting...")
		exit(1)

	# Import data to be classified
	data_clf = import_files(test_irl, args.s)

	# Removes the columns that begin with Prob_sequence as they hold information similar to those that begin with Avg_distance
	data_clf = data_clf.loc[:, ~data_clf.columns.str.startswith('Prob_sequence_')]

	# -> Train the models on the training dataset
	print(data_clf.shape)

	# Separate features and target first
	X_clf, Y_clf, le_clf = drop_columns(data_clf)

	# Predict the classifications
	predictions_xgb = xgbmodel_trained.predict(X_clf)
	predictions_rf = rfmodel_trained.predict(X_clf)
	predictions_knn = knnmodel_trained.predict(X_clf)
	predictions_mlp = mlpmodel_trained.predict(X_clf)

	predicted_labels_xgb = le.inverse_transform(predictions_xgb)
	predicted_labels_rf = le.inverse_transform(predictions_rf)
	predicted_labels_knn = le.inverse_transform(predictions_knn)
	predicted_labels_mlp = le.inverse_transform(predictions_mlp)

	# Write predictions to a TSV file
	write_predictions(data_clf, predicted_labels_xgb, "predictions_xgboost.tsv")
	write_predictions(data_clf, predicted_labels_rf, "predictions_randforest.tsv")
	write_predictions(data_clf, predicted_labels_knn, "predictions_knn.tsv")
	write_predictions(data_clf, predicted_labels_mlp, "predictions_mlp.tsv")




