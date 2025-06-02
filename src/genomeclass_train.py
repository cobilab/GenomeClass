import pickle
import time
import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.ensemble import GradientBoostingRegressor, RandomForestClassifier
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, \
	mean_absolute_percentage_error, accuracy_score
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.neighbors import KNeighborsRegressor, KNeighborsClassifier
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import LabelEncoder


def import_files(filename, name_pickle):  # import the csv file

	chunks = pd.read_csv(filename, sep='\t', low_memory=False, chunksize=500000)
	data = pd.concat(chunks, ignore_index=True)

	print(data)

	#data.to_pickle(name_pickle)

	return data


def print_info(data):  # prints data information

	# check dimensions
	data.shape

	# check the info on the columns - no null values
	data.info()

	# Summary of the numerical attributes
	data.describe()


def correlation(data):
	# Delete less relevant features based on the correlation with the output variable
	data_copy = data
	data_copy.drop("Sequence id", axis=1, inplace=True)
	cor = data_copy.corr()  # calculate correlations

	sns.set(font_scale=1.3)

	# Correlation graph
	plt.figure(figsize=(22, 22))
	sns.heatmap(cor, annot=True, cmap=sns.diverging_palette(20, 220, n=200), vmin=-1, vmax=1)
	plt.savefig('correlation.pdf')
	plt.savefig('correlation.jpg')
	plt.show()


# Correlation with output variable
# list_columns_dropped = remove_low_correlation(data, cor_target)

# return list_columns_dropped




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


def cross_validation_MLPRegressor(X_train, y_train, y_test):
	param_activation = ['tanh', 'relu']
	param_solver = ['sgd', 'adam']
	param_alpha = [0.0001, 0.001, 0.01, 0.1]
	param_learning_rate = ['constant', 'adaptive']
	param_hidden_layer_sizes = [(50, 50, 50), (50, 100, 50), (100,), (5, 5, 5), (50, 50), (5, 5), (2, 5, 2)]

	for activation in param_activation:

		for solver in param_solver:

			for alpha in param_alpha:

				for learning_rate in param_learning_rate:

					for hidden_layer_sizes in param_hidden_layer_sizes:
						model = MLPRegressor(hidden_layer_sizes=hidden_layer_sizes, activation=activation,
											 solver=solver, alpha=alpha, learning_rate=learning_rate, random_state=42)

						model.fit(X_train, y_train)
						y_pred = model.predict(X_test)

						# Evaluate the model's performance
						mse = mean_squared_error(y_test, y_pred)
						r2 = r2_score(y_test, y_pred)
						mae = mean_absolute_error(y_test, y_pred)
						mape = mean_absolute_percentage_error(y_test, y_pred)
						print(hidden_layer_sizes, activation, solver, alpha, learning_rate)
						print(f"Mean squared error: {mse:.2f}")
						print(f"R-squared: {r2:.2f}")
						print(f"Mean absolute error: {mae:.2f}")
						print(f"Mean absolute percentage error: {mape:.2f}")

						info = [hidden_layer_sizes, activation, solver, alpha, learning_rate, mse, r2, mae, mape]

						print_to_files(info)


def cross_validation_MLPRegressor_v2(X_train, y_train, y_test):
	param_grid = {
		'hidden_layer_sizes': [(20,), (15, 30), (15, 15), (10,), (10, 10), (20, 10), (10, 20), (20, 20)],
		# 'hidden_layer_sizes': [(20, 20), (20, 20, 20), (20, 20, 20, 20)],
		# 'hidden_layer_sizes': [(20,), (19,), (18,), (21,), (22,), (20, 5), (20, 20), (20, 30)],
		'activation': ['relu'],
		'solver': ['adam', 'sgd'],
		'alpha': [0.01, 0.005, 0.05],
		'learning_rate': ['constant', 'adaptive']
	}

	model = MLPRegressor(random_state=42)

	cv = GridSearchCV(model, param_grid, n_jobs=6, verbose=10, cv=3)
	cv.fit(X_train, y_train)
	print(cv.best_estimator_)
	print(cv.best_score_)
	print(cv.best_params_)

	return cv


def cross_validation_GradientBoostingRegression(X_train, y_train, y_test):
	param_grid = {
		'loss': ["squared_error", "absolute_error"],
		'learning_rate': [0.1, 0.2, 0.3],
		'criterion': ["friedman_mse"],
		'n_estimators': [15, 30, 50],
		'min_samples_split': [2, 4]
	}

	model = GradientBoostingRegressor(random_state=42)

	cv = GridSearchCV(model, param_grid, n_jobs=6, verbose=10, cv=3)
	cv.fit(X_train, y_train)
	print(cv.best_estimator_)
	print(cv.best_score_)
	print(cv.best_params_)

	return cv


def cross_validation_NNR(X_train, y_train, y_test):
	param_grid = {
		'n_neighbors': [1, 3, 5, 7, 9, 15],
		'weights': ['uniform', 'distance'],
		'algorithm': ['brute', 'auto'],
	}

	model = KNeighborsRegressor()

	cv = GridSearchCV(model, param_grid, n_jobs=-1, verbose=10, cv=2, error_score='raise')
	cv.fit(X_train, y_train)
	print(cv.best_estimator_)
	print(cv.best_score_)
	print(cv.best_params_)

	return cv


def cross_validation_LinearRegression(X_train, y_train, y_test):
	param_grid = {
		'fit_intercept': [True, False],
		'n_jobs': [-1],
		'C': [1, 5, 10, 50, 100],
		'gamma': ["scale", "auto"]
	}

	model = LinearRegression(random_state=42)

	cv = GridSearchCV(model, param_grid, n_jobs=-1, verbose=1, cv=2)
	cv.fit(X_train, y_train)
	print(cv.best_estimator_)
	print(cv.best_score_)
	print(cv.best_params_)

	return cv


def generate_plots(data):
	data.hist(bins=50, figsize=(20, 15))
	plt.show()


def fit_and_predict(model, name, is_test):
	if is_test == True:
		print("Testing the " + name + "...")

		print(data.columns)

		model.fit(X_train, y_train)

		y_pred = model.predict(X_test)

		# Evaluate the model's performance
		mse = mean_squared_error(y_test, y_pred)
		r2 = r2_score(y_test, y_pred)
		mae = mean_absolute_error(y_test, y_pred)
		mape = mean_absolute_percentage_error(y_test, y_pred)
		accuracy = accuracy_score(y_test, y_pred)
		print(f"Mean squared error: {mse:.4f}")
		print(f"Mean absolute error: {mae:.4f}")
		print(f"R-squared: {r2:.4f}")
		print(f"Mean absolute percentage error: {mape:.4f}")
		print(f"Accuracy: {accuracy:.4f}")

	else:
		print("Saving the " + name + "...")

		model.fit(X, Y)

		# save model
		filename = name + '.sav'
		pickle.dump(model, open(filename, 'wb'))


if __name__ == '__main__':

	pd.set_option('display.max_columns', 30)

	filename = "performance_model.tex"
	file_tsv = "performance_model.tsv"
	if os.path.exists(filename):
		os.remove(filename)
	if os.path.exists(file_tsv):
		os.remove(file_tsv)

	f_tsv = open(file_tsv, "w")
	f_tsv.write("hidden_layer_sizes\tactivation\tsolver\talpha\tlearning_rate\tmse\tr2\tmae\tmape\n")
	f_tsv.close()

	init_time = time.perf_counter()

	base_dataset_name = "output"

	# if pickle file exists read from there as it is faster
	if os.path.exists(base_dataset_name + '.pickle'):
		data = pd.read_pickle(base_dataset_name + '.pickle')
	else:
		data = import_files(base_dataset_name + '.tsv', base_dataset_name + '.pickle')

	print("TIME ->", time.perf_counter() - init_time)
	print(data.shape)

	print(data.dtypes)
	print(data.shape)

	# correlation(data)

	X, Y = drop_columns(data)
	X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

	# get parameters
	# print("Starting MLPRegressor")
	# model_mlp = cross_validation_MLPRegressor_v2(X_train, y_train, X_test)

	# print("Starting gradientboosting")
	# cross_validation_GradientBoostingRegression(X_train, y_train, X_test)

	# train and save models
	mlp_model = RandomForestClassifier(random_state=42)
	fit_and_predict(mlp_model, "mlp_model", True)
	#fit_and_predict(mlp_model, "mlp_model", False)

'''	gbr_model = GradientBoostingRegressor(learning_rate=0.3, min_samples_split=4, n_estimators=50, random_state=42)
	fit_and_predict(gbr_model, "gbr_model", True)
	fit_and_predict(gbr_model, "gbr_model", False)'''
