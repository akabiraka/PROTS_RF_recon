from sklearn.ensemble import RandomForestRegressor
from sklearn.datasets import make_regression
from scipy import stats

X_train, y_train = make_regression(n_samples=500, n_features=27, n_informative=20, n_targets=1, shuffle=True)
X_test, y_test = make_regression(n_samples=100, n_features=27, n_informative=20, n_targets=1, shuffle=True)
print(X_train.shape, X_train[0].dtype, y_train.shape, y_train[0].dtype)

regr = RandomForestRegressor()

regr.fit(X_train, y_train)

y_pred = regr.predict(X_test)
r, p = stats.pearsonr(y_test, y_pred)

print(r, p)