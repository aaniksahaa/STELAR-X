#!/usr/bin/env python3
"""
poly_fit_predict.py

Usage:
    python poly_fit_predict.py
    python poly_fit_predict.py --csv data.csv --xcol x --ycol y --max-degree 6

Requirements:
    pip install numpy pandas matplotlib scikit-learn

The script:
 - Loads data (default sample is embedded)
 - Fits polynomial models of degrees 1..max_degree
 - Reports training R² and CV RMSE
 - Chooses best degree by lowest CV RMSE
 - Re-fits chosen model on entire dataset
 - Prints polynomial formula and example predictions
 - Saves plot to ./polynomial_fits.png
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold, cross_val_score
from sklearn.metrics import mean_squared_error, r2_score
import textwrap
import sys

def load_data_from_args(args):
    if args.csv:
        df = pd.read_csv(args.csv)
        if args.xcol not in df.columns or args.ycol not in df.columns:
            raise ValueError(f"CSV must contain columns '{args.xcol}' and '{args.ycol}'")
        return df[[args.xcol, args.ycol]].rename(columns={args.xcol: "x", args.ycol: "y"})
    else:
        # Default data you provided
        data = {
            "x": [1, 2, 5, 10, 15, 20, 25, 30],
            "y": [2.219, 4.886, 99.781, 940.196, 2994.442, 8141.740, 15584.114, 27279.418],
        }
        # data = {
        #     "x": [15, 20, 25],
        #     "y": [18950.02, 27698.211, 35323],
        # }
        # data = {
        #     "x": [10, 20, 25, 50],
        #     "y": [940.196, 1918.991, 2020.297, 3682.449],
        # }
        return pd.DataFrame(data)

def fit_and_evaluate(X, y, max_degree=5, cv_splits=5):
    results = []
    cv = KFold(n_splits=cv_splits, shuffle=True, random_state=0)
    for deg in range(1, max_degree + 1):
        model = make_pipeline(PolynomialFeatures(degree=deg, include_bias=False),
                              LinearRegression())
        # Fit on full data for train metrics
        model.fit(X, y)
        y_pred_train = model.predict(X)
        train_r2 = r2_score(y, y_pred_train)
        #train_rmse = mean_squared_error(y, y_pred_train, squared=False)
        train_rmse = np.sqrt(mean_squared_error(y, y_pred_train))


        # CV RMSE (negative MSE -> convert)
        neg_mse_scores = cross_val_score(model, X, y, scoring="neg_mean_squared_error", cv=cv)
        cv_rmse = np.sqrt(-neg_mse_scores.mean())

        # Extract coefficients from fitted linear reg step
        lin = model.named_steps['linearregression']
        # Coeffs correspond to polynomial feature order when include_bias=False
        coef = lin.coef_.copy()
        intercept = lin.intercept_.item() if np.size(lin.intercept_) == 1 else lin.intercept_

        results.append({
            "degree": deg,
            "model": model,
            "train_r2": train_r2,
            "train_rmse": train_rmse,
            "cv_rmse": cv_rmse,
            "coef": coef,
            "intercept": intercept
        })

    return results

def poly_formula_string(intercept, coefs):
    # coefs is array for x^1, x^2, ...
    parts = [f"{intercept:.6g}"]
    for i, c in enumerate(coefs, start=1):
        parts.append(f"{c:.6g}*x^{i}")
    # join with + / - nicely
    expr = " + ".join(parts)
    # cosmetic: replace '+ -' with '- '
    expr = expr.replace("+ -", "- ")
    return expr

def predict_from_model(model, x_vals):
    arr = np.array(x_vals).reshape(-1, 1)
    return model.predict(arr)

def bootstrap_predict_interval(X, y, degree, x_values, n_bootstrap=1000, alpha=0.05, random_state=1):
    # Simple residual bootstrap: fit on full data, sample residuals
    rng = np.random.RandomState(random_state)
    base_model = make_pipeline(PolynomialFeatures(degree=degree, include_bias=False), LinearRegression())
    base_model.fit(X, y)
    y_hat = base_model.predict(X)
    residuals = y - y_hat
    preds_boot = np.zeros((n_bootstrap, len(x_values)))
    for i in range(n_bootstrap):
        resampled = rng.choice(residuals, size=len(residuals), replace=True)
        y_star = y_hat + resampled
        m = make_pipeline(PolynomialFeatures(degree=degree, include_bias=False), LinearRegression())
        m.fit(X, y_star)
        preds_boot[i, :] = m.predict(np.array(x_values).reshape(-1, 1))
    lower = np.percentile(preds_boot, 100 * (alpha/2), axis=0)
    upper = np.percentile(preds_boot, 100 * (1 - alpha/2), axis=0)
    median = np.median(preds_boot, axis=0)
    return median, lower, upper

def plot_results(df, results, best_degree, out_path="polynomial_fits.png", show=True):
    X = df["x"].values.reshape(-1, 1)
    y = df["y"].values
    x_min, x_max = df["x"].min(), df["x"].max()
    x_grid = np.linspace(x_min - (x_max-x_min)*0.1, x_max + (x_max-x_min)*0.1, 600).reshape(-1, 1)

    plt.figure(figsize=(10, 6))
    plt.scatter(X.flatten(), y, marker='x', s=70, label="Data points")

    for r in results:
        deg = r["degree"]
        m = make_pipeline(PolynomialFeatures(degree=deg, include_bias=False), LinearRegression())
        m.fit(X, y)
        y_grid = m.predict(x_grid)
        lw = 2.5 if deg == best_degree else 1.0
        label = f"deg {deg} (CV RMSE={r['cv_rmse']:.2f})"
        plt.plot(x_grid.flatten(), y_grid, linewidth=lw, label=label)

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(f"Polynomial fits (degrees 1..{results[-1]['degree']}) — selected deg {best_degree}")
    plt.legend()
    plt.grid(axis='y', alpha=0.4)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    if show:
        plt.show()
    plt.close()
    print(f"Saved plot to: {out_path}")

def main(argv=None):
    predict_values = [75, 100]
    predict_values = [35, 40, 45, 50]
    parser = argparse.ArgumentParser(description="Polynomial fitting, comparison, and prediction")
    parser.add_argument("--csv", type=str, default=None, help="Optional CSV file with x and y columns")
    parser.add_argument("--xcol", type=str, default="x", help="Column name for x (if CSV provided)")
    parser.add_argument("--ycol", type=str, default="y", help="Column name for y (if CSV provided)")
    parser.add_argument("--max-degree", type=int, default=5, help="Max polynomial degree to try")
    parser.add_argument("--cv-splits", type=int, default=5, help="K for K-Fold CV")
    parser.add_argument("--predict", type=float, nargs="*", default=predict_values, help="x values to predict (default examples)")
    parser.add_argument("--bootstrap-intervals", action="store_true",
                        help="Compute bootstrap prediction intervals for the chosen model (slower)")
    parser.add_argument("--out-plot", type=str, default="polynomial_fits.png", help="Output plot path")
    args = parser.parse_args(argv)

    df = load_data_from_args(args)
    print("\nData:")
    print(df.to_string(index=False))

    X = df["x"].values.reshape(-1, 1)
    y = df["y"].values

    print("\nFitting polynomials (this may take a second)...")
    results = fit_and_evaluate(X, y, max_degree=args.max_degree, cv_splits=args.cv_splits)
    summary = []
    for r in results:
        summary.append((r["degree"], r["train_r2"], r["train_rmse"], r["cv_rmse"]))
    summary_df = pd.DataFrame(summary, columns=["degree", "train_r2", "train_rmse", "cv_rmse"]).set_index("degree")
    print("\nMetrics for each polynomial degree:")
    print(summary_df.to_string(float_format=lambda v: f"{v:.4g}"))

    # choose best by lowest CV RMSE
    best = min(results, key=lambda r: r["cv_rmse"])
    best_degree = best["degree"]
    print(f"\nSelected best degree by CV RMSE: {best_degree}  (CV RMSE={best['cv_rmse']:.4g})")

    # Refit best model on full data and print formula
    best_model = make_pipeline(PolynomialFeatures(degree=best_degree, include_bias=False), LinearRegression())
    best_model.fit(X, y)
    lin = best_model.named_steps['linearregression']
    coefs = lin.coef_
    intercept = lin.intercept_
    print("\nPolynomial formula (approx):")
    print(textwrap.fill(poly_formula_string(intercept, coefs), width=80))

    # Predictions
    x_to_pred = args.predict
    y_preds = predict_from_model(best_model, x_to_pred)
    print("\nPredictions (using selected best degree):")
    for xv, yv in zip(x_to_pred, y_preds):
        print(f"  x = {xv}  ->  y = {yv:.6g}")

    # Optional bootstrap intervals
    if args.bootstrap_intervals:
        print("\nComputing bootstrap prediction intervals (this will be slower)...")
        median, lower, upper = bootstrap_predict_interval(X, y, best_degree, x_to_pred,
                                                         n_bootstrap=2000, alpha=0.05, random_state=0)
        print("Bootstrap 95% prediction intervals:")
        for xv, m, lo, hi in zip(x_to_pred, median, lower, upper):
            print(f"  x = {xv}  ->  median={m:.6g}, 95% PI = [{lo:.6g}, {hi:.6g}]")

    # Save metrics to CSV
    metrics_out = "poly_fit_metrics.csv"
    summary_df.to_csv(metrics_out)
    print(f"\nSaved metrics CSV to: {metrics_out}")

    # Plot and save figure
    plot_results(df, results, best_degree, out_path=args.out_plot, show=True)

    # Print small usage-snippet the user can copy/paste
    example_snippet = f"""
# Example: reuse the selected model to predict new values
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import numpy as np

# build and fit (same as main)
model = make_pipeline(PolynomialFeatures(degree={best_degree}, include_bias=False), LinearRegression())
X = np.array({df['x'].tolist()}).reshape(-1,1)
y = np.array({df['y'].tolist()})
model.fit(X, y)

# predict
x_new = np.array([{', '.join(map(str, x_to_pred))}]).reshape(-1,1)
y_new = model.predict(x_new)
print(y_new)
"""
    print(textwrap.dedent(example_snippet))

if __name__ == "__main__":
    main()
