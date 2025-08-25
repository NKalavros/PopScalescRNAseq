"""Streamlit dashboard (and optional CLI exporter) to explore scIB benchmark CSV results across tissues/studies.

Run with:
	streamlit run viz/tabulate_scib.py

It auto-discovers files matching pattern:
	<TissueAtlas>/<Study>/figures/scib_results_<Study>.csv

Expected CSV format (flexible):
	Rows: methods (embedding / integration method names)
	Columns: metrics (ARI, ASW, NMI, etc.)
	Optionally one identifying column like 'method' / 'Method' / 'embedding'.

Features (Streamlit UI):
	- Filter by Tissue, Study, Method, Metric
	- Interactive wide & long tables
	- Aggregate (mean / std) over selected studies
	- Heatmap (Methods vs Studies) for a selected metric
	- Per-metric bar charts ranking methods (mean across studies)
	- Download filtered data
	- Basic theming & highlighting of best scores

CLI Export (headless) usage:

	python tabulate_scib.py --root /path/to/results --export out_dir

Produces:
	out_dir/combined_long.csv          (all discovered rows, long format)
	out_dir/aggregate_mean.csv         (Method x Metric table of means)
	out_dir/aggregate_mean_std.csv     (Method x Metric with mean±std string)
	out_dir/heatmap_<metric>.html      (interactive Plotly heatmap for each metric)
	out_dir/ranking_<metric>.html      (interactive Plotly bar ranking per metric)

You can then download these files from HPC (e.g. scp / rsync) without running Streamlit.

Remote usage (HPC -> local) for the dashboard:

	1. On HPC: export SCIB_ATLASES=BreastAtlas,KidneyAtlas && cd /gpfs/scratch/nk4167/PopScalescRNAseq && streamlit run tabulate_scib.py --server.address 0.0.0.0 --server.port 12345 -- --root /gpfs/scratch/nk4167/
	2. SSH tunnel locally:
		 ssh -L 12345:localhost:12345 bigpurple
		 a. ssh to your compute node if needed
		 ssh -L 12345:localhost:12345 nk4167@a100-4005
	3. Open http://localhost:12345 in your browser.

Environment variables (override defaults):
	SCIB_ROOT : default root search path
	SCIB_FILE_GLOB : override file pattern (default '*Atlas/*/figures/scib_results_*.csv')

Robustness:
	- Attempts to auto-detect method name column
	- Converts values to numeric where possible
	- Ignores empty / non-numeric metric columns
"""

from __future__ import annotations

import os
import io
import glob
import argparse
from typing import List, Optional

import pandas as pd
import numpy as np

try:
	import streamlit as st  # type: ignore
	import plotly.express as px  # type: ignore
except Exception as e:  # pragma: no cover - graceful degradation
	raise SystemExit(
		"This script requires streamlit and plotly. Install via: pip install streamlit plotly"
	) from e


# ------------------------------ Configuration ------------------------------ #
DEFAULT_ROOT = os.environ.get("SCIB_ROOT", ".")  # root to start searching from; can be changed in UI / env
FILE_GLOB = os.environ.get("SCIB_FILE_GLOB", "*Atlas/*/figures/scib_results_*.csv")  # still used for CLI pattern override if desired
DEFAULT_ATLASES = [a.strip() for a in os.environ.get("SCIB_ATLASES", "BreastAtlas,KidneyAtlas").split(",") if a.strip()]
METHOD_ID_CANDIDATES = [
	"method",
	"Method",
	"embedding",
	"Embedding",
	"integration",
	"Integration",
]


@st.cache_data(show_spinner=False)
def discover_result_files(root: str, atlases: Optional[List[str]] = None) -> List[str]:
	"""Discover scIB result CSVs quickly by scanning specified atlas folders.

	Instead of a broad recursive glob, we enumerate studies under each atlas
	folder (<Atlas>/<Study>/figures/scib_results_<Study>.csv) for speed on HPC.

	Parameters
	----------
	root : str
		Root directory containing atlas folders (e.g. /gpfs/scratch/nk4167)
	atlases : list[str] | None
		Names of atlas directories to scan. Defaults to DEFAULT_ATLASES.

	Returns
	-------
	list[str]
		Sorted list of discovered CSV paths.
	"""
	atlases = atlases or DEFAULT_ATLASES
	results: List[str] = []
	for atlas in atlases:
		atlas_path = os.path.join(root, atlas)
		if not os.path.isdir(atlas_path):
			continue
		# Each immediate subdirectory is assumed to be a study
		try:
			studies = [d for d in os.listdir(atlas_path) if os.path.isdir(os.path.join(atlas_path, d))]
		except OSError:
			continue
		for study in studies:
			fig_dir = os.path.join(atlas_path, study, "figures")
			if not os.path.isdir(fig_dir):
				continue
			# Expected file name starts with scib_results_<Study>.csv but we allow any scib_results_*.csv inside figures.
			for fname in os.listdir(fig_dir):
				if not fname.endswith('.csv'):
					continue
				if not fname.startswith('scib_results_'):
					continue
				results.append(os.path.join(fig_dir, fname))
	return sorted(set(results))


def _find_method_col(df: pd.DataFrame) -> Optional[str]:
	for c in METHOD_ID_CANDIDATES:
		if c in df.columns:
			return c
	# If first column is object & not purely numeric, treat as method id
	first = df.columns[0]
	if df[first].dtype == object:
		return first
	return None


def _clean_method_names(series: pd.Series) -> pd.Series:
	return (
		series.astype(str)
		.str.strip()
		.str.replace("_embedding", "", regex=False)
		.str.replace(".h5ad", "", regex=False)
		.str.replace("x_", "", regex=False)
		.str.replace("X_", "", regex=False)
	)


@st.cache_data(show_spinner=False)
def load_all_results(files: List[str]) -> pd.DataFrame:
	records = []
	for path in files:
		try:
			df = pd.read_csv(path)
		except Exception as e:  # pragma: no cover
			st.warning(f"Failed to read {path}: {e}")
			continue
		method_col = _find_method_col(df)
		if method_col is None:
			# Assume index are methods
			df = df.copy()
			df.index.name = "method"
			df.reset_index(inplace=True)
			method_col = "method"
		df = df.copy()
		df[method_col] = _clean_method_names(df[method_col])

		# Remove completely empty columns & non-metric placeholders
		metric_cols = [
			c for c in df.columns if c != method_col and df[c].notna().any()
		]
		# Derive tissue & study from path: <root>/<TissueAtlas>/<Study>/figures/...
		parts = path.split(os.sep)
		# Find the index of 'figures'
		try:
			fig_idx = parts.index("figures")
			study = parts[fig_idx - 1]
			tissue = parts[fig_idx - 2]
		except Exception:
			study = os.path.basename(path).replace("scib_results_", "").replace(
				".csv", ""
			)
			tissue = "Unknown"

		# Melt to long form
		long_df = df[[method_col] + metric_cols].melt(
			id_vars=method_col, var_name="Metric", value_name="Value"
		)
		long_df["Tissue"] = tissue
		long_df["Study"] = study
		long_df.rename(columns={method_col: "Method"}, inplace=True)

		# Coerce to numeric silently
		long_df["Value"] = pd.to_numeric(long_df["Value"], errors="coerce")
		long_df = long_df.dropna(subset=["Value"])  # keep only numeric rows
		records.append(long_df)
	if not records:
		return pd.DataFrame(columns=["Tissue", "Study", "Method", "Metric", "Value"])
	all_long = pd.concat(records, ignore_index=True)
	# Standardize tissue names (remove 'Atlas' suffix for cleaner filters)
	all_long["Tissue"] = all_long["Tissue"].str.replace("Atlas", "", regex=False)
	return all_long


def _style_dataframe(df: pd.DataFrame) -> "pd.io.formats.style.Styler":  # type: ignore
	# Highlight maxima per column (metrics) if numeric
	numeric_cols = df.select_dtypes(include=[np.number]).columns
	def highlight_max(s):  # noqa: D401
		is_max = s == s.max()
		return ["background-color: #264653; color: #fff;" if v else "" for v in is_max]
	styler = df.style.format(precision=4)
	for c in numeric_cols:
		# Use lowercase 'viridis' (matplotlib colormap names are lowercase) and guard against errors
		try:
			styler = styler.background_gradient(cmap="viridis", subset=[c])
		except ValueError:
			# Fallback to 'Greys' if somehow unavailable
			styler = styler.background_gradient(cmap="Greys", subset=[c])
	styler = styler.apply(highlight_max, subset=numeric_cols)
	return styler


def main():  # pragma: no cover - Streamlit runtime function
	st.set_page_config(
		page_title="scIB Benchmark Explorer",
		page_icon="🧬",
		layout="wide",
	)
	st.markdown(
		"""
		<style>
		.block-container {padding-top: 1rem;}
		.metric-badge {background:#264653;color:#fff;padding:2px 6px;border-radius:4px;font-size:0.75rem;}
		</style>
		""",
		unsafe_allow_html=True,
	)
	st.title("🧬 scIB Benchmark Explorer")
	st.caption(
		"Interactive dashboard to compare integration / embedding methods across studies and tissues."
	)

	# ---------------------------- Sidebar Controls ---------------------------- #
	with st.sidebar:
		st.header("Filters & Data")
		root = st.text_input("Search root", DEFAULT_ROOT)
		atlas_input = st.text_input("Atlas folders (comma-separated)", ",".join(DEFAULT_ATLASES))
		atlas_list = [a.strip() for a in atlas_input.split(',') if a.strip()]
		files = discover_result_files(root, atlas_list)
		st.write(f"Found {len(files)} result file(s).")
		if st.toggle("Show file list", value=False):
			st.write(files)

	data = load_all_results(files)
	if data.empty:
		st.error("No valid scIB result CSVs found. Adjust the root path.")
		st.stop()

	tissues = sorted(data["Tissue"].unique())
	sel_tissues = st.sidebar.multiselect("Tissues", tissues, default=tissues)
	d1 = data[data["Tissue"].isin(sel_tissues)]

	studies = sorted(d1["Study"].unique())
	sel_studies = st.sidebar.multiselect("Studies", studies, default=studies)
	d2 = d1[d1["Study"].isin(sel_studies)]

	methods = sorted(d2["Method"].unique())
	default_methods = methods  # could filter top-n later
	sel_methods = st.sidebar.multiselect("Methods", methods, default=default_methods)
	d3 = d2[d2["Method"].isin(sel_methods)]

	metrics = sorted(d3["Metric"].unique())
	sel_metrics = st.sidebar.multiselect("Metrics", metrics, default=metrics)
	df_filt = d3[d3["Metric"].isin(sel_metrics)]

	if df_filt.empty:
		st.warning("No data after applying filters.")
		st.stop()

	st.subheader("Filtered Long Table")
	with st.expander("Show / Hide Long Format Data", expanded=False):
		st.dataframe(df_filt, use_container_width=True, height=400)
		csv_buf = io.StringIO()
		df_filt.to_csv(csv_buf, index=False)
		st.download_button(
			"Download filtered (long) CSV",
			data=csv_buf.getvalue(),
			file_name="scib_filtered_long.csv",
			mime="text/csv",
		)

	# ---------------------------- Aggregated Table ---------------------------- #
	st.subheader("Aggregate (Mean ± SD across Studies)")
	agg = (
		df_filt.groupby(["Method", "Metric"], as_index=False)["Value"]
		.agg(["mean", "std", "count"])
		.reset_index()
	)
	agg["mean"] = agg["mean"].astype(float)
	agg_pivot = agg.pivot(index="Method", columns="Metric", values="mean")
	st.dataframe(agg_pivot, use_container_width=True)

	# Styled wide table (optional) if not too large
	if agg_pivot.shape[1] <= 20:
		st.markdown("**Styled Table (highlight best per metric)**")
		st.dataframe(_style_dataframe(agg_pivot), use_container_width=True)

	# ---------------------------- Heatmap ---------------------------- #
	st.subheader("Heatmap: Methods vs Studies")
	metric_for_heatmap: Optional[str] = None
	if len(sel_metrics) == 1:
		metric_for_heatmap = sel_metrics[0]
	else:
		metric_for_heatmap = st.selectbox(
			"Select metric for heatmap", sel_metrics, index=0
		)
	heat_df = (
		df_filt[df_filt["Metric"] == metric_for_heatmap]
		.groupby(["Study", "Method"], as_index=False)["Value"].mean()
	)
	pivot_heat = heat_df.pivot(index="Method", columns="Study", values="Value")
	# Ensure order stable
	pivot_heat = pivot_heat.reindex(index=sorted(pivot_heat.index))
	fig_heat = px.imshow(
		pivot_heat,
		color_continuous_scale="Viridis",
		aspect="auto",
		origin="lower",
		labels=dict(color=metric_for_heatmap),
	)
	fig_heat.update_layout(height=500, margin=dict(l=10, r=10, t=30, b=10))
	st.plotly_chart(fig_heat, use_container_width=True)

	# ---------------------------- Ranking Bar Charts ---------------------------- #
	st.subheader("Method Rankings (Higher = Better)")
	col_rank_metric, col_rank_methods = st.columns([2, 3])
	with col_rank_metric:
		metric_for_rank = st.selectbox(
			"Metric for ranking", sel_metrics, index=0, key="rank_metric"
		)
	rank_df = (
		df_filt[df_filt["Metric"] == metric_for_rank]
		.groupby("Method", as_index=False)["Value"].mean()
		.sort_values("Value", ascending=False)
	)
	fig_rank = px.bar(
		rank_df,
		x="Value",
		y="Method",
		orientation="h",
		color="Value",
		color_continuous_scale="Viridis",
		range_color=(rank_df.Value.min(), rank_df.Value.max()),
	)
	fig_rank.update_layout(
		yaxis=dict(categoryorder="total ascending"),
		height=600,
		margin=dict(l=10, r=10, t=30, b=10),
	)
	col_rank_methods.plotly_chart(fig_rank, use_container_width=True)

	# ---------------------------- Per-Study Drilldown ---------------------------- #
	st.subheader("Per-Study Detailed Table")
	selected_study = st.selectbox("Study", sel_studies, index=0)
	ds = df_filt[df_filt["Study"] == selected_study]
	study_pivot = ds.pivot_table(
		index="Method", columns="Metric", values="Value", aggfunc="mean"
	)
	st.dataframe(study_pivot, use_container_width=True)

	# ---------------------------- Footer ---------------------------- #
	with st.expander("About / Help", expanded=False):
		st.markdown(
			"""
			- Scores are aggregated as simple arithmetic means across selected studies.
			- Best method highlighting is per-metric column (highest value).
			- Adjust root path in sidebar if directory layout differs.
			- Download button provides the filtered long-format dataset for external analysis.
			- Extend this script with statistical tests or additional visualizations as needed.
			"""
		)


def run_cli_export(root: str, out_dir: str, pattern: Optional[str] = None, atlases: Optional[List[str]] = None) -> None:
	"""Headless export aggregations (no Streamlit server)."""
	pattern = pattern or FILE_GLOB  # kept for backward compat; not used unless fallback triggered
	files = discover_result_files(root, atlases)
	if not files:
		# Fallback to legacy broad glob only if pattern explicitly provided
		if pattern:
			legacy_pattern = os.path.join(root, pattern)
			files = sorted(glob.glob(legacy_pattern))
		if not files:
			print("No files found. Nothing to export.")
		return
	data = load_all_results(files)
	if data.empty:
		print("No valid data parsed. Aborting.")
		return
	os.makedirs(out_dir, exist_ok=True)
	# Save long format
	long_path = os.path.join(out_dir, "combined_long.csv")
	data.to_csv(long_path, index=False)
	# Aggregate
	agg = (
		data.groupby(["Method", "Metric"], as_index=False)["Value"].agg(["mean", "std", "count"]).reset_index()
	)
	agg_mean = agg.pivot(index="Method", columns="Metric", values="mean")
	agg_mean.to_csv(os.path.join(out_dir, "aggregate_mean.csv"))
	# Mean±Std string table
	def fmt_row(row):
		return f"{row['mean']:.4f}±{row['std']:.4f}" if not pd.isna(row['std']) else f"{row['mean']:.4f}"
	combo = agg.copy()
	combo["mean_std"] = combo.apply(fmt_row, axis=1)
	combo_pivot = combo.pivot(index="Method", columns="Metric", values="mean_std")
	combo_pivot.to_csv(os.path.join(out_dir, "aggregate_mean_std.csv"))
	print(f"Wrote long + aggregate tables to {out_dir}")
	# Generate simple Plotly HTML artifacts per metric
	metrics = sorted(data.Metric.unique())
	for metric in metrics:
		metric_df = data[data.Metric == metric]
		# Heatmap: Methods (rows) x Studies (columns)
		hm = metric_df.groupby(["Study", "Method"], as_index=False)["Value"].mean()
		pivot_hm = hm.pivot(index="Method", columns="Study", values="Value")
		fig_heat = px.imshow(pivot_hm, color_continuous_scale="Viridis", aspect="auto", origin="lower", labels=dict(color=metric))
		fig_heat.write_html(os.path.join(out_dir, f"heatmap_{metric}.html"), include_plotlyjs="cdn")
		# Ranking bar chart
		rank = metric_df.groupby("Method", as_index=False)["Value"].mean().sort_values("Value", ascending=False)
		fig_rank = px.bar(rank, x="Value", y="Method", orientation="h", color="Value", color_continuous_scale="Viridis")
		fig_rank.update_layout(yaxis=dict(categoryorder="total ascending"))
		fig_rank.write_html(os.path.join(out_dir, f"ranking_{metric}.html"), include_plotlyjs="cdn")
	print(f"Exported per-metric heatmaps & rankings for {len(metrics)} metrics.")


def parse_args():  # pragma: no cover
	parser = argparse.ArgumentParser(description="scIB results dashboard / exporter")
	parser.add_argument("--root", default=DEFAULT_ROOT, help="Root directory to search (default env SCIB_ROOT or '.')")
	parser.add_argument("--export", metavar="OUT_DIR", help="Run in headless export mode to OUT_DIR (no Streamlit)")
	parser.add_argument("--pattern", default=None, help="Override file glob pattern")
	return parser.parse_args()


if __name__ == "__main__":  # pragma: no cover
	args = parse_args()
	if args.export:
		run_cli_export(args.root, args.export, args.pattern)
	else:
		# In Streamlit context extra CLI args after -- are ignored by streamlit itself
		main()
