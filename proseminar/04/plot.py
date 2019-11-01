import os
import sys

import re
import pandas as pd

import plotly.offline as ply
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from plotly.colors import DEFAULT_PLOTLY_COLORS

def find_int_in_string(string):
	ints_in_string = re.findall(r'\d+', string)
	if len(ints_in_string) == 0:
		return None
	return int(ints_in_string[0])

def get_least_ranks(df):
	min_col = None
	min_ranks = None
	for col in df.columns:
		if df[col].isnull().values.any():
			continue
		num_ranks = find_int_in_string(col)
		if num_ranks is not None and (min_ranks is None or num_ranks < min_ranks):
			min_ranks = num_ranks
			min_col = col

	return min_col, min_ranks


def plot_data(dirs, filename):
	df = pd.read_csv(os.path.join(DATA_PATH, filename))

	fig = make_subplots(
		rows=2, cols=2, 
		specs=[[{"colspan": 2}, None], [{}, {}]],
		subplot_titles=("runtime", "speedup", "efficiency"),
		horizontal_spacing=0.1,
		vertical_spacing=0.15)

	problem_size_column = "room_size"
	seq_column = "seq_3D"
	if seq_column in df.columns and not df[seq_column].isnull().values.any():
		# only plot runtime of sequential column without speedup/efficiency
		comparison_column = seq_column
		comparison_num_ranks = 1
		speedup_type = "absolute"

		seq_runtime_trace = go.Scatter(
                    x=df[problem_size_column], y=df[seq_column],
                				legendgroup=seq_column, name=seq_column, marker=dict(color=DEFAULT_PLOTLY_COLORS[0]))
		fig.add_trace(seq_runtime_trace, row=1, col=1)
	else:
		comparison_column, comparison_num_ranks = get_least_ranks(df)
		speedup_type = "relative"

		print("sequential data missing or containing holes!\n--> doing relative speedup instead with '%s' as reference"%comparison_column)


	for i, column in enumerate(df.columns):
		if column in [problem_size_column, seq_column]:
			continue

		num_ranks = find_int_in_string(column)
		print(column)

		runtimes = df[column]
		speedups = (df[comparison_column]*comparison_num_ranks) / runtimes
		efficiencies = speedups / num_ranks

		runtime_trace = go.Scatter(
			x=df[problem_size_column], y=runtimes, 
			legendgroup=column, name=column, marker=dict(color=DEFAULT_PLOTLY_COLORS[1+i]))

		speedup_trace = go.Scatter(
			x=df[problem_size_column], y=speedups, 
			legendgroup=column, marker=dict(color=DEFAULT_PLOTLY_COLORS[1+i]), showlegend=False)

		efficiency_trace = go.Scatter(
			x=df[problem_size_column], y=efficiencies, 
			legendgroup=column, marker=dict(color=DEFAULT_PLOTLY_COLORS[1+i]), showlegend=False)
		
		fig.add_trace(runtime_trace, row=1, col=1)
		fig.add_trace(speedup_trace, row=2, col=1)
		fig.add_trace(efficiency_trace, row=2, col=2)


	fig.update_layout(title_text=filename.split(".")[0], xaxis_type="log")

	fig.update_xaxes(title="room size in all dimensions", type="log")

	fig.update_yaxes(title="runtime in s", rangemode="tozero", type="log", row=1, col=1)
	fig.update_yaxes(title="%s speedup"%speedup_type, rangemode="tozero", row=2, col=1)
	fig.update_yaxes(title="%s efficiency"%speedup_type, range=[0., 1.], row = 2, col = 2)

	ply.plot(fig, filename=os.path.join(PLOTS_PATH, filename.split(".")[0]+".html"))


if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("wrong number of arguments: "+str(sys.argv[1:]))
		sys.exit(-1)

	DATA_PATH = sys.argv[1]
	PLOTS_PATH = sys.argv[2]

	for path, _, files in os.walk(DATA_PATH):
		for filename in files:
			if filename.split(".")[-1] != "csv":
				print("incompatible file: %s"%filename)
				continue
			print("plotting %s"%filename)  
			plot_data(path, filename)
