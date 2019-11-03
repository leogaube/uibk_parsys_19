import os
import sys

from math import log2
import re
import pandas as pd

import plotly.offline as ply
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import colorlover as cl

def fade_color(color, strength):
	values = [i*strength for i in find_ints_in_string(color)]
	return "rgb(%d, %d, %d)" % (values[0], values[1], values[2])

def find_ints_in_string(string):
	ints_in_str = re.findall(r'\d+', string)
	if len(ints_in_str) == 0:
		return None
	return [int(i) for i in ints_in_str]

def find_int_in_string(string):
	ints = find_ints_in_string(string)
	if ints is None:
		return None
	return ints[0]

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

def get_most_ranks(df):
	max_col = None
	max_ranks = None
	for col in df.columns:
		if df[col].isnull().values.any():
			continue
		num_ranks = find_int_in_string(col)
		if num_ranks is not None and (max_ranks is None or num_ranks > max_ranks):
			max_ranks = num_ranks
			max_col = col

	return max_col, max_ranks


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

		seq_color = "rgb(0,0,205)"
		seq_runtime_trace = go.Scatter(
                    x=df[problem_size_column], y=df[seq_column],
					legendgroup=seq_column, name=seq_column, marker=dict(color=seq_color))
		fig.add_trace(seq_runtime_trace, row=1, col=1)
	else:
		comparison_column, comparison_num_ranks = get_least_ranks(df)
		speedup_type = "relative"

		print("sequential data missing or containing holes!\n--> doing relative speedup instead with '%s' as reference"%comparison_column)


	_, max_ranks = get_most_ranks(df)
	next_color_index = 0

	COLORS = cl.scales["9"]["seq"]
	COLOR_NAMES = ["Reds", "Greens", "Purples", "Greys", "Oranges"]
	colors = {}

	mpi_columns = [column for column in df.columns if column not in [problem_size_column, seq_column]]
	for i, column in enumerate(sorted(mpi_columns, key=find_int_in_string)):
		program_group, rank_group = column.split("_fillup_")
		num_ranks = find_int_in_string(column)
		if program_group not in colors:
			colors[program_group] = COLORS[COLOR_NAMES[next_color_index]]
			next_color_index += 1
		color = colors[program_group][int(2+log2(num_ranks))]

		runtimes = df[column]
		speedups = (df[comparison_column]*comparison_num_ranks) / runtimes
		efficiencies = speedups / num_ranks

		runtime_trace = go.Scatter(
			x=df[problem_size_column], y=runtimes, 
			legendgroup=num_ranks, name=column, marker=dict(color=color), visible=True if (num_ranks == max_ranks) else "legendonly")

		speedup_trace = go.Scatter(
			x=df[problem_size_column], y=speedups, 
			legendgroup=num_ranks, marker=dict(color=color), showlegend=False, visible=True if (num_ranks == max_ranks) else "legendonly")

		efficiency_trace = go.Scatter(
			x=df[problem_size_column], y=efficiencies, 
			legendgroup=num_ranks, marker=dict(color=color), showlegend=False, visible=True if (num_ranks == max_ranks) else "legendonly")
		
		fig.add_trace(runtime_trace, row=1, col=1)
		fig.add_trace(speedup_trace, row=2, col=1)
		fig.add_trace(efficiency_trace, row=2, col=2)


	fig.update_layout(title_text=filename.split(".")[0], xaxis_type="log")

	fig.update_xaxes(title="room size in all dimensions", type="log", tickvals=df[problem_size_column])

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
