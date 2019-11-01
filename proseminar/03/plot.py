import os
import sys

import re
import pandas as pd

import plotly.offline as ply
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from plotly.colors import DEFAULT_PLOTLY_COLORS

def find_int_in_string(string):
	return int(re.findall(r'\d+', string)[0])


def plot_data(dirs, filename):
	df = pd.read_csv(os.path.join(DATA_PATH, filename))

	fig = make_subplots(
		rows=2, cols=2, 
		specs=[[{"colspan": 2}, None], [{}, {}]],
		subplot_titles=("runtime", "speedup", "efficiency"),
		horizontal_spacing=0.1,
		vertical_spacing=0.15)

	fig.update_layout(title_text=filename.split(".")[0] + " heat stencil")
	fig.update_xaxes(title="room size in all dimensions", type="log")
	fig.update_yaxes(title="runtime in s", rangemode="tozero", type="log", row=1, col=1)
	fig.update_yaxes(title="absolute speedup", rangemode="tozero", row=2, col=1)
	fig.update_yaxes(title="absolute efficiency", range=[0., 1.], row = 2, col = 2)

	fig.update_layout(xaxis_type="log")

	column = "seq"
	seq_runtime_trace = go.Scatter(
            x=df["room_size"], y=df[column],
            legendgroup=column, name=column, marker=dict(color=DEFAULT_PLOTLY_COLORS[0]))

	fig.add_trace(seq_runtime_trace, row=1, col=1)

	for i, column in enumerate(df.columns):
		if column in ["room_size", "seq"]:
			continue
		num_ranks = find_int_in_string(column)

		runtimes = df[column]
		speedups = df["seq"]/runtimes
		efficiencies = speedups/num_ranks

		runtime_trace = go.Scatter(
			x=df["room_size"], y=runtimes, 
			legendgroup=column, name=column, marker=dict(color=DEFAULT_PLOTLY_COLORS[1+i]))

		speedup_trace = go.Scatter(
			x=df["room_size"], y=speedups, 
			legendgroup=column, marker=dict(color=DEFAULT_PLOTLY_COLORS[1+i]), showlegend=False)

		efficiency_trace = go.Scatter(
			x=df["room_size"], y=efficiencies, 
			legendgroup=column, marker=dict(color=DEFAULT_PLOTLY_COLORS[1+i]), showlegend=False)
		
		fig.add_trace(runtime_trace, row=1, col=1)
		fig.add_trace(speedup_trace, row=2, col=1)
		fig.add_trace(efficiency_trace, row=2, col=2)

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
