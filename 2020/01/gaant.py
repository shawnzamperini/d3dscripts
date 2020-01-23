import plotly.figure_factory as ff

df = [dict(Task='Conferences', Start='2020-06-01', Finish='2020-06-30')]

fig = ff.create_gantt(df)
