# Introduction
This code performs the calculations for the publication <b>Most Undirected Random Graphs are Amplifiers of Selection for Birth-death Dynamics, but Suppressors of Selection for death-Birth Dynamics</b> by <i>Laura Hindersin</i> and <i>Arne Traulsen</i> (Journal link coming soon).

An IPython notebook is available <a href="http://nbviewer.ipython.org/github/hindersin/randomGraphs/blob/master/randomGraphs.ipynb">online here</a>


#Included files
<table>
<tr><td><i>randomGraphs.ipynb</i>: an interactice IPython notebook that performs the calculation and visualization process.</td></tr>

<tr><td><i>randomGraphsNumerical.py</i>: the main Python file that generates the random graphs and calculates the fixation probability.</td></tr>

<tr><td><i>transitionMatrix.py</i>: for a given undirected graph, this file generates the transition matrix for the Moran process on that graph.</td></tr>

<tr><td><i>transitionMatrixDirected.py</i>: same as above for directed graphs.</td></tr>

<tr><td><i>AmplifierQ.py</i>: determines according to the fixation probability, whether a graph is an amplifier or suppressor of selection.</td></tr>

<tr><td><i>Plot.py</i>: visualization of the results.</td></tr>

<tr><td><i>LICENSE</i>: MIT License</td></tr>

<tr><td><i>README.md</i>: this file</td></tr>
</table>

#Dependencies
<b>randomGraphs.ipynb</b> was created using <i>IPython</i> version 3.1.0

The python files were created <i>python</i> 2.7.9, <i>numpy</i> version 1.9.2, <i>networkx</i> version 1.9.1, <i>sympy</i> version 0.7.6, <i>matplotlib</i> version 1.4.3

#License
See <i>LICENSE</i> for details.
