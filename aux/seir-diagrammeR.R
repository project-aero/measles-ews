library(DiagrammeR)


# The neato engine uses random seeds to do diagram orientation,
# I have set the seed to start = 6906 when I call the function
# by changing the seed you will get a different layout
myplot <- "
digraph neato {

# several 'node' statements
node [shape = circle,
style=filled,
color = '#B2A389',
fontname = Helvetica] // sets as circles
Births;

node [shape = circle,
style=filled,
color = '#B2A389',
fontname = Helvetica] // sets as circles
Imports;

node [shape = circle,
style=filled,
color = '#F2E3B6',
fontname = Helvetica] // sets as circles
Susceptible;

node [shape = circle,
style=filled,
color = '#D55D45',
fontname = Helvetica] // sets as circles
Exposed;

node [shape = circle,
style=filled,
color = '#D55D45',
fontname = Helvetica] // sets as circles
Infectious;

node [shape = circle,
style=filled,
color = '#BCDAF2',
fontname = Helvetica] // sets as circles
Recovered;


node [shape = oval,
style=filled,
color = '#026773',fontcolor=white,
fontname = Helvetica] // sets as circles
Reported;

# several 'edge' statements
edge [arrowhead=vee]
Births->Susceptible [label = '&#956;N'] 
Susceptible->Exposed [label = '&#946;SI/N']
Exposed->Infectious [label = '&#957;E']
Imports->Infectious [label = '&#968;']
Infectious->Recovered [label = '&#947;I']

edge [arrowhead=vee,style=dashed]
Infectious->Reported [label = '&#961;&#947;I']
# a 'graph' statement
graph [overlap = true, fontsize = 10,mindist=0.4]
}
"
grViz(myplot,engine='neato',options= list(start = 6906))
