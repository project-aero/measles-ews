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
Births [label = 'Births', fontsize = 30];

node [shape = circle,
style=filled,
color = '#B2A389',
fontname = Helvetica] // sets as circles
Imports [label = 'Imports', fontsize = 30];

node [shape = circle,
style=filled,
color = '#F2E3B6',
fontname = Helvetica] // sets as circles
Susceptible [label = 'Susceptible', fontsize = 30];

node [shape = circle,
style=filled,
color = '#D55D45',
fontname = Helvetica] // sets as circles
Exposed [label = 'Exposed', fontsize = 30];

node [shape = circle,
style=filled,
color = '#D55D45',
fontname = Helvetica] // sets as circles
Infectious [label = 'Infectious', fontsize = 30];

node [shape = circle,
style=filled,
color = '#BCDAF2',
fontname = Helvetica] // sets as circles
Recovered [label = 'Recovered', fontsize = 30];


node [shape = oval,
style=filled,
color = '#026773',fontcolor=white,
fontname = Helvetica] // sets as circles
Reported [label = 'Reported', fontsize = 30];

# several 'edge' statements
edge [arrowhead=vee, fontname = Helvetica]
Births->Susceptible [label = '&#956;N', fontsize = 35] 
Susceptible->Exposed [label = '&#946;SI/N', fontsize = 35]
Exposed->Infectious [label = '&#957;E', fontsize = 35]
Imports->Infectious [label = '&#968;', fontsize = 35]
Infectious->Recovered [label = '&#947;I', fontsize = 35]

edge [arrowhead=vee,style=dashed]
Infectious->Reported [label = '&#961;&#947;I', fontsize = 35]
# a 'graph' statement
graph [overlap = true, fontsize = 10,mindist=0.4]
}
"

png(filename = "../figures/seir-diagram.png", width = 357, height = 615, units = "px")
grViz(myplot,engine='neato',options= list(start = 6906))
dev.off()
