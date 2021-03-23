# Create a Flowchart showing participant counts

# Clear variables and set the path
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# install.packages("DiagrammeR")
library("DiagrammeR")

flowchart <- grViz("digraph flowchart {
      # node definitions with substituted label text
      node [fontname = Helvetica, shape = rectangle, height=1, width=7]        
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
      tab5 [label = '@@5']
      tab6 [label = '@@6']

      # edge definitions with the node IDs
      tab1 -> tab2;
      tab2 -> tab3;
      tab3 -> tab4;
      tab3 -> tab5 -> tab6
      }

      [1]: 'Participants Attended UK Biobank Assessment Centre n=502537'
      [2]: '31 participants withdrawn from the\ cohort n=502506'
      [3]: '15209 participants removed due to\ missing genotype data n=487297'
      [4]: 'Non-Hodgkins Lymphoma (NHL) Cases\ n=3394'
      [5]: 'Controls n=483903'
      [6]: '22605 participants removed due to\ developing non-NHL cancer n=461298'
      ")
png(flowchart, filename="../../figures/flowchart.png")
flowchart
