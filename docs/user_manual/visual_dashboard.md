# Visual Dashboard

The Visual Dashboard is the first explicit illustration layer of the guided ppigFinder workflow.

It shows:

- workflow progression as connected nodes;
- how input data, operations and outputs relate;
- ORF discovery summaries after ORF prediction;
- compact ORF map preview;
- ORF length distribution;
- downstream decision points for annotation, neighbourhood, AlphaFold/PPI and reports.

This dashboard is implemented with Qt/QPainter to avoid additional browser dependencies. Later versions can replace or extend these visuals with HTML/SVG/D3-style panels and lovis4u-like genomic neighbourhood illustrations.


## Interactive visual dashboard

The visual dashboard now supports interactive navigation and figure expansion.

Available interactions:

- click workflow nodes to navigate to Data, Genome, ORFs, Annotation, AlphaFold/PPI, DaVinci/HPC or Reports;
- open the workflow graph in a dedicated full-screen view;
- open the data/process/output diagram in a dedicated full-screen view;
- open ORF map and ORF length distribution figures in a dedicated full-screen view;
- click ORF figure areas inside the dashboard to expand them.

This prepares the interface for future lovis4u-like genomic neighbourhood diagrams and interactive candidate-selection panels.
