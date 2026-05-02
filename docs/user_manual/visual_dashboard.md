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
