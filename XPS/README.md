# XPS Data 
This branch contains processed X-ray Photoelectron Spectroscopy (XPS) data used to assess polypropylene (PP) surface oxidation and microbial interactions in the study:

“Cold-water gut isolate from threespine stickleback (Gasterosteus aculeatus) reveals polypropylene surface oxidation and co-culture inhibition.”

All measurements were performed in triplicate using a Thermo Fisher Scientific K-Alpha XPS with a monochromatic Al Kα source. Data were acquired and deconvoluted using CasaXPS software and subsequently processed in Excel, R, and GraphPad Prism. Dates denote biological replicates.

XPS/ - Folder containing raw and processed XPS spectra for each sample (e.g., control, monoculture, and coculture). Each file includes individual high-resolution C 1s and O 1s scans exported from CasaXPS software. All individual samples were scanned in triplicate.

XPS_Chemical_Composition.csv - Summary table of elemental composition and relative atomic percentages derived from survey spectra for each treatment and replicate.

xps_permanova.r - R script performing pairwise PERMANOVA on peak-area percentages or elemental composition data using the vegan package.

Experimental Details • Instrument: Thermo Fisher Scientific K-Alpha, Al Kα source (1486.6 eV) • Spot size: 400 µm • Energy step size (high-res scans): 0.1 eV • Pass energy: 20 eV • Charge compensation: Flood gun enabled • Replicates: Three biological replicates per condition; three technical scans each • Sample types analysed: 1. Uninoculated PP control 2. KMM 191 monoculture 3. KMM 195 monoculture 4. KMM 191 + KMM 195 coculture
