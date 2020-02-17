# Scripts for analyzing heterochromatin chips<br />

<strong>Separate spike-in and sample reads (sam_jointmapsplit.pl):</strong><br />
Perl script to separate spike-in and sample reads from a sam file where reads are mapped to a combined genome of the two different species. A joint reference can be made by concatenating the two genomes of interest. Mapping should be done with bwa mem. Reads in the sam files are distinguished based on which of the two genomes they map better to. The script requires two .fai files one for each of the two references. (Do not provide .fai file of the concatenated reference). Two files will be generated: file.split1.sam and file.split2.sam. They contain the reads mapping to the first and second .fai files, respectively. <br />

Usage: <br />
perl sam_jointmapsplit.pl genome1.fasta.fai genome2.fasta.fai mapped.sam > mapped.stats <br /> <br />

<strong>Quantile based spike-in normalization (qqenrich.R):</strong><br />
Rscript for quantile based spike-in normalization. Requires a spike-in reference chip, spike-in reference input, actual sample chip, actual sample input, and spike-in chip, and spike-in input. Arguments in the scripts must be changed for the appropriate winows sizes and chromosome (X and autosomes) names and sizes. <br />
Files must have three tab-delimited columns of chromosome name, start position, and coverage. The start position must be in uniform intervals for the entire chromosome/genome. <br />
A pdf of the normalization scheme will be made containing the normalization profile (topleft), how normalization affects the spike-in (topright), and the actual samples (bottom). <br />
Requires the R package zoo. <br />

Usage: <br />
Rscript qqenrich ref.chip ref.input sample.chip sample.input spike.chip spike.input <br /> <br />

<strong>Infer per basepair coverage around a macs2 peak (gcov.macsrange.peak.pl):</strong><br />
Based on macs2 peak calling, determine the enrichment +/- **N** base pairs around each peak. Requires an enrichment profile file that has three tab-delimited columns of chromosome name, position (every basepiar in the genome), and log2(normalized chip/normalized input). To be used with enrich.heatmap.R (see below). <br />
Usage: <br />

perl gcov.macsrange.peak.pl macs2.peak enrichment.file N > output.enrich.range <br /> <br />

<strong>Collapse overlapping peaks from different macs2 peak-calling on replicates (peakcountbywindow.pl):</strong><br />
Determine the number of peaks that are overlapping (by default within 100bp of each other) across separate replicates. MACS2 peakcalling was applied to each individual replicate. <br />
Usage: <br />

perl peakcountbywindow.pl rep1.macs2.peak rep2.macs2.peak rep3.macs2.peak ... > combined.peaks <br /> <br />

<strong>Heatmap of enrichment around peaks (enrich.heatmap.R):</strong><br />
Generates heatmap of the the enrich.range file (From gcov.macsrange.peak.pl). Requires the r color package RColorBrewer. <br />
Usage: <br />

Rscript enrich.heatmap.R file.enrich.range <br /> <br />

