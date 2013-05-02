All files required to run SourceTracker in R are in sourcetracker-X.Y.Z.tar.gz. Data must be integral (counts, not relative abundances).

In the top-level folder, example usage of SourceTracker software is in 'example.r'. The example allows you to run the analysis described in the paper using the OTU table and mapping file contained in the "data" folder. Documentation for all package functions is included in the source code file 'src/SourceTracker.r'.

There is also a convenience wrapper, "sourcetracker_for_qiime.r", for QIIME users that allows you to run SourceTracker from the command line, rather than in R. This requires that you have the path to your SourceTracker.r script stored in the environment variable, "SOURCETRACKER_PATH". One way to do this permanently is to add "export SOURCETRACKER_PATH=/path/to/your/SourceTracker.r" to the end of the .bashrc (Linux) or .bash_profile (Mac) file in your home directory. This script requires a mapping file, and either an OTU table or a taxon table. The input files must be in QIIME format:

Mapping file: tab-delimited, first line contains the column headers, first column header is "#SampleID"
OTU table: tab-delimited, first line is a comment starting with "#", second line contains the column headers, first column header is "#OTU ID".
Taxon table: tab-delimited, first line contains the column headers, first column header is "Taxon" (output from summarize_taxa.py).

Example usage, ">" precedes commands:

To see a listing of command-line parameters:
>R --slave --vanilla --args -h < sourcetracker_for_qiime.r

Run sink predictions using QIIME taxon abundance file:
>R --slave --vanilla --args -t taxa.txt -m map.txt < sourcetracker_for_qiime.r

Run leave-one-out source-sample predictions using QIIME taxon abundance file:
>R --slave --vanilla --args -t taxa.txt -m map.txt -s < sourcetracker_for_qiime.r

Run sink predictions using QIIME OTU table:
>R --slave --vanilla --args -i otutable.txt -m map.txt < sourcetracker_for_qiime.r

Run sink predictions using QIIME OTU table with 1000 burnins, 25 random restarts, and rarefaction depth of 100:
>R --slave --vanilla --args -i otutable.txt -m map.txt -b 1000 -n 25 -r 100 < sourcetracker_for_qiime.r


Change log:
Version 0.9.2:
 - Now allows mapping files without "SourceSink" column when using -s (leave-one-out predictions).

Version 0.9.1:
 - Modified default alpha values to 0.001.
 - Added automatic tuning of alpha values. This is slow, but should be performed before publishing the output.

Version 0.9.0:
 - First 'official' beta release.
 - Includes a script 'sourcetracker_for_qiime.r' for running SourceTracker on QIIME-formatted files directly from the command-line.

