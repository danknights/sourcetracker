All files required to run SourceTracker in R are in sourcetracker-X.Y.Z.tar.gz. Data must be integral (counts, not relative abundances).

In the top-level folder, example usage of SourceTracker software is in 'example.r'. The example allows you to run the analysis described in the paper using the OTU table and mapping file contained in the "data" folder. Documentation for all package functions is included in the source code file 'src/SourceTracker.r'.

There is also a convenience wrapper, "sourcetracker_for_qiime.r", for QIIME users that allows you to run SourceTracker from the command line, rather than in R. This requires that you have the path to your top-level SourceTracker repository folder script stored in the environment variable, "SOURCETRACKER_PATH". One way to do this permanently is to add "export SOURCETRACKER_PATH=/path/to/your/sourcetracker/repository/folder" to the end of the .bash_profile or .Renviron file in your home directory. You can do this with one of these commands:

For most systems:
```bash
echo >> $HOME/.bash_profile
echo "SOURCETRACKER_PATH=$PWD" >> $HOME/.bash_profile
```

If that fails:
```bash
echo >> $HOME/.Renviron
echo "SOURCETRACKER_PATH=$PWD" >> $HOME/.Renviron
```

This script requires a mapping file, and either an OTU table or a taxon table. The input files must be in QIIME format:

Mapping file: tab-delimited, first line contains the column headers, first column header is "#SampleID"
OTU table: tab-delimited, first line is a comment starting with "#", second line contains the column headers, first column header is "#OTU ID".
Taxon table: tab-delimited, first line contains the column headers, first column header is "Taxon" (output from summarize_taxa.py).

Example usage, ">" precedes commands:

To see a listing of command-line parameters:
```bash
Rscript sourcetracker_for_qiime.r -h
```

Run sink predictions using QIIME taxon abundance file:
```bash
Rscript sourcetracker_for_qiime.r -t taxa.txt -m map.txt
```

Run leave-one-out source-sample predictions using QIIME taxon abundance file:
```bash
R sourcetracker_for_qiime.r -t taxa.txt -m map.txt -s
```

Run sink predictions using QIIME OTU table:
```bash
R sourcetracker_for_qiime.r -i otutable.txt -m map.txt
```

Run sink predictions using QIIME OTU table with 1000 burnins, 25 random restarts, and rarefaction depth of 100:
```bash
R sourcetracker_for_qiime.r -i otutable.txt -m map.txt -b 1000 -n 25 -r 100
```

Change log:

Version 1.0:
 - Fixed bug in leave-one-out predictions when suppressing full results; ported to github.

Version 0.9.8:
 - Fixed bug in leave-one-out predictions, removed extraneous output

Version 0.9.7:
 - Now outputs per-taxon or per-OTU source estimates for every sample

Version 0.9.6:
 - Fixed bug when user provides only 1 source or 1 sink, and R's default behavior drops the empty dimension of the source/sink OTU matrix thus converting it to a vector and breaking subsequent matrix operations (http://radfordneal.wordpress.com/2008/08/20/design-flaws-in-r-2-%E2%80%94-dropped-dimensions/).

Version 0.9.5:
 - Now expects SOURCETRACKER_PATH to point to the SourceTracker parent directory. This will allow sys admins to install SourceTracker in an arbitrary location, while the user does not need to know the location. The user can then run with: "Rscript $SOURCETRACKER_PATH/sourcetracker_for_qiime.r". 

Version 0.9.4:
 - Output 'map.txt' now has TRUE/FALSE columns for contamination at thresholds 0.05, 0.10, ..., 0.95.

Version 0.9.3:
 - Improved memory management
 - No fails gracefully when a sample contains zero observations

Version 0.9.2:
 - Now allows mapping files without "SourceSink" column when using -s (leave-one-out predictions).

Version 0.9.1:
 - Modified default alpha values to 0.001.
 - Added automatic tuning of alpha values. This is slow, but should be performed before publishing the output.

Version 0.9.0:
 - First 'official' beta release.
 - Includes a script 'sourcetracker_for_qiime.r' for running SourceTracker on QIIME-formatted files directly from the command-line.

