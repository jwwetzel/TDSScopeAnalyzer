# TDSScopeAnalyzer
TDS Scope Waveform Analyzer

First, use the WFMconverter to make the txt files (.fff) out of the .wfm files.

Second, run the TDSwfmAnalyzer (after making) with:

./TDSwfmAnalyzer /path/to/textdatafiles/*.fff

And then it will ask if you want to tuple them.  Enter 1 to say yes.

Then you can run the analyzer part on the ntuples.

You can also type 6 to open a TBrowser to quickly view your results.
