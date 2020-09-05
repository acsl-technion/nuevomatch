# Information

The official source code for NuevoMatch as described in "A Computational Approach to Packet Classification" [1]. This library contains the source-code for training and running RQ-RMI models (in a stand-alone fashion), the algorithm for generating iSets from ClassBench [2] files, NuevoMatch classifier, and other several tools and benchmarks.

# Prerequisites

General prerequisites:
* A Linux operating system / WSL for Windows operating system
* Intel CPU that supports AVX2 and FMA extensions (the configuration script will check that for you)
* g++ compiler that supports C++14
* GNU Make
* Python 3.5+
* Python3-dev package (must be compatible with Python 3.5+)
* Make sure both hyperthreading and frequency scaling are disabled for maximum performance.

Python prerequisites:
* Tensorflow (tested with version 1.13.1)
* Numpy (tested with version 1.16.2)
* Distutils (tested with version 3.6.7)
* Matplotlib (tested with version 3.0.3)

# Installation

* Run the configuration script. It checks that all prerequisites are met, and generates a Makefile. 
```
./configure
```
* Run GNU make (-f option recommended)
```
make -f
```

A directory named "bin" will be created with executable and scripts.

# Simple Python Examples

The "examples" directory contains three simple Python examples.
* ``simple_rmi.py:`` (1) generates a key-value database according to Zipf law (2) trains a limited RMI [3] model (3) saves the model to a file and exports a chart with some information.
* ``simple_rqrmi.py:`` Does the same as simple_rmi.py but with RQRMI training & theory (responsibilities, dataset sampling, etc.). The saved model can be loaded using *bench_rqrmi.exe* tool for performance measurements.
* ``simple_lookup.py:`` Perform a full key-value lookup using RQRMI and secondary search in Python. Exports the data-structure to a file. The saved data-structure can be loaded using *bench_lookup.exe* tool for performance measurements.

# Tools and Microbenchmarks

All tools and microbenchmarks are available under the "bin" directory. Use "-\-help" for usage instructions.

Tools:
* ``tool_classifier.exe:`` Use this tool to evaluate NuevoMatch against CutSplit [5], NeuroCuts [4], and TupleMerge [6].
* ``tool_trace_generator.exe:`` Generates accurate packet traces (5-tuple + matched priority) from ClassBench files with uniform rule distribution.
* ``tool_locality.exe:`` Locality tool for generating skewed traces. Can be used to extract temporal locality from PCAP files (together with tcpdump), or
to generate Zipf distribution with various parameters.
* ``nuevomatch.py:`` Generates NuevoMatch classifiers from ClassBench files. 
* ``ruleset_analysis.py:`` Analyze ClassBench rulesets. Mainly used for debugging iSets.
* ``pack_neurocuts.py:`` Converts NeuroCuts [4] classifiers to binary files that can be read using our native implementation of NeuroCuts.
* ``classifier_analysis.py:`` Extracts RQRMI models from NuevoMatch classifier files.

Microbenchmarks:
* ``bench_rqrmi.exe``: Tests the performance of an RQRMI model using both serial code and AVX acceleration. The model is loaded from a file (use simple_rqrmi.py for generating such files).
* ``bench_lookup.exe``: Tests both the performance and the correctness of RQRMI using a secondary search. The data-structure is loaded from a file (use simple_lookup.py for generating such files).
* ``bench_echo.exe``: Tests the communication performance between two threads.
* ``bench_reducer.exe``: Tests the communication performance between several threads.

# A Full Example

You can use the tools to run NuevoMatch as follows:

1. Get a ClassBench ruleset. You can generate one using the [original tool](http://www.arl.wustl.edu/classbench/db_generator.tar.gz).
2. Train a model using ``./bin/nuevomatch.py`` (takes about 40 minutes):

```./bin/nuevomatch.py -f my_ruleset.txt -o my_model.data --min-size 64 --max-error 64```

3. Use ``./bin/tool_trace_generator.exe`` to generate a trace from a ClassBench ruleset:

```./bin/tool_trace_generator.exe -f my_ruleset.txt -n 10000000 -o my_trace.txt```

4. You can shuffle the trace according to Zipf temporal locality:

```./bin/tool_locality.exe --trace-file my_trace.txt --output-trace my_trace_modified.txt --zipf --alpha 0.99```

5. Use ``./bin/tool_classifier.exe`` to test NuevoMatch:

```./bin/tool_classifier.exe -m nuevomatch -l -in my_model.data --trace my_trace_modified.txt --parallel 1 --trace-silent --max-subsets 1 --remainder-type tuplemerge --force-remainder-build -v 3```


# License and Credits

MIT license. See LICENSE.md for more details.

____

[1] A computational approach to packet classification (ACM SIGCOMM, 2020)

[2] Classbench: A packet classification benchmark (ACM TON, 2007)

[3] The case for learned index structures (ACM SIGMOD, 2018)

[4] Neural packet classification (ACM SIGCOMM, 2019)

[5] Cutsplit: A decision-tree combining cutting and splitting for scalable packet classification (IEEE INFOCOM, 2018)

[6] TupleMerge: Fast software packet processing for online packet classification (ACM TON, 2019)
