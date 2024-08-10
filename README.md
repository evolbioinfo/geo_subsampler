# GEO subsampler

Geo_subsampler subsamples a given phylogenetic tree to rebalance the samples at different locations 
according to user-specified proportions. Moreover, for each location the kept samples are chosen 
in a balanced way over the sampling intervals (e.g. months).
With these constraints in mind, the script uses phylogenetic diversity [[Faith 1992]](https://www.sciencedirect.com/science/article/pii/0006320792912013) 
to pick the samples to be removed.
Additional options allow to keep all the samples before a certain data, 
and to ensure a minimal number of samples picked by location, despite the other criteria.

### Article

If you find geo_sampler useful, please cite: 

A Zhukova, L Blassel, F Lemoine, M Morel, J Voznica, O Gascuel (2021) __Origin, evolution and global spread of SARS-CoV-2__
CRAS 344(1): 57-75 doi:[10.5802/crbiol.29](https://doi.org/10.5802/crbiol.29).


## Installation
To install geo_subsampler, first install python 3, then run:

```bash
pip3 install geo_subsampler
```



## Input data
As an input, one needs to provide a **NON**-dated phylogenetical tree in [newick](https://en.wikipedia.org/wiki/Newick_format) format,
a metadata table containing tip names, locations and states, 
in tab-delimited (by default) or csv format (to be specified with *'--sep ,'* option).
To subsample according to user-specified proportions, one should also input a location case counts, 
as tab(or comma, see above)-separated table whose first column contains locations and the second case counts.

### Example
The folder [example_data](example_data) contains an example of an input tree ([covid.nwk](example_data/covid.nwk)) 
representing an early SARS-COV-2 epidemic,
the corresponding metadata table ([metadata.tab](example_data/metadata.tab)), and a case count table ([cases.tab](example_data/cases.tab)).

The input tree contains 11 167 sampled tips.


The metadata table is a tab-separated file, containing tip ids in the first column, 
their countries of sampling in the second column, and the sampling dates in the third column:

id	| country	| sampling date
----- |  ----- | -----
EPI_ISL_402119	| China	| 30/12/2019
EPI_ISL_402123	| China	| 24/12/2019
EPI_ISL_403962	| Thailand	| 08/01/2020
... | ... | ...

The case count table contains numbers of declared cases for each country:

country	| cases
----- |  ----- 
China |	84024
Thailand |	3017
... | ...

The following geo_subsampler command subsamples the input tree according to the case proportions and (as much as possible) equally between the months,
in order to keep 1000 tips:

```bash
geo_subsampler --tree example_data/covid.nwk --metadata example_data/metadata.tab \
--location_column country --date_column "sampling date" --cases example_data/cases.tab \
--output_dir example_data/results --size 1000
```

The resulting tree is put into [example_data/results](example_data/results) folder:
([covid.subsampled.0.nwk](example_data/results/covid.subsampled.0.nwk)). This folder also contains the ids of the tips retained in the subsampled tree:
([covid.subsampled.0.ids](example_data/results/covid.subsampled.0.ids)), and two tables with the statistics on the subsampling:
[case_counts.tab](example_data/results/case_counts.tab) and [case_counts_per_time.tab](example_data/results/case_counts_per_time.tab).


## Detailed options
- **--tree TREE**           Path to the input phylogeny (NOT time-scaled) in newick format.
- **--metadata METADATA**   Path to the metadata table containing location and date annotations, in a tab-delimited format.
- **--sep SEP**             Separator used in the metadata and case tables. By default a tab-separated table is assumed.
- **--index_column INDEX_COLUMN**
                        number (starting from zero) of the index column (containing tree tip names) in the metadata table. By default is the first column (corresponding to the number 0)
- **--location_column LOCATION_COLUMN**
                        name of the column containing location annotations in the metadata table.
- **--date_column DATE_COLUMN**
                        name of the column containing date annotations in the metadata table.
- **--cases CASES**         A tab-separated file with two columns. The first column lists the locations, while the second column contains the numbers of declared cases or proportions for the
                        corresponding locations
- **--start_date START_DATE**
                        If specified, all the cases before this date will be included in all the sub-sampled data sets.
- **--size SIZE**           Target size of the sub-sampled data set (in number of samples). By default, will be set to a half of the data set represented by the input tree.
- **--repetitions REPETITIONS** Number of sub-sampled trees to produce. By default 1.
- **--output_dir OUTPUT_DIR**
                        Path to the directory where the sub-sampled results should be saved.
- **--min_cases MIN_CASES**
                        Minimum number of samples to retain for each location.
- **--date_precision {year,month,day}**
                        Precision for homogeneous subsampling over time within each location. By default (month) will aim at distributing selected location samples equally over months.

