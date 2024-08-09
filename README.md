# GEO subsampler

The script subsamples a given phylogenetic tree to rebalance the samples at different locations 
according to user-specified proportions. For each country the cases are picked to be uniformly spread over sampling period.
With these constraints in mind, the script uses phylogenetic diversity to pick the samples to remove.


## Script parameters
- --tree TREE           Path to the input phylogeny (rooted and NOT time-scaled) in newick format.
- --metadata METADATA   Path to the metadata table containing location and date annotations, in a tab-delimited format.
- --sep SEP             Separator used in the metadata and case tables. By default a tab-separated table is assumed.
- --index_column INDEX_COLUMN
                        number (starting from zero) of the index column (containing tree tip names) in the metadata table. By default is the first column (corresponding to the number 0)
- --location_column LOCATION_COLUMN
                        name of the column containing location annotations in the metadata table.
- --date_column DATE_COLUMN
                        name of the column containing date annotations in the metadata table.
- --cases CASES         A tab-separated file with two columns. The first column lists the locations, while the second column contains the numbers of declared cases or proportions for the
                        corresponding locations
- --start_date START_DATE
                        If specified, all the cases before this date will be included in all the sub-sampled data sets.
- --size SIZE           Target size of the sub-sampled data set (in number of samples). By default, will be set to a half of the data set represented by the input tree.
- --repetitions REPETITIONS Number of sub-samplings to perform. By default 1.
- --output_dir OUTPUT_DIR
                        Path to the directory where the sub-sampled results should be saved.
- --min_cases MIN_CASES
                        Minimum number of samples to retain for each location.
- --date_precision {year,month,day}
                        Precision for homogeneous subsampling over time within each location. By default (month) will aim at distributing selected location samples equally over months.


## Installation
To install geo_subsampler:
```bash
pip3 install geo_subsampler
```