# GEO subsampler

The script subsamples a given phylogenetic tree to rebalance the samples at different locations 
according to user-specified proportions. For each country the cases are picked to be uniformly spread over sampling period.
With these constraints in mind, the script uses phylogenetic diversity to pick the samples to remove.


## Script parameters
- --tree TREE.nwk -- the input phylogenetic tree in newick format, must __not__ be dated and must be rooted. 
- --size N -- the number of samples to retain in the subsampled tree
- --repetitions R -- the number of subsampled trees to output
-  --output_dir NAME -- the path to the directory where the subsampled trees will be stored


- 
- pandas, numpy, ete3 и hdx-python-country (например, с помощью pip3 install pandas numpy ete3 hdx-python-country). Дальше можно запускать скрипт вот так:

python3 phylogenetic_diversity.py --nwk MPXV_gisaid_20221111_tree.nwk --cases cases.tab --repetitions 10 --output_dir subsampling
 

## Installation
To install geo_subsampler:
```bash
pip3 install geo_subsampler
```