import logging
import os
from collections import defaultdict

import numpy as np
import pandas as pd
from ete3 import Tree

MIN_CASES = 'min_to_take'

LOC_DATE = 'loc_date'

CASES = 'cases'

SAMPLED_CASES = 'sampled_cases'

SUBSAMPLED_CASES = 'subsampled_cases'

DAY = 'day'

MONTH = 'month'

YEAR = 'year'


def read_tree(tree_path):
    tree = None
    for f in (3, 2, 5, 0, 1, 4, 6, 7, 8, 9):
        try:
            tree = Tree(tree_path, format=f, quoted_node_names=True)
            break
        except:
            continue
    if not tree:
        raise ValueError('Could not read the tree {}. Is it a valid newick?'.format(tree_path))
    return tree


def remove_certain_leaves(tree, to_remove=lambda node: False):
    """
    Removes all the branches leading to leaves identified positively by to_remove function.
    :param tree: the tree of interest (ete3 Tree)
    :param to_remove: a method to check if a leaf should be removed.
    :return: void, modifies the initial tree.
    """

    tips = [tip for tip in tree if to_remove(tip)]
    for node in tips:
        if node.is_root():
            return None
        parent = node.up
        parent.remove_child(node)
        # If the parent node has only one child now, merge them.
        if len(parent.children) == 1:
            brother = parent.children[0]
            brother.dist += parent.dist
            if parent.is_root():
                brother.up = None
                tree = brother
            else:
                grandparent = parent.up
                grandparent.remove_child(parent)
                grandparent.add_child(brother)
    return tree


def calc_size_stats(df, N):
    undersampled_locs = sorted(df[(df['sampled_cases'] / df['frequencies']) < N].index)
    logging.info('{} locations ({}) are undersampled for the selected threshold {}.'
                 .format(len(undersampled_locs), ', '.join(undersampled_locs), N))


def subsample_by_phylogenetic_diversity(tree, df, sampled_case_per_time_df, to_keep):
    bin2tips = defaultdict(list)
    for _ in tree:
        if _.name not in to_keep and sampled_case_per_time_df.loc[df.loc[_.name, LOC_DATE], 'remove_cases']:
            bin2tips[_.dist].append(_)
    while sampled_case_per_time_df['remove_cases'].sum():
        min_bin = min(bin2tips.keys())
        min_tips = bin2tips[min_bin]
        tip = np.random.choice(min_tips, 1)[0]
        min_tips.remove(tip)
        if not min_tips:
            del bin2tips[min_bin]
        if sampled_case_per_time_df.loc[df.loc[tip.name, LOC_DATE], 'remove_cases'] > 0:
            sampled_case_per_time_df.loc[df.loc[tip.name, LOC_DATE], 'remove_cases'] -= 1
            parent = tip.up
            parent.remove_child(tip)
            if len(parent.children) == 1 and not parent.is_root():
                grandparent = parent.up
                child = parent.remove_child(parent.children[0])
                if child.is_leaf() and child.dist in bin2tips and child in bin2tips[child.dist]:
                    bin2tips[child.dist].remove(child)
                    if not bin2tips[child.dist]:
                        del bin2tips[child.dist]
                    bin2tips[parent.dist + child.dist].append(child)
                grandparent.add_child(child, dist=parent.dist + child.dist)
                grandparent.remove_child(parent)
    while len(tree.children) == 1:
        tree = tree.children[0]
        tree.dist = 0
        tree.up = None
    return tree


def main():
    """
    Entry point for tree subsampling with command-line arguments.
    :return: void
    """
    import argparse

    parser = argparse.ArgumentParser(description='Subsampling of phylogenetic trees.')

    parser.add_argument('--tree', required=True, type=str,
                        help='Path to the input phylogeny (NOT time-scaled) in newick format.')
    parser.add_argument('--metadata', required=True, type=str,
                        help='Path to the metadata table containing location and date annotations, '
                             'in a tab-delimited format.')
    parser.add_argument('--index_column', type=int, default=0,
                        help='Number (starting from zero) of the index column (containing tree tip names) '
                             'in the metadata table. '
                             'By default is the first column (corresponding to the number 0)')
    parser.add_argument('--location_column', type=str, default='country',
                        help='name of the column containing location annotations in the metadata table.')
    parser.add_argument('--date_column', type=str, default='date',
                        help='name of the column containing date annotations in the metadata table.')
    parser.add_argument('--cases', required=True, type=str,
                        help='Path to the case count table, in a tab-separated format, with two columns. '
                             'The first column lists the locations, '
                             'while the second column contains the numbers of declared cases or proportions '
                             'for the corresponding locations')
    parser.add_argument('--sep', type=str, default='\t',
                        help='Separator used in the metadata and case tables. '
                             'By default tab-separated tables are assumed.')
    parser.add_argument('--start_date', default=None, type=str,
                        help='If specified, all the cases before this date '
                             'will be included in all the sub-sampled data sets.')
    parser.add_argument('--size', default=0, type=int,
                        help='Target size of the sub-sampled data set (in number of samples). '
                             'By default, will be set to a half of the data set represented by the input tree.')
    parser.add_argument('--repetitions', type=int, default=1,
                        help='Number of sub-samplings to perform. By default 1.')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Path to the directory where the sub-sampled results should be saved.')
    parser.add_argument('--min_cases', type=int, default=0,
                        help='Minimum number of samples to retain for each location.')
    parser.add_argument('--date_precision', default=MONTH, choices=[YEAR, MONTH, DAY],
                        help='Precision for homogeneous subsampling over time within each location. '
                             'By default (month) will aim at distributing selected location samples equally over months.')
    params = parser.parse_args()

    logging.getLogger().setLevel(level=logging.INFO)

    os.makedirs(params.output_dir, exist_ok=True)

    tree = read_tree(params.tree)
    tree_ids = {_.name for _ in tree}

    if not params.size:
        params.size = int(len(tree) * 0.5)

    df, tree_ids = parse_metadata(tree_ids, csv=params.metadata, sep=params.sep, index_col=params.index_column,
                                  loc_col=params.location_column, date_col=params.date_column,
                                  date_precision=params.date_precision)

    min_date_to_subsample = pd.to_datetime(params.start_date) if params.start_date is not None \
        else min(df[params.date_column])
    to_keep = list(df[df[params.date_column] < min_date_to_subsample].index)
    kept_loc_df = df.loc[to_keep, :]

    sampled_case_per_time_df = df[[LOC_DATE, params.date_column]].groupby([LOC_DATE]).count()
    sampled_case_per_time_df.columns = [SAMPLED_CASES]

    sampled_case_df = df[[params.location_column, params.date_column]].groupby([params.location_column]).count()
    sampled_case_df.columns = [CASES]

    case_df = parse_cases(params.cases, params.sep, sampled_case_df, params.size)

    case_df[MIN_CASES] = np.minimum(case_df[SAMPLED_CASES], params.min_cases)

    kept_case_df = df.loc[to_keep, [params.location_column, params.date_column]] \
        .groupby([params.location_column]).count()
    kept_case_df.columns = [CASES]

    case_df.loc[kept_case_df.index, MIN_CASES] = (
        np.maximum(case_df.loc[kept_case_df.index, MIN_CASES], kept_case_df[CASES]))
    case_df[SUBSAMPLED_CASES] = case_df[MIN_CASES]

    all_locs = set(case_df.index)
    processed_locs = set()

    while True:
        todo_locs = list(all_locs - processed_locs)
        freq_df = case_df.loc[todo_locs, [CASES]]
        freq_df['frequencies'] = freq_df[CASES] / freq_df[CASES].sum()
        case_df.loc[todo_locs, SUBSAMPLED_CASES] = \
            np.round(freq_df['frequencies']
                     * (params.size - case_df.loc[list(processed_locs), SUBSAMPLED_CASES].sum()), 0).astype(int)
        # If the min_cases is larger than the number of cases corresponding to this location's proportion,
        #   let's just take the min number of cases and be done with it.
        # If the number of target subsampled cases is larger than the total number of sampled cases for this location,
        #   let's just take them all and be done with it.
        new_processed_locs = (set(case_df[case_df[SUBSAMPLED_CASES] <= case_df[MIN_CASES]].index) |
                              set(case_df[case_df[SUBSAMPLED_CASES] >= case_df[SAMPLED_CASES]].index))
        case_df[SUBSAMPLED_CASES] = np.minimum(np.maximum(case_df[SUBSAMPLED_CASES], case_df[MIN_CASES]),
                                               case_df[SAMPLED_CASES])
        if processed_locs == new_processed_locs or new_processed_locs == all_locs \
                or case_df[SUBSAMPLED_CASES].sum() >= params.size:
            break
        processed_locs = new_processed_locs

    while case_df[SUBSAMPLED_CASES].sum() < params.size:
        loc = np.random.choice([_ for _ in all_locs
                                if case_df.loc[_, SUBSAMPLED_CASES] < case_df.loc[_, SAMPLED_CASES]], 1)[0]
        case_df.loc[loc, SUBSAMPLED_CASES] += 1

    while case_df[SUBSAMPLED_CASES].sum() > max(params.size, len(to_keep)):
        loc = np.random.choice([_ for _ in all_locs
                                if case_df.loc[_, SUBSAMPLED_CASES] > case_df.loc[_, MIN_CASES]], 1)[0]
        case_df.loc[loc, SUBSAMPLED_CASES] -= 1

    case_df[[CASES, SAMPLED_CASES, SUBSAMPLED_CASES]] \
        .to_csv(os.path.join(params.output_dir, 'case_counts.tab'), sep='\t', index_label=params.location_column)

    sampled_case_per_time_df[SUBSAMPLED_CASES] = 0
    sampled_case_per_time_df[MIN_CASES] = 0

    for location in case_df.index:
        loc_date_count_df = df[df[params.location_column] == location]
        loc_date_count_df = loc_date_count_df[[LOC_DATE, params.date_column]].groupby([LOC_DATE]).count()
        loc_date_count_df.columns = [CASES]
        loc_date_count_df[MIN_CASES] = (
            loc_date_count_df.index.map(lambda _: len(kept_loc_df[kept_loc_df[LOC_DATE] == _])))
        all_locdates = list(loc_date_count_df.index)
        target = case_df.loc[location, SUBSAMPLED_CASES]

        if case_df.loc[location, SAMPLED_CASES] == target:
            sampled_case_per_time_df.loc[all_locdates, SUBSAMPLED_CASES] \
                = sampled_case_per_time_df.loc[all_locdates, SAMPLED_CASES]
            continue

        sampled_case_per_time_df.loc[loc_date_count_df.index, MIN_CASES] = loc_date_count_df[MIN_CASES]
        sampled_case_per_time_df[SUBSAMPLED_CASES] \
            = np.maximum(sampled_case_per_time_df[SUBSAMPLED_CASES], sampled_case_per_time_df[MIN_CASES])

        already_took = loc_date_count_df[MIN_CASES].sum()
        if already_took != target:
            n_ld = len(loc_date_count_df)

            all_locdates = set(all_locdates)
            processed_locdates = set()

            while True:
                todo_locdates = list(all_locdates - processed_locdates)
                freq = 1 / n_ld
                sampled_case_per_time_df.loc[todo_locdates, SUBSAMPLED_CASES] = (
                    np.round(freq
                             * (target - sampled_case_per_time_df.loc[
                        list(processed_locdates), SUBSAMPLED_CASES].sum()),
                             0).astype(int))
                new_processed_locdates = \
                    ((set(sampled_case_per_time_df[sampled_case_per_time_df[SUBSAMPLED_CASES]
                                                   <= sampled_case_per_time_df[MIN_CASES]].index) |
                      set(sampled_case_per_time_df[sampled_case_per_time_df[SUBSAMPLED_CASES]
                                                   >= sampled_case_per_time_df[SAMPLED_CASES]].index))
                     & all_locdates)
                sampled_case_per_time_df[SUBSAMPLED_CASES] = (
                    np.minimum(np.maximum(sampled_case_per_time_df[SUBSAMPLED_CASES],
                                          sampled_case_per_time_df[MIN_CASES]),
                               sampled_case_per_time_df[SAMPLED_CASES]))
                if processed_locdates == new_processed_locdates or new_processed_locdates == all_locdates \
                        or sampled_case_per_time_df.loc[list(all_locdates), SUBSAMPLED_CASES].sum() >= target:
                    break
                processed_locdates = new_processed_locdates
            while sampled_case_per_time_df.loc[list(all_locdates), SUBSAMPLED_CASES].sum() < target:
                locdate = np.random.choice([_ for _ in all_locdates
                                            if sampled_case_per_time_df.loc[_, SUBSAMPLED_CASES]
                                            < sampled_case_per_time_df.loc[_, SAMPLED_CASES]], 1)[0]
                sampled_case_per_time_df.loc[locdate, SUBSAMPLED_CASES] += 1

            while sampled_case_per_time_df.loc[list(all_locdates), SUBSAMPLED_CASES].sum() > target:
                locdate = np.random.choice([_ for _ in all_locdates
                                            if sampled_case_per_time_df.loc[_, SUBSAMPLED_CASES]
                                            > sampled_case_per_time_df.loc[_, MIN_CASES]], 1)[0]
                sampled_case_per_time_df.loc[locdate, SUBSAMPLED_CASES] -= 1

    sampled_case_per_time_df[[SAMPLED_CASES, SUBSAMPLED_CASES]] \
        .to_csv(os.path.join(params.output_dir, 'case_counts_per_time.tab'), sep='\t', index_label=LOC_DATE)

    tree_name = os.path.splitext(os.path.basename(params.tree))[0]
    for i in range(params.repetitions):
        sampled_case_per_time_df['remove_cases'] = sampled_case_per_time_df[SAMPLED_CASES] \
                                                   - sampled_case_per_time_df[SUBSAMPLED_CASES]
        tree = remove_certain_leaves(read_tree(params.tree), lambda _: _.name not in tree_ids)
        tree = subsample_by_phylogenetic_diversity(tree, df, sampled_case_per_time_df, to_keep)

        with open(os.path.join(params.output_dir, '{}.subsampled.{}.ids'.format(tree_name, i)), 'w+') as f:
            f.write('\n'.join(_.name for _ in tree))
        tree.write(outfile=os.path.join(params.output_dir, '{}.subsampled.{}.nwk'.format(tree_name, i)))


def parse_cases(csv, sep, sampled_case_df, size):
    if csv:
        case_df = pd.read_csv(csv, sep=sep, index_col=0)
        case_df.columns = [CASES]

        # make sure that the locations in the declared and sampled case tables are the same
        for _ in set(sampled_case_df.index) - set(case_df.index):
            logging.warning('No cases declared for {}, though some samples are present in the tree'.format(_))
            case_df.loc[_, CASES] = 0
        case_df = case_df.loc[sampled_case_df.index, :]
        if case_df[CASES].sum() < 2:
            logging.warning('Seems like the proportions are given instead of declared case counts, '
                            'we then gonna generate a case count according to these proportions.')
            case_df[CASES] *= len(sampled_case_df) * sampled_case_df[CASES].sum()
            case_df[CASES] = case_df[CASES].astype(int)
        # Update the declared cases with the sampled case figure if there are more sampled cases than the declared ones
        case_df = pd.concat([case_df, sampled_case_df]).groupby(level=0).max()
    else:
        case_df = pd.DataFrame(index=sampled_case_df.index)
        case_df[CASES] = sampled_case_df[CASES].sum()
    case_df['frequencies'] = case_df[CASES] / case_df[CASES].sum()
    case_df[SAMPLED_CASES] = sampled_case_df.loc[case_df.index, CASES]
    calc_size_stats(case_df, size)
    return case_df


def parse_metadata(tree_ids, csv, sep, index_col, loc_col, date_col, date_precision):
    df = pd.read_csv(csv, sep=sep, index_col=index_col)[[loc_col, date_col]]
    df.index = df.index.map(str)
    ids_wo_metadata = tree_ids - set(df.index)
    if len(ids_wo_metadata) == len(tree_ids):
        raise ValueError('Could not find annotations for the tree tips in the metadata table, '
                         'please check their format.')
    if len(ids_wo_metadata):
        logging.warning('{} of the tree tips do not have annotations in the metadata table, we will remove them'
                        .format(len(ids_wo_metadata)))
        tree_ids -= ids_wo_metadata
    df = df.loc[list(tree_ids), :]
    df[date_col] = pd.to_datetime(df[date_col].astype(str).str.replace('.0', '', regex=False),
                                  infer_datetime_format=True, errors='coerce')
    ids_wo_dates = set(df[pd.isna(df[date_col])].index)
    if len(ids_wo_dates) == len(tree_ids):
        raise ValueError('Could not parse the dates in the metadata table. Consider changing the format, '
                         'e.g. to YYYY-MM-DD, or YYYY-MM, or YYYY.')
    if ids_wo_dates:
        logging.warning('{} tree tips do not have date information, we will remove them'
                        .format(len(ids_wo_dates)))
        tree_ids -= ids_wo_dates
        df = df.loc[list(tree_ids), :]
    ids_wo_locs = set(df[pd.isna(df[loc_col])].index)
    if len(ids_wo_locs) == len(tree_ids):
        raise ValueError('Could not parse the locations in the metadata table, '
                         'check the metadata and location column name ({}).'.format(loc_col))
    if ids_wo_locs:
        logging.warning('{} tree tips do not have location information, we will remove them'
                        .format(len(ids_wo_locs)))
        tree_ids -= ids_wo_locs
        df = df.loc[list(tree_ids), :]

    def stratify_date(d):
        d = d.strftime('%Y-%m-%d')
        if date_precision == MONTH:
            d = d[:7]
        elif date_precision == YEAR:
            d = d[:4]
        return d

    df['loc_date'] = df[loc_col] + df[date_col].apply(stratify_date)

    return df, tree_ids


if '__main__' == __name__:
    main()
