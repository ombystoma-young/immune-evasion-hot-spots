import argparse

import numpy as np
import pandas as pd


PRELIM_INFO = {
    'prot': ['NC_001604.1.5822',
             'NC_001604.1.3100',
             'NC_001604.1.1972',
             'NC_001604.1.1797',
             'NC_001604.1.1639',
             'NC_001604.1.1433',
             'NC_001604.1.1278',

             'NC_003298.1.5630',
             'NC_003298.1.2905',
             'NC_003298.1.1781',
             'NC_003298.1.1627',
             'NC_003298.1.1359',
             'OL964749.1.1549',
             'OL964748.1.1549'
             ],
    'prelim_info': ['T7pRNAP',
                    'T7pSTPK',
                    'T7p0.6B',
                    'T7p0.6A',
                    'T7p0.5',
                    'T7p0.4',
                    'Ocr',

                    'T3pRNAP',
                    'T3pSTPK',
                    'T3p0.65',
                    'T3p0.6',
                    'SAMase',
                    'Ocr::STPK',
                    'Ocr::STPK'
                    ]
}


def parse_args():
    parser = argparse.ArgumentParser(description='Part of "snake_process_clustering_results" pipeline. '
                                                 'Read information about clustering by MMSeqs2')
    parser.add_argument('--clu', default=None, type=str, nargs='?',
                        help='path to input table file with clusters')
    parser.add_argument('--phrogs', default=None, type=str, nargs='?',
                        help='path to tsv file with found phrogs domains')
    parser.add_argument('--pfam', default=None, type=str, nargs='?',
                        help='path to tsv file with found pfam domains')
    parser.add_argument('--apis', default=None, type=str, nargs='?',
                        help='path to tsv file with found apis domains')
    parser.add_argument('--phys', default=None, type=str, nargs='?',
                        help='path to tsv file with info about '
                             'physical properties of sequences')
    parser.add_argument('--clans', default=None, type=str, nargs='?',
                        help='path to tsv file with information about clans')
    parser.add_argument('--output', default=None, type=str, nargs='?',
                        help='path to output table file')
    parser.add_argument('--long', default=None, type=str, nargs='?',
                        help='path to output table file without grouping by clusters')
    return parser.parse_args()


def read_table(in_path: str) -> pd.DataFrame:
    READ_PARAMS = {
        clu_descr: {'skiprows': 1,
                    'colnames': ['clu', 'prot'],
                    'sep': '\t'},
        phrogs: {'skiprows': 1,
                 'colnames': ['prot', 'phrog', 'annot', 'phrog_category'],
                 'sep': '\t'},
        pfam: {'skiprows': 1,
               'colnames': ['prot', 'pfam_dom', 'pfam_id', 'pfam_full_seq_e', 'pfam_best_one_domain_e',
                            'pfam_description_of_target'],
               'sep': '\t'},
        apis: {'skiprows': 1,
               'colnames': ['prot', 'apis', 'dbapis_clan_id', 'apis_genes', 'defense_systems'],
               'sep': '\t'},
        phys_char: {'skiprows': 0,
                    'colnames': ['prot', 'seq', 'length', 'pi'],
                    'sep': '\t'},
        clans: {'skiprows': 0,
                'colnames': ['clu', 'clan'],
                'sep': '\t'}
    }
    df = pd.read_csv(in_path, sep=READ_PARAMS[in_path]['sep'],
                     skiprows=READ_PARAMS[in_path]['skiprows'],
                     names=READ_PARAMS[in_path]['colnames'])
    return df


def preprocess_pfam_df(source_df: pd.DataFrame) -> pd.DataFrame:
    df = source_df.copy()
    df = df[['prot', 'pfam_dom', 'pfam_id']]
    df = df.groupby('prot').agg(lambda x: ','.join(x))
    return df.reset_index()


def unite_info_about_each_protein(clu_descr_path: str,
                                  phrogs_path: str,
                                  pfam_path: str,
                                  apis_path: str,
                                  phys_char_path: str) -> pd.DataFrame:
    # read tables
    prelim_info = pd.DataFrame.from_dict(PRELIM_INFO).set_index('prot')
    clu_descr_df = read_table(clu_descr_path)
    phrogs_df = read_table(phrogs_path)
    pfam_df = preprocess_pfam_df(read_table(pfam_path))
    apis_df = read_table(apis_path)
    phys_char_df = read_table(phys_char_path)
    # perform chain of left joins, based on clu_descr info dataframe
    agg_df = clu_descr_df.set_index('prot').join(prelim_info)
    agg_df = (agg_df
              .join(phrogs_df.set_index('prot'))
              .join(pfam_df.set_index('prot'))
              .join(apis_df.set_index('prot'))
              .join(phys_char_df.set_index('prot'))
              ).reset_index()
    return agg_df


def add_clans_info(clans_path: str, source_agg_df: pd.DataFrame) -> pd.DataFrame:
    clans_df = read_table(clans_path)
    clans_df['clu'] = clans_df.apply(lambda x: int(x['clu'].split('_')[1]), axis=1)
    agg_df = source_agg_df.set_index('clu').join(clans_df.set_index('clu'))
    return agg_df


def process_prelim_info(x):
    x = set(x)
    if all(el is np.nan for el in x):
        return np.nan
    else:
        return {el for el in x if el is not np.nan}


def process_cols2counts(x):
    counts = x.value_counts(dropna=False).to_dict()
    return counts


def process_clan_info(x):
    x = set(x)
    if all(el is np.nan for el in x):
        return np.nan
    else:
        return list(x)[0]


def group_by_clu(source_agg_df: pd.DataFrame) -> pd.DataFrame:
    grouped_agg_df = source_agg_df.groupby(['clu'], dropna=False).agg(
        clan=pd.NamedAgg(column='clan',
                         aggfunc=process_clan_info),
        num_reprs=pd.NamedAgg(column="prot",
                              aggfunc="count"),
        prelim_info=pd.NamedAgg(column="prelim_info",
                                aggfunc=process_prelim_info),
        annotations=pd.NamedAgg(column='annot',
                                aggfunc=process_cols2counts),
        pfam_domains=pd.NamedAgg(column='pfam_dom',
                                 aggfunc=process_cols2counts),
        mean_length=pd.NamedAgg(column='length',
                                aggfunc='mean'),
        std_length=pd.NamedAgg(column='length',
                               aggfunc='std'),
        mean_pi=pd.NamedAgg(column='pi',
                            aggfunc='mean'),
        std_pi=pd.NamedAgg(column='pi',
                           aggfunc='std'),
        apises=pd.NamedAgg(column='apis',
                           aggfunc=process_cols2counts),
        dbapis_clans=pd.NamedAgg(column='dbapis_clan_id',
                                 aggfunc=process_cols2counts),
        defense_systems=pd.NamedAgg(column='defense_systems',
                                    aggfunc=process_cols2counts),
        phrogs=pd.NamedAgg(column='phrog',
                           aggfunc=process_cols2counts),
        pfam_ids=pd.NamedAgg(column='pfam_id',
                             aggfunc=process_cols2counts)
    )
    grouped_agg_sorted = (grouped_agg_df
                          .reset_index().
                          sort_values(by='num_reprs',
                                      ascending=False,
                                      ignore_index=True))
    grouped_agg_sorted.index = np.arange(1, len(grouped_agg_sorted) + 1)
    grouped_agg_sorted.reset_index(inplace=True)
    grouped_agg_sorted = grouped_agg_sorted.rename(columns={'clu': 'old_clu', 'index': 'clu'})
    grouped_agg_sorted['old_clu'] = grouped_agg_sorted['old_clu'].apply(lambda x: f'clu_{x}')
    grouped_agg_sorted['new_clan'] = grouped_agg_sorted.apply(
        lambda x: x['clan'] if x['clan'] is not np.nan else f'clan_00{x["clu"]}', axis=1)
    return grouped_agg_sorted


if __name__ == '__main__':
    clu_descr = parse_args().clu
    phrogs = parse_args().phrogs
    pfam = parse_args().pfam
    apis = parse_args().apis
    phys_char = parse_args().phys
    clans = parse_args().clans

    res_tab = parse_args().output
    res_tab_long = parse_args().long
    # add protein properties
    agg_prot_df = unite_info_about_each_protein(clu_descr_path=clu_descr,
                                                phrogs_path=phrogs,
                                                pfam_path=pfam,
                                                apis_path=apis,
                                                phys_char_path=phys_char
                                                )
    # add clusters properties
    agg_prot_df = add_clans_info(clans, agg_prot_df)
    agg_prot_df.to_csv(res_tab_long, sep='\t')
    agg_grouped_prot_df = group_by_clu(agg_prot_df)
    agg_grouped_prot_df.to_csv(res_tab, sep='\t', index=False)
