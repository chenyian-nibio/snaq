import pandas as pd
import json
import click

rank_empty_tid = {
    'k': '', 'p':'', 'c':'', 'o':'', 'f':'','g':'','s':''
}

def map_name_to_taxonbid(taxon_name_arr_list, bacteria, non_bacteria, taxonpath):
    # to prevent some identical names between bacteria and fungi, we need two different dictionary
    taxon_id_arr_list = [[bacteria.get(x[3:]) if y[0] == 'd__Bacteria' else non_bacteria.get(x[3:]) for x in y] for y in taxon_name_arr_list]
    ret = []
    for tid_list in taxon_id_arr_list:
        if all(id is None for id in tid_list):
            continue
        else:
            max_rank_tid = tid_list[max([i for i in range(len(tid_list)) if tid_list[i] is not None])]
            txp = taxonpath.get(max_rank_tid, rank_empty_tid)
            item = []

            for r in ['k', 'p', 'c', 'o', 'f', 'g', 's']:
                _item = txp[r]
                if _item == "":
                    _item = "uc"
                item.append(_item)
            ret.append(item)
    return ret

def top_taxons(df):
    ret = df.sort_values(9, ascending=False)
    ret['cumsum'] = ret[9].cumsum()
    # Include at least the first row even if it exceeds 90
    ret = ret[(ret['cumsum'] < 90) | (ret['cumsum'] == ret[9].iloc[0])]
    ret.drop(['cumsum', 8,9,10], axis="columns", inplace=True)
    ret = ret.melt(id_vars=[0, 11])
    ret.drop_duplicates(inplace=True)
    ret.columns = ['sample_id', 'method_id','rank_id', 'taxonomy_id']
    ret = ret[ret['taxonomy_id']!= 'uc']
    return ret[['sample_id', 'rank_id', 'taxonomy_id', 'method_id']]

@click.command()
@click.option("-i", "input_file", required=True, type=str)
@click.option("-o", "output_file", required=True, type=str)
@click.option("-t", "taxonpath", required=True, type=str)
@click.option("-s", "sample_file_name", required=True, type=str)
@click.option("-a", "abundant_taxonomy", required=True, type=str)
@click.option("-n", "names", required=True, type=str)
@click.option("-d", "database", required=True, type=str)
@click.option("-x", "output_taxonomy", required=True, type=str)
def manta(input_file, output_file, taxonpath, abundant_taxonomy, sample_file_name, names, database, output_taxonomy):
    df = pd.read_csv(input_file, sep="\t", skiprows=[0])
    # get rid of Mitochondria and Chloroplast 
    df = df[~df['#OTU ID'].str.contains('g__Mitochondria') & ~df['#OTU ID'].str.contains('g__Chloroplast')].reset_index(drop=True)

    with open(taxonpath) as f:
        taxonpath=json.load(f)
    with open(names) as f:
        names=json.load(f)

    all_reads = {}
    for sid in df.columns.drop("#OTU ID"):
        all_reads[sid] = df[sid].sum()

    taxon_name_arr_list = [x.split(";") for x in df['#OTU ID']]

    # we probably don't need those taxon ids that don't belong to any hierachy
    subset = {k:v for k, v in names.items() if k in taxonpath}
    # to prevent some identical names between bacteria and fungi, we need two different dictionary
    bacteria = {v:k for k, v in subset.items() if taxonpath[k]['k'] == '2'}
    non_bacteria = {v:k for k, v in subset.items() if taxonpath[k]['k'] != '2'}

    # synonyms for the legacy taxon name
    bacteria['Paludicola'] = '2038676'

    df_taxonid = pd.DataFrame(map_name_to_taxonbid(taxon_name_arr_list, bacteria, non_bacteria, taxonpath))
    df_taxonid.columns = ['0', '1', '2', '3', '4', '5', '6']
    df.drop("#OTU ID", axis="columns", inplace=True)
    df3 = df_taxonid.join(df)

    df3m = df3.melt(id_vars=['0', '1', '2', '3', '4', '5', '6'])
    df4 = df3m.loc[:,['0', '1', '2', '3', '4', '5', '6']].copy().melt()

    df3m = df3m.loc[:,['variable', '0', '1', '2', '3', '4', '5', '6', 'value']]

    df3m['pct'] = (df3m['value']/df3m['variable'].apply(lambda x: all_reads[x]))*100

    df3m['db'] = int(database)
    df3m['method_id'] = 1
    df3m = df3m[df3m['value']!= 0]
    df3m.rename(columns={'0':"kingdom_id", '1':'phylum_id', '2':'class_id', '3':'order_id', '4':'family_id', '5':'genus_id', '6':'species_id', 
                         'value':'read_num', 'variable':'sample_id', 'pct':'read_pct', 'db':'reference_db_id'}, inplace=True)
    df3m['read_num'] = df3m['read_num'].astype(int)
    df3m.to_csv(output_file, index=False)

    df3m['sample_id'].rename("id").drop_duplicates().to_csv(sample_file_name, index=False)

    df4['variable'] = df4['variable'].astype(int) + 1
    df4 = df4[df4['value']!= "uc"]
    df4['names'] = [names.get(x) for x in df4['value']]
    df4=df4.iloc[:,[1,0,2]].drop_duplicates()
    df4.columns = ['id', "rank_id", "name"]
    df4.to_csv(output_taxonomy, index=False)

    # abundant_taxons
    df5 = df3m.copy()
    df5.columns = range(12) # type: ignore
    abundant_taxons = pd.DataFrame(df5.groupby(0).apply(lambda x: top_taxons(x))).reset_index(drop=True)
    abundant_taxons.to_csv(abundant_taxonomy, index=False)


if __name__ == "__main__":
    manta()
