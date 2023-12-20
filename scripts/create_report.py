import pandas as pd
import json
import click

@click.group
def manta():
    pass

@manta.command
@click.option("-i", "input_file", prompt="the input file", help="The csv file for manta", required=True, type=click.Path(exists=True))
@click.option("-a", "adiv_file", prompt="the alpha diversity file", help="The tsv file for alpha diversity", required=True, type=click.Path(exists=True))
@click.option("-p", "taxonpath", prompt="the taxonpath file", help="The path to the taxonpath.json file", required=True, type=click.Path(exists=True))
@click.option("-n", "names", prompt="the taxonomy names file", help="The path to the names.json file", required=True, type=click.Path(exists=True))
@click.option("-o", "output_file", prompt="the output file", help="The path to the report file", required=True, type=click.Path(exists=False))
def summarize(input_file, taxonpath, names, adiv_file, output_file):
    '''Generate ready-for-manta format data'''
    with open(taxonpath) as f:
        taxonpath=json.load(f)

    with open(names) as f:
        names=json.load(f)

    rank_codes = ['p','c','o','f','g']
    def get_taxonstring(tid):
        tp = taxonpath[tid]
        cr = tp['rank'][0:1]
        list = []
        for r in rank_codes:
            if tp[r] == "":
                list.append("")
            else:
                list.append(f'{r}_{names[tp[r]]}')
            if r == cr:
                break
        return ";".join(list)

    df = pd.read_csv(input_file)
    all_ranks = ['phylum', 'class', 'order', 'family', 'genus']
    with pd.ExcelWriter(output_file) as writer:
        pd.read_csv(adiv_file, sep="\t", skiprows=[1], index_col=['Sample ID']).to_excel(writer, sheet_name='alpha_div')
        for rank in all_ranks:
            out = df.pivot_table(index='sample_id', columns=f'{rank}_id', values='read_pct', aggfunc=sum, fill_value=0)
            cdic = {t:get_taxonstring(t) for t in out.columns.drop('uc')}
            cdic['uc'] = 'unclassified'
            out.rename(columns=cdic, inplace=True)
            clist = out.columns.to_list()
            clist.sort()
            out[clist].to_excel(writer, sheet_name=rank)

if __name__ == '__main__':
    manta()
