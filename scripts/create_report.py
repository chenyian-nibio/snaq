import pandas as pd
import pandas.io.formats.excel
import json
import click
import xlsxwriter
import xlsxwriter.worksheet

rank_codes = ['p','c','o','f','g']
def get_taxonstring(tid, taxonpath, names):
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

@click.group
def create():
    pass

@create.command
@click.option("-i", "input_file", prompt="the input file", help="The csv file for manta", required=True, type=click.Path(exists=True))
@click.option("-a", "adiv_file", prompt="the alpha diversity file", help="The tsv file for alpha diversity", required=True, type=click.Path(exists=True))
@click.option("-p", "taxonpath", prompt="the taxonpath file", help="The path to the taxonpath.json file", required=True, type=click.Path(exists=True))
@click.option("-n", "names", prompt="the taxonomy names file", help="The path to the names.json file", required=True, type=click.Path(exists=True))
@click.option("-o", "output_file", prompt="the output file", help="The path to the report file", required=True, type=click.Path(exists=False))
def summarize(input_file, taxonpath, names, adiv_file, output_file):
    '''Generate report for participantants'''
    with open(taxonpath) as f:
        taxonpath=json.load(f)

    with open(names) as f:
        names=json.load(f)

    df = pd.read_csv(input_file)
    
    pandas.io.formats.excel.ExcelFormatter.header_style = None
    with pd.ExcelWriter(output_file) as writer:
        wb: xlsxwriter.Workbook = writer.book
        font_fmt = wb.add_format({'font_name': 'Times New Roman', 'font_size': 11})
        # pd.read_csv(adiv_file, sep="\t", skiprows=[1], index_col=['Sample ID']).to_excel(writer, sheet_name='alpha_div')
        adv_df = pd.read_csv(adiv_file, sep='\t', comment='#', index_col='Sample ID')
        adv_df.index.name = None
        adv_df.rename(columns={'shannon': 'Shannon', 'simpson': 'Simpson', 'chao1': 'Chao1', 'observed_features': 'Observed features'}, inplace=True)
        adv_df[['Observed features','Chao1','Shannon','Simpson']].to_excel(writer, sheet_name='alpha_diversity', engine='xlsxwriter')
        ws: xlsxwriter.worksheet.Worksheet = writer.sheets['alpha_diversity']
        ws.set_column(0, len(adv_df.columns), None, font_fmt)
        ws.freeze_panes(1, 1)
        ws.autofilter(0, 0, len(adv_df), len(adv_df.columns))

        all_ranks = ['phylum', 'class', 'order', 'family', 'genus']
        for rank in all_ranks:
            out = df.pivot_table(index='sample_id', columns=f'{rank}_id', values='read_pct', aggfunc='sum', fill_value=0)
            cdic = {t:get_taxonstring(t, taxonpath, names) for t in out.columns.drop('uc')}
            cdic['uc'] = 'unclassified'
            taxid_list = out.columns.to_list()
            out.rename(columns=cdic, inplace=True)

            new_row_df = pd.DataFrame([taxid_list], columns=out.columns, index=['taxonomy_id'])
            out = pd.concat([new_row_df, out])

            clist = out.columns.to_list()
            clist.sort()
            out[clist].to_excel(writer, sheet_name=rank, engine='xlsxwriter')
            ws: xlsxwriter.worksheet.Worksheet = writer.sheets[rank]
            ws.set_column(0, len(out.columns), None, font_fmt)
            # ws.write_string(1, 0, '')
            ws.freeze_panes(2, 1)
            ws.autofilter(1, 0, len(out), len(out.columns))

if __name__ == '__main__':
    create()
