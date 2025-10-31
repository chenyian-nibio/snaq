import pandas as pd
import pandas.io.formats.excel
import json
import click
import xlsxwriter
import xlsxwriter.worksheet

@click.group
def create():
    pass

@create.command
@click.option("-i", "input_file", prompt="the input file", help="The csv file for manta", required=True, type=click.Path(exists=True))
@click.option("-a", "adiv_file", prompt="the alpha diversity file", help="The tsv file for alpha diversity", required=True, type=click.Path(exists=True))
@click.option("-t", "taxonomy", prompt="the manta taxonomy file", help="The path to the taxonomy", required=True, type=click.Path(exists=True))
@click.option("-o", "output_file", prompt="the output file", help="The path to the report file", required=True, type=click.Path(exists=False))
def summarize(input_file, taxonomy,adiv_file, output_file):
    '''Generate report for participantants'''
    df = pd.read_csv(input_file)

    taxon_df = pd.read_csv(taxonomy)
    names = dict(zip(taxon_df.iloc[:, 0], taxon_df.iloc[:, 2]))
    
    pandas.io.formats.excel.ExcelFormatter.header_style = None
    with pd.ExcelWriter(output_file) as writer:
        wb: xlsxwriter.Workbook = writer.book
        font_fmt = wb.add_format({'font_name': 'Times New Roman', 'font_size': 11})
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

            if rank == 'phylum':
                df['name']=f'{rank[0]}_' + df[f'{rank}_id'].apply(lambda x: names[str(x)])
            else:
                df['name']=df['name'] + f';{rank[0]}_' + df[f'{rank}_id'].apply(lambda x: names[str(x)])
            
            out = df.pivot_table(index='sample_id', columns='name', values='read_pct', aggfunc='sum', fill_value=0)

            taxid_dic = df[[f'{rank}_id', 'name']].drop_duplicates().set_index('name').to_dict()[f'{rank}_id']
            taxid_list = []
            for tname in out.columns:
                id = taxid_dic[tname]
                taxid_list.append('-' if id == 'uc' else id)

            new_row_df = pd.DataFrame([taxid_list], columns=out.columns, index=['taxonomy_id'])
            out = pd.concat([new_row_df, out])

            clist = out.columns.to_list()
            clist.sort()
            out[clist].to_excel(writer, sheet_name=rank, engine='xlsxwriter')
            ws: xlsxwriter.worksheet.Worksheet = writer.sheets[rank]
            ws.set_column(0, len(out.columns), None, font_fmt)
            ws.freeze_panes(2, 1)
            ws.autofilter(1, 0, len(out), len(out.columns))

if __name__ == '__main__':
    create()
