import pandas as pd
import click

@click.command()
@click.option("-a", "adv_file", required=True, type=str)
@click.option("-d", "biom_file", required=True, type=str)
@click.option("-o", "report_file", required=True, type=str)
def prepare_return_report(adv_file, biom_file, report_file):

    with pd.ExcelWriter(report_file, mode='w') as writer:
        adv_df = pd.read_csv(adv_file, sep='\t', comment='#', index_col='Sample ID')
        adv_df.index.name = None
        adv_df.rename(columns={'shannon_entropy': 'Shannon', 'simpson': 'Simpson', 'chao1': 'Chao1', 'observed_features': 'Observed features'}, inplace=True)
        adv_df.to_excel(writer, sheet_name='alpha_diversity')

        df = pd.read_csv(biom_file, sep='\t', skiprows=1)
        df = df[df['#OTU ID'].str.startswith('d__Bacteria;p__')]
        ranks = ['phylum', 'class', 'order', 'family', 'genus']
        for i, r in enumerate(ranks):
            df[r] = df['#OTU ID'].apply(lambda x: ';'.join(x.split(';')[1:(i + 2)]))
            df[~df[r].str.endswith(';__')].drop(['#OTU ID'], axis=1).groupby([r]).sum().T.to_excel(writer, sheet_name=f'{r}')
            df.drop([r], axis=1, inplace=True)

if __name__ == "__main__":
    prepare_return_report()

