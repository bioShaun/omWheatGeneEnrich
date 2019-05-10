import click
import os
import delegator
import EnrichConfig
import re
import sys
import pysnooper
import pandas as pd


GENE_PATTERN = re.compile('TraesCS(\w+)1G(\w+)')
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def transfer_id(gene_id):
    """
    TraesCS3D01G355600 -> TraesCS3D02G355600
    """
    if GENE_PATTERN.match(gene_id):
        id_1, id_2 = GENE_PATTERN.match(gene_id).groups()
        return 'TraesCS{id_1}2G{id_2}'.format(**locals())
    else:
        return gene_id


def check_gene_id(gene_list, outdir):
    gene_list_df = pd.read_csv(gene_list, header=None,
                               names=['gene_id'], sep='\t')
    gene_list_df.loc[:, 'gene_id'] = gene_list_df.gene_id.map(transfer_id)
    gene_list_out = os.path.join(outdir, 'gene.list')
    gene_list_df.to_csv(gene_list_out, sep='\t',
                        header=False, index=False)
    return gene_list_df


def genelist_blasttab(species, gene_df, outdir):
    blast_file = os.path.join(
        EnrichConfig.BLAST_DIR,
        'wheat.vs.{}.pep.fasta.blasttab'.format(species))
    blast_df = pd.read_csv(blast_file,
                           header=None,
                           index_col=0, sep='\t')
    overlap_gene = [gene_i for gene_i in gene_df.gene_id
                    if gene_i in blast_df.index]
    gene_blast_df = blast_df.loc[overlap_gene]
    gene_blast_file = os.path.join(outdir, 'gene.blasttab')
    gene_blast_df.to_csv(gene_blast_file, sep='\t', header=None)
    return gene_blast_file


def run_kobas(species, gene_blast_file, outdir):
    enrich_out = os.path.join(outdir, 'enrich.txt')
    kobas_cmd = '/usr/bin/python {run_kobas_py} -i {blast} -t blastout:tab -s {sp} -d K/G -o {out}'.format(
        blast=gene_blast_file, sp=species,
        out=enrich_out, run_kobas_py=EnrichConfig.KOBAS_PY
    )
    delegator.run(kobas_cmd)
    return enrich_out


def format_enrich_file(enrich_out):
    kegg_out_dir, kegg_out_name = os.path.split(enrich_out)
    kegg_tmp_file = os.path.join(kegg_out_dir, 'tmp.%s' % kegg_out_name)

    def check_KOBAS_out(kobas_out):
        kobas_out_info = open(kobas_out, 'r').readlines()
        flag = True
        for eachline in kobas_out_info:
            if not eachline.startswith("#") and len(eachline.strip().split('\t')) == 9:
                flag = True
                break
        else:
            flag = False
        return flag

    if check_KOBAS_out(enrich_out):
        os.system('mv %s %s' % (enrich_out, kegg_tmp_file))
        kegg_out_info = open(enrich_out, 'w')
        with open(kegg_tmp_file, 'r') as kegg_tmp_file_info:
            count = 0
            for eachline in kegg_tmp_file_info:
                if len(eachline.strip().split('\t')) == 9:
                    if count == 0 and eachline.startswith("#"):
                        kegg_out_info.write(eachline)
                        count += 1
                    elif not eachline.startswith("#"):
                        kegg_out_info.write(eachline)
        kegg_out_info.close()
        os.system('rm %s' % (kegg_tmp_file))


def full_term_id(term_id, sp):
    if term_id.startswith('GO'):
        return term_id
    else:
        return '{sp}{term}'.format(
            term=term_id, sp=sp
        )


def add_term_cat(enrich_file, sp):
    enrich_df = pd.read_csv(enrich_file, sep='\t')
    term_df = pd.read_csv(EnrichConfig.TERM_CAT, sep='\t')
    term_df.loc[:, 'ID'] = [full_term_id(term, sp) for
                            term in term_df.ID]
    enrich_df = enrich_df.merge(term_df)
    enrich_df.to_csv(enrich_file, sep='\t', index=False)


def species_abbr(species):
    name_file = os.path.join(SCRIPT_DIR, 'wheatdb.organism.avail.txt')
    name_df = pd.read_csv(name_file, sep='\t', index_col=2,
                          header=None, names=['id', 'abbr', 'kingdom'])
    if species in name_df.abbr.values:
        return species
    else:
        if species in name_df.index:
            return name_df.loc[species].abbr
        else:
            sys.exit('species [{}] not found!'.format(species))


@click.command()
@click.option(
    '-s',
    '--species',
    type=click.STRING,
    help='KEGG species abbr or full name [hsa/Homo sapiens (human)]',
    required=True
)
@click.option(
    '-g',
    '--gene_list',
    help='gene list file.',
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '-o',
    '--outdir',
    help='output directory.',
    type=click.Path(),
    required=True
)
def main(species, gene_list, outdir):
    species = species_abbr(species)
    gene_list = os.path.abspath(gene_list)
    outdir = os.path.abspath(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    gene_list_df = check_gene_id(
        gene_list, outdir)
    gene_blast_file = genelist_blasttab(
        species, gene_list_df, outdir)
    enrich_out = run_kobas(
        species, gene_blast_file, outdir)
    format_enrich_file(enrich_out)
    add_term_cat(enrich_out, species)


if __name__ == '__main__':
    main()
