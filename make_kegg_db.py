import click
import os
import sys
import delegator
import EnrichConfig


def place_wheat_pep(wheat_pep):
    cp_flag = 0
    if os.path.exists(EnrichConfig.WHEAT_PEP_DB):
        if wheat_pep is not None:
            if click.confirm(
                    '{} file exists, overwrite?'.format(EnrichConfig.WHEAT_PEP_DB)):
                cp_flag = 1
    else:
        cp_flag = 1
        if wheat_pep is None:
            sys.exit(
                '{} file not exist. --wheat_pep is required!'.format(
                    EnrichConfig.WHEAT_PEP_DB))
    if cp_flag:
        os.system('cp {wheat_pep} {wheat_pep_db}'.format(
            wheat_pep=wheat_pep,
            wheat_pep_db=EnrichConfig.WHEAT_PEP_DB
        ))


def download_file(url, outfile):
    dl_cmd = 'wget -O {out}.gz "{url}"'.format(
        out=outfile, url=url
    )
    uncomp_cmd = 'gunzip -f {out}.gz'.format(
        out=outfile
    )
    delegator.run(dl_cmd)
    delegator.run(uncomp_cmd)


def makeblastdb(pepfile):
    cmd = 'makeblastdb -in {pep} -dbtype prot'.format(
        pep=pepfile
    )
    delegator.run(cmd)


def blast2wheat(pepfile, threads):
    pep_name = os.path.basename(pepfile)
    blastout = os.path.join(EnrichConfig.BLAST_DIR,
                            'wheat.vs.{}.blasttab'.format(pep_name))
    if os.path.exists(blastout):
        if not click.confirm(
            '{} file exists, overwrite?'.format(blastout)
        ):
            return 1
    blast_cmd = 'blastp -query {wheat} -db {db} -num_threads {threads} -out {out} {default}'.format(
        wheat=EnrichConfig.WHEAT_PEP_DB,
        db=pepfile, threads=threads,
        default=EnrichConfig.BLAST_DEFAULT,
        out=blastout
    )
    delegator.run(blast_cmd)


def prepare_sp_files(species, threads):
    organism_db = os.path.join(EnrichConfig.DB_DIR,
                               'organism.db')
    if not os.path.exists(organism_db):
        organism_db_url = EnrichConfig.DB_URL_TEMP.format('organism')
        download_file(organism_db_url, organism_db)
    pep_file = os.path.join(EnrichConfig.PEP_DIR,
                            '{}.pep.fasta'.format(species))
    db_file = os.path.join(EnrichConfig.DB_DIR,
                           '{}.db'.format(species))
    dl_flat = 1
    if os.path.exists(pep_file):
        if not click.confirm(
                '{} file exists, overwrite?'.format(pep_file)):
            dl_flat = 0
    if dl_flat:
        pep_url = EnrichConfig.PEP_URL_TEMP.format(species)
        db_url = EnrichConfig.DB_URL_TEMP.format(species)
        download_file(pep_url, pep_file)
        download_file(db_url, db_file)
        makeblastdb(pep_file)
    blast2wheat(pep_file, threads)


@click.command()
@click.option(
    '-s',
    '--species',
    help='KEGG species abbrs, sep with comma.',
    required=True
)
@click.option(
    '-w',
    '--wheat_pep',
    help='wheat pep file.',
    default=None
)
@click.option(
    '-t',
    '--threads',
    help='blast threads.',
    default=16
)
def main(species, wheat_pep, threads):
    place_wheat_pep(wheat_pep)
    species_list = species.split(',')
    [prepare_sp_files(sp_i, threads)
     for sp_i in species_list]


if __name__ == '__main__':
    main()
