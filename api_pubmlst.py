import requests
import base64
import gzip
import bz2
from pathlib import Path
import pandas as pd
from multiprocessing import Pool

magic_dict = {
    b"\x1f\x8b\x08": (gzip.open, 'rb'),
    b"\x42\x5a\x68": (bz2.BZ2File, 'r'),
    }

max_len = max(len(x) for x in magic_dict)

def open_by_magic(filename):
    with open(filename,  "rb") as f:
        file_start = f.read(max_len)
    for magic, (fn, flag) in magic_dict.items():
        if file_start.startswith(magic):
            return fn(filename, flag)
    return open(filename, 'r')

equiv = {key.split('\t')[0]: (key.split('\t')[1], key.split('\t')[2].strip()) for key in open('fastmlst_pubmlst_schnum.tsv').readlines()}

def get_st_pubmlst(assembly, scheme):
    fasta = open_by_magic(assembly).read()

    data = {
            'base64': 'true',
            'sequence': base64.b64encode(fasta).decode()
        }
    pasteur = ['kingella', 'staphlugdunensis', 'listeria']
    if equiv[scheme][0] in pasteur:
        endpoint = f'https://bigsdb.pasteur.fr/api/db/pubmlst_{equiv[scheme][0]}_seqdef/schemes/{equiv[scheme][1]}/sequence'
    else:
        endpoint = f'https://rest.pubmlst.org/db/pubmlst_{equiv[scheme][0]}_seqdef/schemes/{equiv[scheme][1]}/sequence'
    r = requests.post(endpoint, json=data)
    if 'fields' in r.json().keys(): # Found a ST, make a tuple (Genome, Scheme, ST, {alleles})
        rjson = r.json()
        response = {
            'Genome':assembly.name,
            'Scheme': scheme,
            'ST': rjson['fields']['ST'],
        }
        for allel, var in rjson['exact_matches'].items():
            response[allel] = var[0]['allele_id']
        return(response)
    elif 'exact_matches' in r.json().keys():
        rjson = r.json()
        response = {
            'Genome':assembly.name,
            'Scheme': scheme,
            'ST': '-',
        }
        for allel, var in rjson['exact_matches'].items():
            response[allel] = var[0]['allele_id']
        return(response)
    else:
        print('wtf:',assembly , scheme)
        print(endpoint)
        print(r)


def process_dir(dirname):
    directory = Path(str(dirname))
    output = Path('pubmlst_output')
    fastas = list(directory.glob('*.gz'))
    schemename = dirname.name
    outfile = f'{schemename}.csv'
    datadic = []
    i = 1
    if not (output/outfile).is_file():
        for fasta in fastas:
            datadic.append(get_st_pubmlst(fasta, schemename))
            print(f'{schemename}: {i} de {len(fastas)}')
            i += 1
        dfout = pd.DataFrame(datadic)
        dfout.to_csv(output/outfile, sep=',', index=False)
        return pd.DataFrame(datadic)
    else:
        print(outfile)
        return True


if __name__ == "__main__":
    genomes = Path('genomes')
    output = Path('pubmlst_output')
    species = list(genomes.glob('*'))
#    for specie in genomes.glob('*'):
#        print(specie.name)
#        if not (output/specie.name).is_file():
#            process_dir(specie)
    with Pool(12) as p:
        p.map(process_dir, species)
