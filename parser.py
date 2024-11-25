# Description: a script that parses and converts to JSON an annotation format (designed by Professor Andrew Martin) for antibody-based drugs


import sys
import os
import re
import json


# regular expression patterns to match specific records or groups of records
# TODO: move to another file and import?
p_request = re.compile(r'^Request:\s(\d+);?\n*') # at least one example (12082.txt) is missing a ';' delimiter, so the request record needs to be extracted separately
p_note = re.compile(r'(.+?):\s(.+?);\n*?(Note):\s(.+?);')
p_chain = re.compile(r'((?:(?:Heavy|Light)\s)?Chain)\[([0-9-,]+)\]:\n*?([\s\S]+?)')


# this function carries out the details of the parsing
def txt_to_json(txtpath, outpath):
    
    # initialise dictionary where the parsed records will go
    annotations = {}
    
    with open(txtpath, 'r') as f:
        txt = f.read()

    request = p_request.search(txt)
    if request:
        annotations['Request'] = request.group(1) # only the digits of the matching pattern
    else:
        annotations['Request'] = os.path.basename(txtpath).rstrip('.txt')
        annotations['Request-Note'] = 'from file name'

    # remove the Request record from the txt string
    txt = p_request.split(txt)[-1]

    # append the related record to each of the Note keywords
    # using a custom replace function with re.sub()
    def noterepl(match):
        related_key = match.group(1)
        related_value = match.group(2)
        note_value = match.group(4)
        return f'{related_key}: {related_value};\n{related_key}-Note: {note_value};\n'
    txt = p_note.sub(noterepl, txt)
    
    # create temp list using ';' as record delimiter and removing new line characters if they exist
    records = txt.split(';')
    records = [record.strip('\n') for record in records]

    # extract and remove the chain records as different parsing is needed
    chains = records.pop()

    # populate the dictionary with the records (except chains)
    for record in records:
        key, value = record.split(':', maxsplit=1)
        annotations[key] = value

    # parse the chain records
    chains = chains.split('//')[-1:]
    chains = [chain.strip('\n') for chain in chains]
    for chain in chains:
        match = p_chain.search(chain)
        chain_key = match.group(1).replace(' ', '')
        chain_instance = match.group(2)
        chain_seq = match.group(3)
        chain_seq = re.sub(r'[^A-Z]', '', chain_seq)
        # when some instances have identical chain sequences
        if ',' in chain_instance:
            instances = [int(i) for i in chain_instance]
            n_instances = max(instances)
            # so that array indices correspond to the instance number
            chain_seqs = [None] * n_instances
            for i in instances:
                chain_seqs[i] = chain_seq
            annotations[chain_key] = chain_seqs
        else:
            chain_key = chain_key + chain_instance
            annotations[chain_key] = chain_seq

    with open(outpath, 'w') as out:
        json.dump(annotations, out)


# this function gets command-line argurments, sets the output path, and calls txt_to_json()
def run_parser():
    
    outdir = os.getcwd() + '/json/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for fpath in sys.argv[1:]:
        inpath = os.path.abspath(fpath)
        outpath = outdir + os.path.basename(inpath).rstrip('.txt') + '.json'
        txt_to_json(inpath, outpath)


if __name__ == "__main__":   
    run_parser()


# TODO: flags/options for different types of user input: file path (one or many) vs directory path
# TODO: split some of the records into two e.g. Antigen-Gene for the gene name in parentheses, and CDRKabatH1-Range for the residue range in parentheses
