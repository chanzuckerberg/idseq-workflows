import argparse
import json
import logging
import time
import xml.etree.ElementTree as ET
from subprocess import run, PIPE
from typing import TypedDict, Iterable


class Metadata(TypedDict):
    name: str
    country: str
    collection_date: str


def fetch_ncbi(accession):
    query = accession
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    search_url = f"{base}/esearch.fcgi?db=nuccore&term={query}&usehistory=y"
    output = run(["curl", search_url], stdout=PIPE, check=True).stdout
    time.sleep(1)
    root = ET.fromstring(output)
    web = root.find('WebEnv').text
    key = root.find('QueryKey').text
    fetch_url = f"{base}/efetch.fcgi?db=nuccore&query_key={key}&WebEnv={web}&rettype=gb&retmode=xml"
    genbank_xml = run(["curl", fetch_url], stdout=PIPE, check=True).stdout
    return {
        'search_url': search_url,
        'fetch_url': fetch_url,
        'genbank_xml': genbank_xml
    }


def get_accession_metadata(accession):
    '''
    Retrieve metadata of an NCBI accession (e.g. name, country, collection date)
    TODO: Put this data in S3 instead and get it from there.
    '''
    accession_metadata = {}
    fetch_ncbi_result = fetch_ncbi(accession)
    genbank_xml = fetch_ncbi_result['genbank_xml']

    root = ET.fromstring(genbank_xml).find('GBSeq')
    time.sleep(1)
    if not root:
        logging.warn(f"{fetch_ncbi_result} did not return a result")
        return accession_metadata
    accession_metadata['name'] = root.find('GBSeq_definition').text
    qualifiers_needed = {'country', 'collection_date'}
    for entry in root.find('GBSeq_feature-table')[0].find('GBFeature_quals'):
        if all(key in accession_metadata for key in qualifiers_needed):
            break
        for key in qualifiers_needed - accession_metadata.keys():
            if entry.find('GBQualifier_name').text == key:
                accession_metadata[key] = entry.find('GBQualifier_value').text
    return accession_metadata

def main(reference_accession_ids: Iterable[str], output_ncbi_metadata: str):
    with open(output_ncbi_metadata, "w") as f:
        json.dump({ a: get_accession_metadata(a) for a in reference_accession_ids }, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-accession-ids")
    parser.add_argument("--output-ncbi-metadata")
    args = parser.parse_args()

    with open(args.reference_accession_ids) as f:
        reference_accession_ids: Iterable[str] = json.load(f)

    main(reference_accession_ids, args.output_ncbi_metadata)
