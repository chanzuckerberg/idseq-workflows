import argparse
import json
from typing import TypedDict, Iterable


class Metadata(TypedDict):
    name: str

def main(reference_accession_ids: Iterable[str], output_ncbi_metadata: str):
    with open(output_ncbi_metadata, "w") as f:
        json.dump({ s: s for s in reference_accession_ids }, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-accession-ids")
    parser.add_argument("--output-ncbi-metadata")
    args = parser.parse_args()

    with open(args.reference_accession_ids) as f:
        reference_accession_ids: Iterable[str] = json.load(f)

    main(reference_accession_ids, args.output_ncbi_metadata)
