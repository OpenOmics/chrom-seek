#!/usr/bin/env python
import os
import argparse
import json

def main(args):
    base_query = {
        "feature": ["gene"],
        "filter_attribute": "gene_type",
        "attribute_values": ["protein_coding"], 
        "feature_anchor": ["start"],
        "relative_location": 
            ["PeakInsideFeature", "FeatureInsidePeak", "Upstream",
             "Downstream", "OverlapStart", "OverlapEnd"],
        "strand": "ignore"
    }

    for i, peak_type in enumerate(args.peak_types):
        if not os.path.exists(os.path.dirname(args.output_json[i])):
            os.makedirs(os.path.abspath(os.path.dirname(args.output_json[i])))
        json_construct = dict()
        json_construct['queries'] = []
        json_construct['show_attributes'] = ["gene_id", "gene_name", "gene_type"]
        json_construct["priority"] = "Yes"
        json_construct["outdir"] = os.path.dirname(args.output_json[i])
        json_construct['gtf'] = args.gtf
        json_construct['bed'] = args.bed

        if args.assay == 'cfchip':
            if peak_type == 'protTSS':
                for ii, _d in enumerate([[3000], [10000], [100000]], start=1):
                    this_q = base_query.copy()
                    this_q['distance'] = _d
                    this_q['name'] = f'query_{str(ii)}'
                    json_construct['queries'].append(this_q)
        else:
            if peak_type == 'prot':
                for ii, _d in enumerate([[5000], [100000]], start=1):
                    this_q = base_query.copy()
                    del this_q["feature_anchor"]
                    this_q['distance'] = _d
                    this_q['name'] = f'query_{str(ii)}'
                    json_construct['queries'].append(this_q)
            elif peak_type == 'genes':
                this_query = {}
                this_query['feature'] = 'gene'
                for ii, _d in enumerate([[5000], [100000]], start=1):
                    this_q = base_query.copy()
                    del this_q["feature_anchor"]
                    del this_q["filter_attribute"]
                    del this_q["attribute_value"]
                    this_q['distance'] = _d
                    this_q['name'] = f'query_{str(ii)}'
                    json_construct['queries'].append(this_q)
            elif peak_type == 'protSEC':
                query_values = (
                    ([3000, 1000], ["start"]), 
                    ([3000], ["end"]), 
                    ([100000], ["center"]),
                    ([100000], None)
                )
                for ii, (_distance, feature_anchor) in enumerate(query_values, start=1):
                    this_q = base_query.copy()
                    del this_q["feature_anchor"]
                    if feature_anchor:
                        this_q["feature_anchor"] = feature_anchor
                    this_q['distance'] = _distance
                    this_q['name'] = f'query_{str(ii)}'
                    json_construct['queries'].append(this_q)
            elif peak_type == 'protTSS':
                for ii, _d in enumerate([[3000, 1000], [10000], [100000]], start=1):
                    this_q = base_query.copy()
                    this_q['distance'] = _d
                    this_q['name'] = f'query_{str(ii)}'
                    json_construct['queries'].append(this_q)

        with open(args.output_json[i], 'w') as jo:
            json.dump(json_construct, jo, indent=4)
            jo.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to prepare the uropa json input')
    parser.add_argument('-g', dest='gtf', required=True, help='Gene GTF used in uropa annotation lookup')
    parser.add_argument('-o', dest='output_json', nargs="*", required=True, help='Path to output UROPA input JSON')
    parser.add_argument('-a', dest='assay', required=True, help='Type of assay being run')
    parser.add_argument('-b', dest='bed', required=True, help='Bed used for UROPA annotation')
    parser.add_argument('--types', '-t', dest='peak_types', nargs="+", required=True, help='Peak types: prot, protTSS, genes, protSEC')
    args = parser.parse_args()
    if isinstance(args.output_json, str):
        args.output_json = [args.output_json]
    if len(args.peak_types) != len(args.output_json):
        raise ValueError('More peaks types than output file paths!')
    main(args)