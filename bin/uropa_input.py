#!/usr/bin/env python
import argparse
import json

def main(args):
    base_query = {
        "feature": "gene",
        "filter.attribute": "gene_type",
        "attribute.value": "protein_coding", 
        "feature.anchor": "start" 
    }

    for i, peak_type in enumerate(args.peak_types):
        json_construct = dict()
        json_construct['queries'] = []
        json_construct['show_attributes'] = ["gene_id", "gene_name", "gene_type"]
        json_construct["priority"] = "Yes"
        json_construct['gtf'] = args.gtf
        json_construct['bed'] = args.bed

        if args.assay == 'cfchip':
            if peak_type == 'protTSS':
                for _d in (3000, 10000, 100000):
                    this_q = base_query.copy()
                    this_q['distance'] = _d
                    json_construct['queries'].append(this_q)
        else:
            if peak_type == 'prot':
                for _d in (5000, 100000):
                    this_q = base_query.copy()
                    del this_q["feature.anchor"]
                    this_q['distance'] = _d
                    json_construct['queries'].append(this_q)
            elif peak_type == 'genes':
                this_query = {}
                this_query['feature'] = 'gene'
                for _d in (5000, 100000):
                    this_q = base_query.copy()
                    del this_q["feature.anchor"]
                    del this_q["filter.attribute"]
                    del this_q["attribute.value"]
                    this_q['distance'] = _d
                    json_construct['queries'].append(this_q)
            elif peak_type == 'protSEC':
                # distance, feature.anchor
                query_values = (
                    ([3000, 1000], "start"), 
                    (3000,         "end"), 
                    (100000,       "center"), 
                    (100000,       None)
                )
                for _distance, feature_anchor in query_values:
                    this_q = base_query.copy()
                    del this_q["feature.anchor"]
                    if feature_anchor: 
                        this_q["feature.anchor"] = feature_anchor
                    this_q['distance'] = _distance
                    json_construct['queries'].append(this_q)
            elif peak_type == 'protTSS':
                for _d in ([3000, 1000], 10000, 100000):
                    this_q = base_query.copy()
                    this_q['distance'] = _d
                    json_construct['queries'].append(this_q)

        with open(args.output_json[i], 'w') as jo:
            json.dump(json_construct, jo, indent=4)
            jo.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to prepare the uropa json input')
    parser.add_argument('-g', dest='gtf', required=True, help='')
    parser.add_argument('-o', dest='output_json', nargs="+", required=True, help='')
    parser.add_argument('-a', dest='assay', required=True, help='')
    parser.add_argument('-b', dest='bed', required=True, help='')
    parser.add_argument('--types', '-t', dest='peak_types', nargs="+", required=True, help='')
    args = parser.parse_args()