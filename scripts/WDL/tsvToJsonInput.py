#!/usr/bin/env python3

import argparse
import json

def convertTsvToWdlInput(tsvFile, workflowName="liftAll"):
    pairs = []
    
    with open(tsvFile) as f:
        for line in f:
            if line.strip():
                genome, annotation = line.strip().split('\t')
                pairs.append({
                    "Left": genome,
                    "Right": annotation
                })
    
    wdlInput = {
        f"{workflowName}.genomeAnnotationPairs": pairs
    }
    
    return json.dumps(wdlInput, indent=2)

def main():
    args = parseArgs()
    jsonInput = convertTsvToWdlInput(args.tsvFile, args.workflowName)
    print(jsonInput)

def parseArgs():
    parser = argparse.ArgumentParser(description='Convert TSV file to WDL input JSON')
    parser.add_argument('tsvFile', 
                       help='Input TSV file with genome and annotation paths')
    parser.add_argument('--workflowName', 
                       default='liftAll',
                       help='Name of the workflow (default: liftAll)')
    return parser.parse_args()

if __name__ == "__main__":
    main()