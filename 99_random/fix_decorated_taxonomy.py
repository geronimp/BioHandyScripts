#!/usr/bin/env python
import argparse
import logging
import os
'''
Created on 5 Jun 2015

@author: joelb_000
'''

BAD_TAX=['d__', 'k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

class Correct:
    
    def __init__(self): pass
    
    def getOrigTax(self, taxonomy):
        '''Read in original taxonomy'''
        logging.info("Reading in original taxonomy")
        ancestory={}
        for entry in taxonomy:
            tax_string=[x.replace(' ', '_') for x in entry.rstrip().split('\t')[1].split(';')]

            for idx, entry in enumerate(tax_string):
                if entry in BAD_TAX: continue
                elif entry not in ancestory:
                    if entry.startswith(BAD_TAX[0]): continue
                    ancestory[entry]=tax_string[idx-1]
        return ancestory
    
    def fixDecoTax(self, deco_tax, orig_tax):
        logging.info("Fixing decorated taxonomy")
        new_taxonomy={}
        try:
            for line in [x.rstrip().split('\t') for x in deco_tax]:
                taxonomy=[x.replace(' ', '_') for x in line[1].split('; ')]
                if taxonomy[0] == BAD_TAX[0]:
                    while taxonomy[0] == BAD_TAX[0]:
                        filled_ranks=[x for x in taxonomy if x not in BAD_TAX]
                        if any(filled_ranks):
                            md_taxonomy=filled_ranks[0]
                            taxonomy[taxonomy.index(md_taxonomy)-1]=orig_tax[md_taxonomy]
                        else:
                             logging.warning("%s found to have no taxonomy" % (line[0]))
                             break
                    new_taxonomy[line[0]] = taxonomy
                else:
                    mr_tax           = [x for x in taxonomy if x not in BAD_TAX][-1]
                    mr_tax_ancestors = taxonomy[:taxonomy.index(mr_tax)]
                    if any([x for x in mr_tax_ancestors if x in BAD_TAX]):
                        logging.error('Poorly formed tax sting encountered. You will have to re-root your tree then re-decorate')
                        logging.error('See taxonomy: %s' % line[1])                    
                        raise
                    else:
                        new_taxonomy[line[0]] = taxonomy
        except:
            logging.error("Failed at fixing tax: %s" % ('; '.joing(taxonomy)))
        return new_taxonomy
                    
    def writeTax(self, taxonomy, output):   
        try: 
            logging.info("Writing corrected taxonomy to file: %s" % output)
            with open(output, 'w') as out:
                for key, tax in taxonomy.iteritems():
                    out.write('%s\t%s\n' % (key, '; '.join(tax)))
        except:
            logging.error("Failed writing to file")
            raise
        
    def main(self, orig_tax, deco_tax, out_tax):
        try:
            new_taxonomy = self.fixDecoTax(open(deco_tax), self.getOrigTax(open(orig_tax)))
            if os.path.isdir(out_tax):
                logging.error('File exists: %s' % out_tax)    
                raise
            else:
                self.writeTax(new_taxonomy, out_tax)
            logging.info('Exiting with status: %i' % 0)
            return 0
        except:
            logging.error('Exiting with status: %i' % 1)
            return 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Extraction from GG database using accession ids specified''')
    parser.add_argument('--original_taxonomy', type=str, help='Taxonomy used to decorate', required=True)
    parser.add_argument('--decorated_taxonomy', type=str, help='Decorated taxonomy', required=True)
    parser.add_argument('--output_taxonomy', type=str, help='taxonomy', required=True)
    
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    args = parser.parse_args()
    
    c=Correct()
    exit(c.main(args.original_taxonomy, args.decorated_taxonomy, args.output_taxonomy))

        
        