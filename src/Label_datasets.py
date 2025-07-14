#!/usr/bin/env python3
"""
Optimized NCBI Taxonomy Retrieval from Accession Number

This script retrieves the full taxonomic lineage for an organism
given its NCBI accession number with minimal server calls.
"""

import requests
import xml.etree.ElementTree as ET
from typing import Optional, Dict, List
import time

class NCBITaxonomyRetriever:
    def __init__(self):
        """Initialize the taxonomy retriever."""
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        
    def get_taxonomy_from_accession(self, accession: str, verbose: bool = False) -> Optional[str]:
        """
        Get full taxonomic lineage from NCBI accession number with minimal API calls.
        
        Args:
            accession: NCBI accession number (e.g., "NM_000546")
            verbose: If True, print debug information
            
        Returns:
            Formatted taxonomy string like "Domain|Kingdom|Phylum|Class|Order|Family|Genus|Species"
            or None if not found
        """
        try:
            # Single optimized call to get both sequence info and taxonomy
            tax_id = self._get_tax_id_optimized(accession, verbose)
            if not tax_id:
                if verbose:
                    print(f"Could not find taxonomy ID for accession: {accession}")
                return None
            
            # Single call to get full taxonomic lineage
            lineage = self._get_taxonomic_lineage(tax_id, verbose)
            if not lineage:
                if verbose:
                    print(f"Could not retrieve taxonomic lineage for tax_id: {tax_id}")
                return None
            
            # Format the lineage
            formatted_lineage = self._format_lineage(lineage)
            return formatted_lineage
            
        except Exception as e:
            if verbose:
                print(f"Error retrieving taxonomy: {e}")
            return None
    
    def _get_tax_id_optimized(self, accession: str, verbose: bool = False) -> Optional[str]:
        """Get taxonomy ID with minimal server calls using esummary."""
        # Try nucleotide database first with esearch + esummary in sequence
        for db in ['nucleotide', 'protein']:
            try:
                # Step 1: Search for the accession
                search_url = f"{self.base_url}esearch.fcgi"
                search_params = {
                    'db': db,
                    'term': accession,
                    'retmode': 'xml',
                    'usehistory': 'y'  # Use history for efficiency
                }
                
                search_response = requests.get(search_url, params=search_params)
                search_response.raise_for_status()
                search_root = ET.fromstring(search_response.content)
                
                # Check if we found the accession
                id_list = search_root.find('.//IdList')
                if id_list is None or len(id_list) == 0:
                    continue
                
                seq_id = id_list.find('Id').text
                if verbose:
                    print(f"Found sequence ID: {seq_id} in {db} database")
                
                # Step 2: Get summary with taxonomy info (more efficient than efetch)
                time.sleep(0.1)  # Rate limiting
                summary_url = f"{self.base_url}esummary.fcgi"
                summary_params = {
                    'db': db,
                    'id': seq_id,
                    'retmode': 'xml'
                }
                
                summary_response = requests.get(summary_url, params=summary_params)
                summary_response.raise_for_status()
                summary_root = ET.fromstring(summary_response.content)
                
                # Look for TaxId in summary (most reliable method)
                for item in summary_root.findall('.//Item'):
                    if item.get('Name') == 'TaxId':
                        tax_id = item.text
                        if verbose:
                            print(f"Found taxonomy ID: {tax_id}")
                        return tax_id
                
                # If TaxId not in summary, try one efetch call as fallback
                time.sleep(0.1)
                fetch_url = f"{self.base_url}efetch.fcgi"
                fetch_params = {
                    'db': db,
                    'id': seq_id,
                    'retmode': 'xml'
                }
                
                fetch_response = requests.get(fetch_url, params=fetch_params)
                fetch_response.raise_for_status()
                fetch_root = ET.fromstring(fetch_response.content)
                
                # Single comprehensive search for taxonomy ID
                tax_id = self._extract_tax_id_from_xml(fetch_root)
                if tax_id:
                    if verbose:
                        print(f"Found taxonomy ID via efetch: {tax_id}")
                    return tax_id
                    
            except Exception as e:
                if verbose:
                    print(f"Error with {db} database: {e}")
                continue
        
        return None
    
    def _extract_tax_id_from_xml(self, root: ET.Element) -> Optional[str]:
        """Extract taxonomy ID from XML using comprehensive patterns."""
        # Define all known patterns for taxonomy ID in NCBI XML
        tax_patterns = [
            './/Org-ref_db/Dbtag[Dbtag_db="taxon"]/Dbtag_tag/Object-id/Object-id_id',
            './/BioSource/BioSource_org/Org-ref/Org-ref_db/Dbtag[Dbtag_db="taxon"]/Dbtag_tag/Object-id/Object-id_id',
            './/Seq-feat/Seq-feat_data/BioSource/BioSource_org/Org-ref/Org-ref_db/Dbtag[Dbtag_db="taxon"]/Dbtag_tag/Object-id/Object-id_id',
            './/Textseq-id_accession/../../..//BioSource/BioSource_org/Org-ref/Org-ref_db/Dbtag[Dbtag_db="taxon"]/Dbtag_tag/Object-id/Object-id_id'
        ]
        
        for pattern in tax_patterns:
            elements = root.findall(pattern)
            if elements:
                return elements[0].text
        
        return None
    
    def _get_taxonomic_lineage(self, tax_id: str, verbose: bool = False) -> Optional[List[Dict]]:
        """Get full taxonomic lineage for a taxonomy ID with single API call."""
        fetch_url = f"{self.base_url}efetch.fcgi"
        params = {
            'db': 'taxonomy',
            'id': tax_id,
            'retmode': 'xml'
        }
        
        time.sleep(0.1)  # Rate limiting
        response = requests.get(fetch_url, params=params)
        response.raise_for_status()
        
        root = ET.fromstring(response.content)
        lineage = []
        
        # Get the main taxon info
        taxon = root.find('.//Taxon')
        if taxon is None:
            return None
        
        # Get scientific name and rank
        sci_name_elem = taxon.find('.//TaxonName[NameClass="scientific name"]/Name')
        rank_elem = taxon.find('Rank')
        
        if sci_name_elem is not None:
            rank = rank_elem.text if rank_elem is not None else 'species'
            lineage.append({
                'rank': rank,
                'name': sci_name_elem.text
            })
        
        # Get complete lineage from parent taxa
        lineage_list = taxon.find('.//LineageEx')
        if lineage_list is not None:
            for taxon_elem in lineage_list.findall('Taxon'):
                rank_elem = taxon_elem.find('Rank')
                name_elem = taxon_elem.find('ScientificName')
                
                if rank_elem is not None and name_elem is not None:
                    lineage.append({
                        'rank': rank_elem.text,
                        'name': name_elem.text
                    })
        
        return lineage
    
    def _format_lineage(self, lineage: List[Dict]) -> str:
        """Format lineage into the requested format."""
        # Define the taxonomic hierarchy we want
        hierarchy = [
            'domain', 'kingdom', 'phylum', 'class', 
            'order', 'family', 'genus', 'species'
        ]
        
        # Create a mapping of ranks to names
        rank_map = {}
        for item in lineage:
            rank = item['rank'].lower()
            # Handle alternative names (NCBI sometimes uses 'superkingdom' instead of 'domain')
            if rank == 'superkingdom':
                rank_map['domain'] = item['name']
            else:
                rank_map[rank] = item['name']
        
        # Build the formatted string
        formatted_parts = []
        for rank in hierarchy:
            if rank in rank_map:
                formatted_parts.append(rank_map[rank])
            else:
                formatted_parts.append('')  # Empty for missing ranks
        
        # Remove trailing empty parts
        while formatted_parts and formatted_parts[-1] == '':
            formatted_parts.pop()
        
        return '|'.join(formatted_parts)
    
    def get_multiple_taxonomies(self, accessions: List[str]) -> Dict[str, Optional[str]]:
        """
        Get taxonomies for multiple accessions efficiently.
        
        Args:
            accessions: List of NCBI accession numbers
            
        Returns:
            Dictionary mapping accession to taxonomy string
        """
        results = {}
        
        for accession in accessions:
            print(f"\nProcessing: {accession}")
            taxonomy = self.get_taxonomy_from_accession(accession)
            results[accession] = taxonomy
            
            # Rate limiting between accessions
            time.sleep(0.2)
        
        return results

def main():
    """Command-line interface for taxonomy retrieval"""
    import sys
    
    if len(sys.argv) != 2:
        print("Usage: python3 script.py <accession>", file=sys.stderr)
        sys.exit(1)
    
    accession = sys.argv[1]
    retriever = NCBITaxonomyRetriever()
    
    taxonomy = retriever.get_taxonomy_from_accession(accession)
    
    if taxonomy:
        print(f">{taxonomy}")
    else:
        print("Could not retrieve taxonomy", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
