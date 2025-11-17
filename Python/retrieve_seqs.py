import argparse
import os
import sys
from Bio import Entrez
from Bio import SeqIO

# Enter Entrez email for better results
Entrez.email = "bcmpswork@gmail.com"

def get_org_id(org_name):
    """
    Summary:
        Get the taxonomy ID for a given organism.
    
    Parameters:
        org_name: Name of the organism which taxonomy ID will be fetched (Ex.: 'Homo Sapiens').
    
    Returns:
        record["IdList"][0]: Taxonomy ID of the species.
    """
    try:
        # Search for the organism taxonomy ID in taxonomy database
        handle = Entrez.esearch(db="taxonomy", term=org_name.lower().capitalize())
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"][0]
    except Exception as e:
        print(f"An error occurred while searching for {org_name}: {e}")

def get_taxon_id(org_id, taxon):
    """
    Summary:
        Get the specified taxonomy level ID of the specified organism. 
    
    Parameters:
        org_id: Taxonomic ID of the searched organism.
        taxon: Taxonomic level specified.
    
    Returns: 
        rank["TaxId"]: Taxonomic ID of the taxonomic level specified of the specified organism.
        or
        org_id: Taxonomic ID of the organism if the same belongs to the taxonomic level specified.
    """
    try:
        # Fetch the results for the organism ID obtained
        handle = Entrez.efetch(db="taxonomy", id=org_id)
        record = Entrez.read(handle)
        handle.close()

        # Get the lists conatianing every taxon and taxonomic level of the organism
        lineage = record[0].get("LineageEx", [])
        org_rank = record[0].get("Rank", "Unknown")
        
        # If the specified taxonomic level is the same as the organism it will just retrun the organism taxonomic ID again
        if org_rank.lower() == taxon.lower():
            return org_id
        
        # Returns the match for the taxonomic ID specified
        for rank in lineage:
            if rank["Rank"].lower() == taxon.lower():
                return rank["TaxId"]
        
        print(f"Taxonomic level '{taxon}' not found in lineage.")
        return None
    except Exception as e:
        print(f"An error occurred while searching for {taxon}: {e}")
        sys.exit(1)

def get_genes_list(txid):
    """
    Summary:
        Get a dictionary mapping gene IDs to gene names for a specific organism.
    
    Parameters:
        txid: Taxonomic ID of the specified taxonomic level for the specified organism.
    
    Returns:
        gene_dict: A dictionary containing every gene sequence ID and gene name for the txid.
    """
    try:
        # Creates an empty dictionary to hold gene IDs and names
        gene_dict = {}

        # Performs an intial search to get the number of results to match with retmax
        handle = Entrez.esearch(db="gene", term=f"txid{txid}[Organism]", retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        total_records = int(record["Count"]) 

        # Search for all gene IDs associated with the txid
        handle = Entrez.esearch(db="gene", term=f"txid{txid}[Organism]", retmode="xml", retmax=total_records)
        record = Entrez.read(handle)
        handle.close()

        # Get a list with the gene IDs
        gene_ids = record.get("IdList", [])
        
        # If no genes were found, exit the program
        if not gene_ids:
            print(f"No genes found for taxonomic ID {txid}.")
            sys.exit(1)
        
        # Fetch the information of every gene ID
        print(f"Total genes found: {len(gene_ids)}")
        handle = Entrez.efetch(db="gene", id=",".join(gene_ids), retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        # If the records are in a list format, parse the data
        if isinstance(records, list):
            try:
                # Get the gene ID and gene name
                for gene in records:
                    gene_id = gene["Entrezgene_track-info"]["Gene-track"]["Gene-track_geneid"]
                    gene_name = gene["Entrezgene_gene"]["Gene-ref"]["Gene-ref_locus"]
                    gene_dict[gene_id] = gene_name
            # Give informtion about gene data that couldn't be parsed
            except KeyError as e:
                print(f"Error parsing gene data: {e}")
        
        print(f"Processed {len(gene_dict)} genes for taxonomic ID {txid}.")

        # Return the dictionary with the gene IDs to gene names
        return gene_dict
    
    # If an error occurs during the process, exit the program
    except Exception as e:
        print(f"An error occurred while fetching genes for taxID {txid}: {e}")
        sys.exit(1)

def download_seqs(gene_dict, gene_sequences):
    """
    Summary:
        Fetch nucleotide sequences for each gene in the provided gene dictionary.
    
    Parameters:
        gene_dict: Dictionary containing gene IDs and gene names.
        gene_sequences: Empty dictionary to store gene sequence names and respective gene sequences.
    """
    try:
        # Create a directory to store the downloaded sequences, if it doesn't already exists
        os.makedirs("organism_seqs", exist_ok=True)
        print("Directory created!")
        
        print("Downloading sequences...")
        
        # Loop through each gene in the gene_dict
        for gene_id, gene_name in gene_dict.items():

            # Standardize the gene name to uppercase
            gene_name = gene_name.upper()
            
            # Search for the gene in the nucleotide database using its GeneID
            handle = Entrez.esearch(db="nucleotide", term=f"GeneID:{gene_id}[Gene]", retmode="xml")
            record = Entrez.read(handle)
            handle.close()
            
            # Loop through each sequence ID found for the gene
            for seq_id in record.get("IdList", []):

                # Fetch the sequence record from GenBank using the sequence ID
                handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
                records = list(SeqIO.parse(handle, "genbank"))
                handle.close()
                
                # Process each sequence
                for record in records:

                    # Extract the species name and replace spaces with underscores
                    species_name = record.annotations.get("organism", "Unknown_species").replace(" ", "_")
                    
                    # Create a copy of the record to modify it without altering the original
                    copy = record[:]
                    copy.id = species_name

                    # Delete the description
                    copy.description = ""  
                    
                    # Standardize the sequence to uppercase
                    copy.seq = copy.seq.upper()
                    
                    # Add the gene name to the sequence if it isn't already there
                    if gene_name not in gene_sequences:
                        gene_sequences[gene_name] = []
                    
                    # Append the downloaded gene sequence to the list under the corresponding gene name
                    gene_sequences[gene_name].append(copy)

        print("Download complete!")
    
    except Exception as e:
        # If any error occurs during the downloading process, exit the program
        print(f"Error downloading the sequences: {e}")
        sys.exit(1)


def save_sequences(gene_sequences):
    """
    Summary:
        Save all sequences to FASTA files.
    
    Parameters:
        gene_sequences: Dictionary containing gene sequence names and their respective sequences.
    """
    try:
        # Loop through each gene and its associated sequence records in the gene_sequences dictionary
        for gene_name, records in gene_sequences.items():

            # Open a FASTA file for each gene, named after the gene
            with open(f"organism_seqs/{gene_name}.fasta", "w") as outfile:

                # Write the sequence records to the FASTA file using SeqIO
                SeqIO.write(records, outfile, "fasta")
                print(f"{gene_name}.fasta saved with {len(records)} sequences!")
    
    except Exception as e:
        # If an error occurs while saving the sequences, exit the program
        print(f"Error saving sequences: {e}")
        sys.exit(1)

if __name__ == "__main__":
    # Define input arguments
    parser = argparse.ArgumentParser(description="Fetch gene sequences from NCBI for a given organism and taxonomic level.")
    parser.add_argument("organism", type=str, help="Name of the organism (Ex.: 'Homo sapiens')")
    parser.add_argument("tax_level", type=str, help="Taxonomic level to fetch sequences from (Ex.: 'genus')")
    parser.add_argument("--outgroup", type=str, help="Outgroup organism (optional)", required=False)
    parser.add_argument("--email", type=str, help="Email address (optional)", required=False)
    args = parser.parse_args()

    # Use an email address for better results if given
    if args.email:
        Entrez.email = args.email

    # Define gene dictionary
    gene_sequences = {}
    
    # Get the organism ID
    org_id = get_org_id(args.organism)

    # Get the taxonomic level ID
    txid = get_taxon_id(org_id, args.tax_level)

    # Retrieve sequences
    gene_dict = get_genes_list(txid)

    # Fetch sequences
    download_seqs(gene_dict, gene_sequences)
    
    # Do the same process to the outgroup organism (if specified), but without getting taxonomic level ID
    if args.outgroup:
        out_id = get_org_id(args.outgroup)
        out_gene_dict = get_genes_list(out_id)
        download_seqs(out_gene_dict, gene_sequences)
    
    # Download sequences to FASTA files
    save_sequences(gene_sequences)