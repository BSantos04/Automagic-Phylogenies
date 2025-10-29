// Define parameters
params.o = ""
params.t = ""
params.outgroup = ""
params.email = ""
params.docker = ""

// Process to download the sequences using EntrezAPI
process dowloadSeqs {

    container "${params.docker}"
    // Define input parameters
    input:
    val(organism)
    val(taxonomy)

    // Define expected output directory
    output:
    path "organism_seqs/", type: "dir"

    // Copies output to a customized directory
    def outdir = params.outgroup ? "${params.o.replace(' ', '_')}_${params.t}_${params.outgroup.replace(' ', '_')}_Results" : "${params.o.replace(' ', '_')}_${params.t}_Results"
    publishDir outdir, mode : "copy"

    // Executes the python script to download the needed data
    script:
    def outgroup = params.outgroup ? "--outgroup '${params.outgroup}'" : ""
    def email = params.email ? "--email '${params.email}'" : ""
    """
    python3 /auto-phylo/Python/retrieve_seqs.py '${organism}' '${taxonomy}' ${outgroup} ${email}
    """
}

// Process to check FASTA files for one sequence and delete if so
process deleteSingleSeqs {

    container "${params.docker}"
    // Define input parameters
    input:
    path fasta_dir

    // Overwrites output directory
    output:
    path fasta_dir, emit: cleaned_fasta_dir

    // Copies output to a customized directory
    def outdir = params.outgroup ? "${params.o.replace(' ', '_')}_${params.t}_${params.outgroup.replace(' ', '_')}_Results" : "${params.o.replace(' ', '_')}_${params.t}_Results"
    publishDir outdir, mode : "copy"

    // Run command-lines to check if the files only contains one sequence
    script:
    """
    # Find and delete any FASTA files with only one sequence
    for file in ${fasta_dir}/*.fasta; do
        seq_count=\$(grep -c "^>" "\$file")
        if [ \$seq_count -eq 1 ]; then
            echo "File \$file contains only one sequence, deleting it."
            rm "\$file"
        fi
    done

    # Check if any FASTA files remain in the directory
    remaining_files=\$(ls ${fasta_dir}/*.fasta 2>/dev/null | wc -l)

    if [ \$remaining_files -eq 0 ]; then
        echo "No FASTA files remaining in the directory. Exiting the pipeline."
        exit 1
    fi

    # If files remain, just pass the cleaned directory
    echo "Remaining FASTA files found, continuing with the pipeline."
    """
}

// Process to align FASTA sequences using KAlign2
process alignSeqs {

    container "${params.docker}"
    // Define input parameters
    input:
    path fasta

    // Define expected output directory
    output:
    path "aln_data/", type: "dir"

    // Copies output to a customized directory
    def outdir = params.outgroup ? "${params.o.replace(' ', '_')}_${params.t}_${params.outgroup.replace(' ', '_')}_Results" : "${params.o.replace(' ', '_')}_${params.t}_Results"
    publishDir outdir, mode: "copy"

    // Align every sequence in the folder previously created and saves the alignments in another folder
    script:
    """
    # Create output directory if it doesn't exist
    mkdir -p "aln_data/"

    # Loop through all FASTA files in the input directory
    for file in ${fasta}/*.fasta; do
        # Extract the filename without path and extension
        filename=\$(basename "\$file" .fasta)

        # Run KAlign2 and save the aligned file in the output directory
        kalign -i "\$file" -o "aln_data/\${filename}_aligned.fasta" -f 0
    done
    """
}

// Process to concatenate aligned sequences using AMAS
process concatenateSeqs {

    container "${params.docker}"
    // Define input parameters
    input:
    path aln_dir

    // Define expected output directory
    output:
    path "concatenated_seqs.fasta"

    // Copies output to a customized directory
    def outdir = params.outgroup ? "${params.o.replace(' ', '_')}_${params.t}_${params.outgroup.replace(' ', '_')}_Results" : "${params.o.replace(' ', '_')}_${params.t}_Results"
    publishDir outdir, mode: "copy"

    // Run AMAS command-line to concatenate every sequences into one
    script:
    """
    AMAS.py concat -i ${aln_dir}/*.fasta -f fasta -d dna -u fasta -t concatenated_seqs.fasta
    """
}

// Process to find the best model for the aligned sequences to build a ML tree
process findBestModel {

    container "${params.docker}"
    // Define input parameter
    input:
    path(concat_file)

    // Define expected output file
    output:
    path "model.out"

    // Copies output to current directory
    def outdir = params.outgroup ? "${params.o.replace(' ', '_')}_${params.t}_${params.outgroup.replace(' ', '_')}_Results" : "${params.o.replace(' ', '_')}_${params.t}_Results"
    publishDir outdir, mode: "copy"

    // Executes Modeltest-NG command-line to get the best model
    script:
    """
    modeltest-ng -d nt -i ${concat_file} -o model -T raxml
    """
}

// Process to run RAxML-NG to get the best ML tree
process bestTree {

    container "${params.docker}"
    // Define input parameters
    input:
    path(model_out)
    path(concat)

    // Define expected output file inside a tuple
    output:
    path("output.raxml.bestTree")

    // Copies output to current directory
    def outdir = params.outgroup ? "${params.o.replace(' ', '_')}_${params.t}_${params.outgroup.replace(' ', '_')}_Results" : "${params.o.replace(' ', '_')}_${params.t}_Results"
    publishDir outdir, mode: "copy"

    // Extracts the best model for RAxML-NG (the last value registered) and execute the command-line to build the best tree
    script:
    def outgroup = params.outgroup ? "--outgroup \$(grep '^>${params.outgroup.replace(' ', '_')}' ${concat} | sed 's/^>//' | tr '\n' ' ' | xargs)" : ""
    """
    # Extract best model for RAxML-NG
    best_model=\$(awk '/> raxml-ng --msa/ {print \$NF}' ${model_out} | tail -n 1)

    # Get total number of CPUs
    cpu=\$(nproc)

    # Run RAxML-NG command-line using the extracted model and aligned file
    raxml-ng --msa ${concat} --model \${best_model} --threads \${cpu} --prefix output --force perf_threads ${outgroup}
    """
}


// Run RAxML-NG command-line to get bootstrap trees
process bootstrapTrees {

    container "${params.docker}"
    // Define input parameters
    input:
    path(model_out)
    path(concat)

    // Define expected output file inside a tuple
    output:
    path("output.raxml.bootstraps")

    // Copies output to current directory
    def outdir = params.outgroup ? "${params.o.replace(' ', '_')}_${params.t}_${params.outgroup.replace(' ', '_')}_Results" : "${params.o.replace(' ', '_')}_${params.t}_Results"
    publishDir outdir, mode: "copy"

    // Extracts the best model for RAxML-NG (the last value registered) and execute the command-line to build the bootstrap trees
    script:
    def outgroup = params.outgroup ? "--outgroup \$(grep '^>${params.outgroup.replace(' ', '_')}' ${concat} | sed 's/^>//' | tr '\n' ' ' | xargs)" : ""
    """
    # Extract best model for RAxML-NG
    best_model=\$(awk '/> raxml-ng --msa/ {print \$NF}' ${model_out} | tail -n 1)

    # Get total number of CPUs
    cpu=\$(nproc)

    # Run RAxML-NG command-line using the extracted model and aligned file
    raxml-ng --bootstrap --msa ${concat} --model \${best_model} --threads \${cpu} --prefix output --bs-trees autoMRE --force perf_threads ${outgroup}
    """
}

// Run RAxML-NG command-line to get support tree
process supportTree {

    container "${params.docker}"
    // Define input parameters
    input:
    path(best_tree)
    path(bs_trees)

    // Define expected output file inside a tuple
    output:
    path("output.raxml.bestTree.raxml.support")

    // Copies output to current directory
    def outdir = params.outgroup ? "${params.o.replace(' ', '_')}_${params.t}_${params.outgroup.replace(' ', '_')}_Results" : "${params.o.replace(' ', '_')}_${params.t}_Results"
    publishDir outdir, mode: "copy"

    // Execute RAxML-NG command-line to create support tree using best tree and bootstrap trees
    script:
    """
    # Get total number of CPUs
    cpu=\$(nproc)

    # Run RAxML-NG command-line using the extracted model and aligned file
    raxml-ng --support --tree ${best_tree} --bs-trees ${bs_trees} --threads \${cpu} --force perf_threads 
    """
}

// Convert the tree into an SVG file using Toytree
process buildTree {

    container "${params.docker}"
    // Define input parameters
    input:
    path(support_tree)

    // Define expected output file
    output:
    path "phylo_tree.svg"

    // Copies output to current directory
    def outdir = params.outgroup ? "${params.o.replace(' ', '_')}_${params.t}_${params.outgroup.replace(' ', '_')}_Results" : "${params.o.replace(' ', '_')}_${params.t}_Results"
    publishDir outdir, mode: "copy"

    // Build phylogenetic tree using the support file provided by RAxML-NG, converting it into a SVG file
    script:
    """
    python3 /auto-phylo/Python/tree_generator.py ${support_tree}
    """
}

workflow {
    // Raise an error message for missing parameters
    if (!params.o || !params.t || !params.docker) {
        error "Error: No input file provided. Use --o '{organism name}' --t '{taxonomy level} --docker {Docker image}'"
    }

    // Convert inputs into channels
    o_ch = channel.of(params.o)
    t_ch = channel.of(params.t)

    // Download sequences
    seqs_dir = dowloadSeqs(o_ch, t_ch)

    // Delete files with one sequence only
    clean_dir = deleteSingleSeqs(seqs_dir)

    // Align sequences
    aln_dir = alignSeqs(clean_dir)

    // Concatenate sequences using AMAS
    concat_seqs = concatenateSeqs(aln_dir)

    // Execute Modeltest-NG process
    model = findBestModel(concat_seqs)

    // Execute RAxML-NG processes
    best_tree = bestTree(model, concat_seqs)
    bs_trees = bootstrapTrees(model, concat_seqs)
    support_tree = supportTree(best_tree, bs_trees)

    // Builds phylogenetic tree using Toytree
    phylo_tree = buildTree(support_tree)
}