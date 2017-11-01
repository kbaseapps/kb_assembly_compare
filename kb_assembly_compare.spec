/*
** A KBase module: kb_assembly_compare
**
** This module contains Apps for comparing, combining, and benchmarking assemblies
*/

module kb_assembly_compare {

    /* 
    ** The workspace object refs are of form:
    **
    **    objects = ws.get_objects([{'ref': params['workspace_id']+'/'+params['obj_name']}])
    **
    ** "ref" means the entire name combining the workspace id and the object name
    ** "id" is a numerical identifier of the workspace or object, and should just be used for workspace
    ** "name" is a string identifier of a workspace or object.  This is received from Narrative.
    */
    typedef string workspace_name;
    typedef string sequence;
    typedef string data_obj_name;
    typedef string data_obj_ref;
    typedef int    bool;


    /* filter_contigs_by_length()
    **
    **  Remove Contigs that are under a minimum threshold
    */
    typedef structure {
        workspace_name workspace_name;
	data_obj_ref   input_assembly_refs;   /* Assemblies or AssemblySets */
        data_obj_name  output_name;
    } Filter_Contigs_by_Length_Params;

    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
    } Filter_Contigs_by_Length_Output;

    funcdef run_filter_contigs_by_length (Filter_Contigs_by_Length_Params params)  returns (Filter_Contigs_by_Length_Output) authentication required;


    /* contig_distribution_compare()
    **
    **  Compare Assembly Contig Length Distributions
    */
    typedef structure {
        workspace_name workspace_name;
	data_obj_ref   input_assembly_refs;   /* Assemblies or AssemblySets */
        /*data_obj_name  output_name;*/
    } Contig_Distribution_Compare_Params;

    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
    } Contig_Distribution_Compare_Output;

    funcdef run_contig_distribution_compare (Contig_Distribution_Compare_Params params)  returns (Contig_Distribution_Compare_Output) authentication required;


    /* benchmark_assemblies_against_genomes_with_MUMmer4()
    **
    **  Align benchmark genomes to assembly contigs
    */
    typedef structure {
        workspace_name workspace_name;
	data_obj_ref   input_genome_refs;     /* Genomes or GenomeSets */
	data_obj_ref   input_assembly_refs;   /* Assemblies or AssemblySets */
        /*data_obj_name  output_name;*/
	string         desc;
    } Benchmark_assemblies_against_genomes_with_MUMmer4_Params;

    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
    } Benchmark_assemblies_against_genomes_with_MUMmer4_Output;

    funcdef run_benchmark_assemblies_against_genomes_with_MUMmer4 (Benchmark_assemblies_against_genomes_with_MUMmer4_Params params)  returns (Benchmark_assemblies_against_genomes_with_MUMmer4_Output) authentication required;

};

