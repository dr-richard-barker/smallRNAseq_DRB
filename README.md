# smallRNAseq plan A, B and C

As we explore different small RNAseq analysis pipelines we start by breaking down the NF_Core small RNAseq pipeline to help us understand how it's working. 
Writing up these early stages are good practice while developing bioinformatics processing pipelines even if they don't make it into the final code.

Initially, the NF_core_smallRNA_pipeline was assessed and then broken down into its different components, (For more information see [smallRNAseq_analysis_pipeline_v1_NF_smRNAseq document](https://github.com/dr-richard-barker/smrRNAseq/blob/main/smallRNAseq_analysis_pipeline_v1_NF_smRNAseq.png))

The document (NF_core_smallRNA_pipeline) holds a summary of the NF_cor_small RNAseq pipeline.
The example commands were created to help test its effect on accuracy using synthetic small RNAseq data from SRA. 

The second document contains instructions for installing and running parts of the NF_smRNAseq workflow.
[link to OSDR and SLURM code integration plan](https://github.com/dr-richard-barker/smrRNAseq/blob/main/smallRNAseq_nf_cor_slurm_v2_for_OSDR)

After testing these tools i've come to believe this is pipeline  SRA-> mirnaQC-> sRNAbench-> sRNAde-> mirNOVO is optimal for cross-species analysis.
However, more work needs to be done such as making a snakemate wrapper to help parse data through the pipeline to make high though put. 
([pipeline smallRNA_DREAM_Pipeline_v3.md](https://github.com/dr-richard-barker/smrRNAseq/blob/main/smallRNA_DREAM_Pipeline_v3.md)) 

See the Repository Links descriptions below for more information. 

**The image below outlines the initial pipeline used to align and quantify microRNA**

![deconstructing the nf_core_smRNAseq pipeline](/smallRNAseq_analysis_pipeline_v1_NF_smRNAseq.png)


Developed and maintained by:
Richard John Barker
